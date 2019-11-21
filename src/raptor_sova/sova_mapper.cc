/*
 * sova_mapper.cc
 *
 *  Created on: Nov 13, 2019
 *      Author: Ivan Sovic
 */

#include <raptor_sova/sova_mapper.h>
#include <raptor/mapper_tools.h>

#include <limits>
#include <iostream>
#include <functional>
#include <unordered_map>
#include <tuple>
#include <thread>
#include <cmath>
#include <iomanip>

#include <log/log_tools.h>
#include <algorithm/lis_is.hpp>
#include <writer/output_formatter.h>
#include <index/kxsort.h>
#include <raptor/dp_chain.h>
#include <raptor/interval_tree_builder.h>
#include <debug_tools/write_seed_hit_1.h>
#include <utility/tictoc.h>

#include <raptor_sova/seed_hit_diag.hpp>
#include <raptor_sova/ses_distance_banded.h>

// #define USE_LIS_FILTER

namespace raptor {
namespace sova {

std::unique_ptr<raptor::sova::SovaMapper> createSovaMapper(const mindex::IndexPtr index,
                                         const std::shared_ptr<raptor::sova::ParamsSovaMapper> params) {
    return std::unique_ptr<raptor::sova::SovaMapper>(new raptor::sova::SovaMapper(index, params));
}

SovaMapper::SovaMapper(const mindex::IndexPtr index, const std::shared_ptr<raptor::sova::ParamsSovaMapper> params)
    : index_(index), params_(params) {}

SovaMapper::~SovaMapper() {}

// template <class T>
// std::vector<std::shared_ptr<raptor::TargetHits<T>>> raptor::sova::SovaMapper::FilterForOverlapping_(
//     const std::vector<std::shared_ptr<raptor::TargetHits<T>>>& target_hits, const mindex::SequencePtr& qseq) {
//     std::vector<std::shared_ptr<raptor::TargetHits<T>>> ret;

//     for (size_t i = 0; i < target_hits.size(); i++) {
//         // Hashes of strings could have collisions, depending on the STL implementation.
//         // auto target_val = index_->seqs()->GetSeqByID(target_hits[i]->env()->t_id)->header_hash();
//         // auto query_val = qseq->header_hash();
//         auto target_internal_id = target_hits[i]->env()->t_id;
//         auto target_val = index_->seqs()->GetSeqByID(target_internal_id)->abs_id();
//         auto query_val = qseq->abs_id();
//         if ((params_->overlap_skip_self_hits && target_val == query_val) ||
//             (params_->overlap_single_arc && target_val > query_val)) {
//             continue;
//         }
//         ret.emplace_back(target_hits[i]);
//     }

//     return ret;
// }

/*
 * Forms anchors simply by binning seeds to narrow diagonals.
 * sorted_dh = sorted diagonal hits
*/
std::vector<std::shared_ptr<raptor::TargetAnchorType>> FormDiagonalAnchors(
            const std::vector<mindex128_t>& sorted_dh,
            const mindex::SequencePtr& qseq,
            const mindex::IndexPtr& index,
            int32_t chain_bandwidth, double align_bandwidth, double align_max_diff,
            int32_t min_num_seeds, int32_t min_span,
            bool overlap_skip_self_hits, bool overlap_single_arc) {
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> target_anchors;

    if (sorted_dh.empty()) {
        return target_anchors;
    }

    auto internal_dh = sorted_dh;

    // Just a lookup to identify where to place an anchor.
    std::unordered_map<int32_t, int32_t> tid_map;
    int32_t num_anchors = 0;

    int32_t begin_id = 0;
    auto prev_shp = mindex::SeedHitDiagPacked(internal_dh[begin_id]);
    int32_t begin_diag = prev_shp.TargetPos() - prev_shp.QueryPos();
    // int32_t begin_tid = prev_shp.TargetId();
    // int32_t begin_trev = prev_shp.TargetRev();
    // int32

    for (int32_t i = 0; i < static_cast<int32_t>(internal_dh.size()); ++i) {
        auto curr_shp = mindex::SeedHitDiagPacked(internal_dh[i]);
        int32_t curr_diag = curr_shp.TargetPos() - curr_shp.QueryPos();
        int32_t diag_diff = abs(curr_diag - begin_diag);

        if (curr_shp.TargetId() != prev_shp.TargetId() ||
                curr_shp.TargetRev() != prev_shp.TargetRev() ||
                diag_diff > chain_bandwidth) {

            // Sort the seeds on the same diagonal bandwidth.
            // Since the diagonal has been set to the same value for all seeds
            // in this stretch, they will be sorted by the TargetPos and then QueryPos.
            kx::radix_sort(internal_dh.begin() + begin_id, internal_dh.begin() + i);
            auto begin_shp = mindex::SeedHitDiagPacked(internal_dh[begin_id]);
            auto end_shp = mindex::SeedHitDiagPacked(internal_dh[i-1]);
            int32_t tid = begin_shp.TargetId();
            int32_t num_seeds = i - begin_id;
            int32_t qspan = end_shp.QueryPos() - begin_shp.QueryPos();
            int32_t tspan = end_shp.TargetPos() - begin_shp.TargetPos();
            begin_id = i;
            begin_diag = curr_diag;

            if ((overlap_skip_self_hits && tid == qseq->abs_id()) ||
                (overlap_single_arc && tid > qseq->abs_id())) {
                continue;
            }

            // if (num_seeds < min_num_seeds || qspan < min_span || tspan < min_span) {
            //     continue;
            // }

            if (num_seeds >= min_num_seeds && qspan > min_span && tspan > min_span &&
                    (overlap_skip_self_hits == false || (overlap_skip_self_hits && tid != qseq->abs_id())) &&
                    (overlap_single_arc == false || (overlap_single_arc && tid < qseq->abs_id()))) {
                /////////////////////////
                /// Add a new anchor. ///
                /////////////////////////
                // The t_id of a chain could have clashes with rev cmp. Encoding
                // it again with the rev cmp flag will avoid this.
                int32_t tkey = prev_shp.TargetIdRev() << 1;
                auto it = tid_map.find(tkey);
                if (it == tid_map.end()) {
                    tid_map[tkey] = target_anchors.size();
                    auto new_env = raptor::createMappingEnv(
                            tid, 0, index->len(tid), begin_shp.TargetRev(),
                            qseq->abs_id(), qseq->len(), false);
                    auto anchor = std::shared_ptr<raptor::TargetAnchorType>(
                        new raptor::TargetAnchorType(new_env));
                    target_anchors.emplace_back(anchor);
                    it = tid_map.find(tkey);
                }
                auto& tanchors_ref = target_anchors[it->second];
                int32_t cov_bases_q = 0, cov_bases_t = 0;
                int32_t edit_dist = -1, score = 0;

                int32_t tstart = tanchors_ref->env()->t_rev ? (tanchors_ref->env()->t_len - end_shp.TargetPos()) : begin_shp.TargetPos();
                int32_t tend = tanchors_ref->env()->t_rev ? (tanchors_ref->env()->t_len - begin_shp.TargetPos()) : end_shp.TargetPos();

                auto tseq = index->FetchSeqAsString(tanchors_ref->env()->t_id, tstart, tend, tanchors_ref->env()->t_rev);
                auto qseq_str = qseq->GetSubsequenceAsString(begin_shp.QueryPos(), end_shp.QueryPos());
                int32_t qspan = end_shp.QueryPos() - begin_shp.QueryPos();
                int32_t tspan = end_shp.TargetPos() - begin_shp.TargetPos();
                // if (tanchors_ref->env()->t_rev) {
                //     std::cerr << qseq_str << "\n" << tseq << "\n";
                // }
                int32_t num_diffs = -1;
                num_diffs = raptor::ses::BandedSESDistance(
                            qseq_str.c_str(),
                            qspan,
                            tseq.c_str(),
                            tspan,
                            align_max_diff,
                            align_bandwidth);
                edit_dist = num_diffs;
                score = num_seeds;

                tanchors_ref->hits().emplace_back(
                    raptor::createRegionMapped(
                                num_anchors, -1, tanchors_ref->env(),
                                begin_shp.QueryPos(), end_shp.QueryPos(),
                                begin_shp.TargetPos(), end_shp.TargetPos(),
                                cov_bases_q, cov_bases_t, num_seeds,
                                edit_dist, score, -1, -1, -1, -1));

                ++num_anchors;
            }
        }

        // Set the same diagonal value to any seed within the bandwidth.
        // This will enable sorting by TargetPos.
        mindex::SeedHitDiagPacked::EncodeDiagonal(internal_dh[i], begin_diag);
        prev_shp = curr_shp;
    }

    return target_anchors;
}

std::vector<std::shared_ptr<raptor::TargetAnchorType>> AlignAnchors(
            const std::vector<std::shared_ptr<raptor::TargetAnchorType>> ta,
            const mindex::SequencePtr& qseq,
            const mindex::IndexPtr& index,
            int32_t bandwidth, int32_t min_num_seeds, int32_t min_span,
            bool overlap_skip_self_hits, bool overlap_single_arc) {

    std::vector<std::shared_ptr<raptor::TargetAnchorType>> aligned_ta;
    return aligned_ta;
}

// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_DIAG_0) << "\n";
// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_DIAG_1) << "\n";
// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_DIAG_2) << "\n";
// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_DIAG_3) << "\n";
// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_TARGET_ID) << "\n";
// std::cerr << PrintInt128(mindex::MASK_INV_SEED_HIT_TARGET_REV) << "\n";
std::string PrintInt128(mindex128_t val) {
    std::ostringstream oss;
    std::vector<std::string> vals;
    for (int32_t i = 0; i < 4; ++i) {
        int32_t low_int = (val >> (i * 32)) & mindex::MASK_SEED_HIT_DIAG_0;
        char buf[512];
        sprintf(buf, "%08X", low_int);
        vals.emplace_back(std::string(buf));
    }
    std::reverse(vals.begin(), vals.end());
    for (size_t i = 0; i < vals.size(); ++i) {
        if (i > 0) {
            oss << " ";
        }
        oss << vals[i];
    }
    return oss.str();
}



/*
 * The Mapping workflow. Collects all hits from the index, filters
 * colinear hits in a narrow band of diagonals, forms anchors from the
 * colinear hits (an anchor is a start/end position for each set of colinear
 * hits), and chains the anchors using dynamic programming. It outputs a
 * vector of TargetHits, where each TargetHits contains all chained anchors
 * for a particular target. There can be more than one TargetHits object
 * per target ID, in case mappings lie on different diagonals.
 */
std::shared_ptr<raptor::LinearMappingResult> raptor::sova::SovaMapper::Map(const mindex::SequencePtr& qseq) {
    // std::thread::id thread_id = std::this_thread::get_id();

    TicToc tt_total;
    tt_total.start();

    std::shared_ptr<raptor::LinearMappingResult> result = raptor::createMappingResult(
        qseq->abs_id(), qseq->data().size(), qseq->header(), index_);

    if (result->qseq_len() == 0) {
        result->return_value(raptor::MapperReturnValueBase::QlenIsZero);
        tt_total.stop();
        result->timings()["map_total"] = tt_total.get_microsecs();
        return result;
    }

    if (result->qseq_len() < params_->min_qlen) {
        result->return_value(raptor::MapperReturnValueBase::QseqTooShort);
        tt_total.stop();
        result->timings()["map_total"] = tt_total.get_microsecs();
        return result;
    }

    TicToc tt_collect_hits;
    tt_collect_hits.start();
    std::vector<mindex::SeedHitPacked> seed_hits = index_->CollectHits(
        (const int8_t*) &(qseq->data()[0]), qseq->data().size(), qseq->abs_id());
    tt_collect_hits.stop();

    TicToc tt_convert;
    auto seed_hits_diag = ConvertToSeedHitDiagPacked128Vector(seed_hits);
    tt_convert.stop();

    TicToc tt_sort;
    kx::radix_sort(seed_hits_diag.begin(), seed_hits_diag.end());
    tt_sort.stop();

    TicToc tt_chain;
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> anchors =
                        FormDiagonalAnchors(seed_hits_diag, qseq, index_,
                                params_->chain_bandwidth, params_->align_bandwidth, params_->align_max_diff,
                                params_->min_num_seeds, params_->chain_min_span,
                                params_->ref_and_reads_path_same && params_->overlap_skip_self_hits,
                                params_->ref_and_reads_path_same && params_->overlap_single_arc);
    tt_chain.stop();

    result->timings()["map_hits_collect"] = tt_collect_hits.get_microsecs();
    result->timings()["map_hits_sort"] = tt_sort.get_microsecs();
    result->timings()["map_hits_convert"] = tt_convert.get_microsecs();
    result->timings()["map_hits_chain"] = tt_chain.get_microsecs();

    // Store the results.
    result->target_anchors(anchors);
    result->return_value(raptor::MapperReturnValueBase::OK);

    TicToc tt_debug_out;
    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
            LOG_ALL("Showing debug info for read %d: %s\n", qseq->id(), qseq->header().c_str());

            raptor::WriteSeedHits("temp/debug/mapping-0-seed_hits.csv", seed_hits, index_->params()->k,
                        qseq->header(), qseq->data().size(), std::string("ref"), 0);
            LOG_ALL("Collected seed hits: %ld\n", seed_hits.size());
            LOG_ALL("Anchors grouped on %ld targets:\n", anchors.size());
            for (size_t i = 0; i < anchors.size(); i++) {
                LOG_ALL("   [anchor i = %ld] %s\n", i, anchors[i]->VerbosePointers().c_str());
            }
            LOG_NEWLINE;
        }
    #endif
    tt_debug_out.stop();
    tt_total.stop();

    result->timings()["map_debug_out"] = tt_debug_out.get_microsecs();
    result->timings()["map_total"] = tt_total.get_microsecs();

    return result;
}

}  // namespace sova
}  // namespace raptor
