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

#include <raptor_sova/sova_overlap.h>

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

raptor::sova::OverlapPtr MakeOverlap(std::vector<mindex128_t>& hits,
            const mindex::SequencePtr& qseq,
            const mindex::IndexPtr& index,
            int32_t begin_id, int32_t end_id,
            int32_t min_tpos_id, int32_t max_tpos_id) {

    auto begin_shp = mindex::SeedHitDiagPacked(hits[min_tpos_id]);
    auto end_shp = mindex::SeedHitDiagPacked(hits[max_tpos_id]);

    int32_t tid = begin_shp.TargetId();
    int32_t num_seeds = end_id - begin_id;
    int32_t qspan = end_shp.QueryPos() - begin_shp.QueryPos();
    int32_t tspan = end_shp.TargetPos() - begin_shp.TargetPos();

    if (end_shp.TargetId() != tid) {
        return nullptr;
    }

    float score = num_seeds;
    float idt = 0.0;
    int32_t edit_dist = -1;
    return raptor::sova::createOverlap(
        qseq->abs_id(), tid,
        score, idt,
        0, begin_shp.QueryPos(), end_shp.QueryPos(), qseq->len(),
        begin_shp.TargetRev(), begin_shp.TargetPos(), end_shp.TargetPos(), index->len(tid),
        edit_dist, num_seeds
    );
}

/*
 * Forms anchors simply by binning seeds to narrow diagonals.
 * sorted_dh = sorted diagonal hits
*/
std::vector<raptor::sova::OverlapPtr> FormDiagonalAnchors(
            const std::vector<mindex128_t>& sorted_dh,
            const mindex::SequencePtr& qseq,
            const mindex::IndexPtr& index,
            int32_t chain_bandwidth,
            int32_t min_num_seeds, int32_t min_span,
            bool overlap_skip_self_hits, bool overlap_single_arc) {
    std::vector<raptor::sova::OverlapPtr> overlaps;

    if (sorted_dh.empty()) {
        return overlaps;
    }

    auto internal_dh = sorted_dh;

    int32_t begin_id = 0;
    auto prev_shp = mindex::SeedHitDiagPacked(internal_dh[begin_id]);
    int32_t begin_diag = prev_shp.TargetPos() - prev_shp.QueryPos();
    int32_t min_tpos_id = 0;
    int32_t max_tpos_id = 0;
    int32_t num_hits = static_cast<int32_t>(internal_dh.size());
    uint64_t min_pos = internal_dh[begin_id] & mindex::MASK_SEED_HIT_DIAG_LOWER_64;
    uint64_t max_pos = internal_dh[begin_id] & mindex::MASK_SEED_HIT_DIAG_LOWER_64;

    for (int32_t i = 0; i < num_hits; ++i) {
        auto curr_shp = mindex::SeedHitDiagPacked(internal_dh[i]);
        int32_t curr_diag = curr_shp.TargetPos() - curr_shp.QueryPos();
        auto curr_tpos = curr_shp.TargetPos();
        int32_t diag_diff = abs(curr_diag - begin_diag);
        uint64_t tq_pos = internal_dh[i] & mindex::MASK_SEED_HIT_DIAG_LOWER_64;

        if (curr_shp.TargetId() != prev_shp.TargetId() ||
                curr_shp.TargetRev() != prev_shp.TargetRev() ||
                diag_diff > chain_bandwidth) {

            auto ovl = MakeOverlap(internal_dh, qseq, index, begin_id, i, min_tpos_id, max_tpos_id);
            begin_id = i;
            begin_diag = curr_diag;

            // Add a new overlap.
            if (ovl->num_seeds >= min_num_seeds && ovl->ASpan() > min_span && ovl->BSpan() > min_span &&
                    (overlap_skip_self_hits == false || (overlap_skip_self_hits && ovl->b_id != qseq->abs_id())) &&
                    (overlap_single_arc == false || (overlap_single_arc && ovl->b_id < qseq->abs_id()))) {

                overlaps.emplace_back(std::move(ovl));
            }

            min_tpos_id = max_tpos_id = i;
            min_pos = max_pos = tq_pos;
        }
        // Track the minimum and maximum target positions for each diagonal.
        if (tq_pos < min_pos) {
            min_tpos_id = i;
            min_pos = tq_pos;
        }
        if (tq_pos > max_pos) {
            max_tpos_id = i;
            max_pos = tq_pos;
        }
        // Set the same diagonal value to any seed within the bandwidth.
        // This will enable sorting by TargetPos.
        mindex::SeedHitDiagPacked::EncodeDiagonal(internal_dh[i], begin_diag);
        prev_shp = curr_shp;
    }

    if ((num_hits - begin_id) > 0) {
        auto ovl = MakeOverlap(internal_dh, qseq, index, begin_id, num_hits, min_tpos_id, max_tpos_id);
        // Add a new overlap.
        if (ovl->num_seeds >= min_num_seeds && ovl->ASpan() > min_span && ovl->BSpan() > min_span &&
                (overlap_skip_self_hits == false || (overlap_skip_self_hits && ovl->b_id != qseq->abs_id())) &&
                (overlap_single_arc == false || (overlap_single_arc && ovl->b_id < qseq->abs_id()))) {

            overlaps.emplace_back(std::move(ovl));
        }
    }

    return overlaps;
}

/*
 * This function will take an overlap, and perform banded global alignment
 * of the overlap from it's start to end.
 * The produced overlap will have the same start/end coordinates, and only
 * different scores/edit distance.
*/
raptor::sova::OverlapPtr AlignOverlap(
                const mindex::IndexPtr& index,
                const mindex::SequencePtr& qseq,
                const raptor::sova::OverlapPtr& ovl,
                double align_bandwidth, double align_max_diff) {
    if (ovl == nullptr) {
        return nullptr;
    }

    raptor::sova::OverlapPtr ret = raptor::sova::createOverlap(ovl);
    int32_t edit_dist = -1, score = 0;
    auto qseq_str = qseq->GetSequenceAsString();
    auto qseq_str_rev = ReverseComplement(qseq_str);

    ///////////////////////////
    /// Align forward pass. ///
    ///////////////////////////
    {
        int32_t qstart = ovl->a_start;
        int32_t qend = ovl->a_len;
        int32_t qspan = qend - qstart;

        int32_t tstart_fwd = ovl->b_rev ? (ovl->b_len - ovl->b_end) : ovl->b_start;
        int32_t tend_fwd = ovl->b_rev ? (ovl->b_len - ovl->b_start) : ovl->b_end;
        std::string tseq = (ovl->b_rev) ?
                                index->FetchSeqAsString(ovl->b_id, 0, tend_fwd, ovl->b_rev) :
                                index->FetchSeqAsString(ovl->b_id, tstart_fwd, ovl->b_len, ovl->b_rev);
        int32_t tspan = tseq.size();

        auto ses_result = raptor::ses::BandedSESDistanceAdvanced(
                                qseq_str.c_str() + qstart,
                                qspan,
                                tseq.c_str(),
                                tspan,
                                align_max_diff,
                                align_bandwidth, 2, -1);

        // ret->a_end = (ses_result.valid) ? ses_result.last_q : ses_result.max_q;
        // ret->b_end = (ses_result.valid) ? ses_result.last_t : ses_result.max_t;
        ret->a_end = ses_result.last_q;
        ret->b_end = ses_result.last_t;

        ret->a_end += ovl->a_start;
        ret->b_end += ovl->b_start; // I think this is fine for both fwd and rev. For fwd it's definitely good. For rev, the alignment begins at "tend_fwd" (because it's reverse complemented), and the distance from tend_fwd to ovl->b_len is (ovl->b_len - tend_fwd) == ovl->b_start;

        int32_t num_diffs = ses_result.diffs;
        ret->edit_dist = num_diffs;
        ret->score = ret->num_seeds;

        // std::cerr << "[1] " << qseq->header() << "\t"
        //     << "last(" << ses_result.last_q + ovl->a_start << ", " << ses_result.last_t + ovl->b_start << ", " << ses_result.last_score << ")" << "\t"
        //     << "max(" << ses_result.max_q + ovl->a_start << ", " << ses_result.max_t + ovl->b_start << ", " << ses_result.max_score << ", " << ses_result.max_score_diffs << ")" << "\t"
        //     << ovl->a_name << "\t" << ovl->b_name << "\t" << ovl->score << "\t" << ovl->identity << "\t"
        //     << ovl->a_rev << "\t" << ovl->a_start << "\t" << ovl->a_end << "\t" << ovl->a_len << "\t"
        //     << ovl->b_rev << "\t" << ovl->b_start << "\t" << ovl->b_end << "\t" << ovl->b_len << "\t"
        //     << "\n";
    }

    // raptor::sova::SovaMapper::PrintOverlapAsM4(stderr, index, qseq, ovl);
    // raptor::sova::SovaMapper::PrintOverlapAsM4(stderr, index, qseq, ret);

    ///////////////////////////
    /// Align reverse pass. ///
    ///////////////////////////
    {
        // Reverse query coordinates.
        int32_t qstart = ret->a_len - ret->a_end;
        int32_t qend = ret->a_len;
        int32_t qspan = qend - qstart;

        int32_t tstart_fwd = ret->b_rev ? (ret->b_len - ret->b_end) : ret->b_start;
        int32_t tend_fwd = ret->b_rev ? (ret->b_len - ret->b_start) : ret->b_end;
        std::string tseq = (ovl->b_rev) ?
                                index->FetchSeqAsString(ret->b_id, tstart_fwd, ret->b_len, !ret->b_rev) :
                                index->FetchSeqAsString(ret->b_id, 0, tend_fwd, !ret->b_rev);
        int32_t tspan = tseq.size();

        auto ses_result = raptor::ses::BandedSESDistanceAdvanced(
                                qseq_str_rev.c_str() + qstart,
                                qspan,
                                tseq.c_str(),
                                tspan,
                                align_max_diff,
                                align_bandwidth, 2, -1);

        // std::cerr << "    tspan = " << tspan << "\n";
        // std::cerr << "[2] " << qseq->header() << "\t"
        //     << "last(" << ses_result.last_q << ", " << ses_result.last_t << ", " << ses_result.last_score << ", " << ses_result.diffs << ")" << "\t"
        //     << "max(" << ses_result.max_q << ", " << ses_result.max_t << ", " << ses_result.max_score << ", " << ses_result.max_score_diffs << ")" << "\t"
        //     << ovl->a_name << "\t" << ovl->b_name << "\t" << ovl->score << "\t" << ovl->identity << "\t"
        //     << ovl->a_rev << "\t" << ovl->a_start << "\t" << ovl->a_end << "\t" << ovl->a_len << "\t"
        //     << ovl->b_rev << "\t" << ovl->b_start << "\t" << ovl->b_end << "\t" << ovl->b_len << "\t"
        //     << "\n";

        // ret->a_start = ret->a_end - ((ses_result.valid) ? ses_result.last_q : ses_result.max_q);
        // ret->b_start = ret->b_end - ((ses_result.valid) ? ses_result.last_t : ses_result.max_t);
        ret->a_start = ret->a_end - ses_result.last_q;
        ret->b_start = ret->b_end - ses_result.last_t;

        int32_t num_diffs = (ses_result.valid) ? ses_result.diffs : ses_result.max_score_diffs;
        ret->edit_dist = num_diffs;
        ret->score = ret->num_seeds;
    }

    // raptor::sova::SovaMapper::PrintOverlapAsM4(stderr, index, qseq, ret);
    // std::cerr << "\n";

    return ret;
}

/*
 * ExtendAlignOverlap will take an overlap and attempt to align it beyond
 * the start and end positions, "extending" it on both ends.
 * The new produced overlap will have modified start/end coordinates.
*/
raptor::sova::OverlapPtr ExtendAlignOverlap(
                const mindex::IndexPtr& index,
                const mindex::SequencePtr& qseq,
                const raptor::sova::OverlapPtr& ovl,
                double align_bandwidth, double align_max_diff) {
    if (ovl == nullptr) {
        return nullptr;
    }
    raptor::sova::OverlapPtr ret = raptor::sova::createOverlap(ovl);

    int32_t edit_dist = -1, score = 0;

    int32_t tstart = ovl->b_rev ? (ovl->b_len - ovl->b_end) : ovl->b_start;
    int32_t tend = ovl->b_rev ? (ovl->b_len - ovl->b_start) : ovl->b_end;

    auto tseq = index->FetchSeqAsString(ovl->b_id, tstart, tend, ovl->b_rev);
    auto qseq_str = qseq->GetSubsequenceAsString(ovl->a_start, ovl->a_end);
    int32_t qspan = ovl->a_end - ovl->a_start;
    int32_t tspan = ovl->b_end - ovl->b_start;
    int32_t num_diffs = raptor::ses::BandedSESDistance(
                            qseq_str.c_str(),
                            qspan,
                            tseq.c_str(),
                            tspan,
                            align_max_diff,
                            align_bandwidth);
    ret->edit_dist = num_diffs;
    ret->score = ret->num_seeds;

    return ret;
}

void AlignOverlaps(
                const mindex::IndexPtr& index,
                const mindex::SequencePtr& qseq,
                std::vector<raptor::sova::OverlapPtr>& overlaps,
                double align_bandwidth, double align_max_diff) {
    for (size_t i = 0; i < overlaps.size(); ++i) {
        overlaps[i] = AlignOverlap(index, qseq, overlaps[i], align_bandwidth, align_max_diff);
    }
}

std::vector<std::shared_ptr<raptor::TargetAnchorType>> OverlapsToTargetAnchors(const std::vector<raptor::sova::OverlapPtr>& overlaps) {
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> ta;

    if (overlaps.empty()) {
        return ta;
    }

    std::unordered_map<int32_t, int32_t> tid_map;
    for (int32_t i = 0, curr_ta = -1; i < static_cast<int32_t>(overlaps.size()); ++i) {
        const auto& curr = overlaps[i];
        if (i == 0 || overlaps[i]->b_id != overlaps[i-1]->b_id ||
                        overlaps[i]->b_rev != overlaps[i-1]->b_rev) {
            // Create a new target anchor.
            // The t_id of a chain could have clashes with rev cmp. Encoding
            // it again with the rev cmp flag will avoid this.
            // const auto& prev = (i > 0) ? overlaps[i-1] : overlaps[i];
            int32_t tkey = (curr->b_id << 1) | static_cast<int32_t>(curr->b_rev);
            auto it = tid_map.find(tkey);
            if (it == tid_map.end()) {
                tid_map[tkey] = ta.size();
                auto new_env = raptor::createMappingEnv(
                        curr->b_id, 0, curr->b_len, curr->b_rev,
                        curr->a_id, curr->a_len, curr->a_rev);
                auto new_ta_container = std::shared_ptr<raptor::TargetAnchorType>(
                        new raptor::TargetAnchorType(new_env));
                ta.emplace_back(new_ta_container);
                it = tid_map.find(tkey);
                curr_ta = it->second;
            }
        }
        auto& ta_container = ta[curr_ta];

        int32_t cov_bases_q = curr->a_end - curr->a_start;
        int32_t cov_bases_t = curr->b_end - curr->b_start;
        ta_container->hits().emplace_back(
            raptor::createRegionMapped(
                        i, -1, ta_container->env(),
                        curr->a_start, curr->a_end,
                        curr->b_start, curr->b_end,
                        cov_bases_q, cov_bases_t,
                            curr->num_seeds, curr->edit_dist, curr->score,
                        -1, -1, -1, -1));
    }

    return ta;
}

// std::vector<std::shared_ptr<raptor::TargetAnchorType>> AlignAnchors(
//             const std::vector<std::shared_ptr<raptor::TargetAnchorType>> ta,
//             const mindex::SequencePtr& qseq,
//             const mindex::IndexPtr& index,
//             int32_t bandwidth, int32_t min_num_seeds, int32_t min_span,
//             bool overlap_skip_self_hits, bool overlap_single_arc) {

//     std::vector<std::shared_ptr<raptor::TargetAnchorType>> aligned_ta;

//     for (auto& one_ta: ta) {
// 		if (one_ta == nullptr) {
// 			continue;
// 		}

// 		// const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result->regions();


//         // for (size_t i = 0; i < regions_to_write.size(); i++) {
//         //     auto aln = regions_to_write[i];
//         //     bool is_secondary = (aln->PathId() > 0);
//         //     bool is_supplementary = (aln->SegmentId() > 0);
//         //     auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());
// 		// 	PrintAsM4(stdout, index, qseq, aln);
//         // }
//     }

//     return aligned_ta;
// }

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
    auto overlaps = FormDiagonalAnchors(seed_hits_diag, qseq, index_,
                        params_->chain_bandwidth,
                        params_->min_num_seeds, params_->chain_min_span,
                        params_->ref_and_reads_path_same && params_->overlap_skip_self_hits,
                        params_->ref_and_reads_path_same && params_->overlap_single_arc);
    AlignOverlaps(index_, qseq, overlaps, params_->align_bandwidth, params_->align_max_diff);
    auto anchors = OverlapsToTargetAnchors(overlaps);
    // std::vector<std::shared_ptr<raptor::TargetAnchorType>> anchors;
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

void SovaMapper::PrintOverlapAsM4(FILE *fp_out, const mindex::IndexPtr& index, const mindex::SequencePtr& qseq,
                            const raptor::sova::OverlapPtr& ovl) {
    int32_t edit_dist = ovl->edit_dist;

    double qspan = ovl->ASpan();
    double tspan = ovl->BSpan();
	double identity = 0.0;

    if (edit_dist >= 0) {
        double edit_dist_double = edit_dist;
        double identity_q = (qspan != 0) ? ((qspan  - edit_dist_double) / qspan) : -2.0;
        double identity_t = (tspan != 0) ? ((tspan - edit_dist_double) / tspan) : -2.0;
        identity = std::min(identity_q, identity_t);
    }

    int32_t t_start = ovl->b_start;
    int32_t t_end = ovl->b_end;
    int32_t t_is_rev = ovl->b_rev;
    int32_t t_len = ovl->b_len;

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

	fprintf(fp_out, "%s %s %d %.2lf %d %d %d %d %d %d %d %d\n",
			qseq->header().c_str(), index->header(ovl->b_id).c_str(),
            static_cast<int32_t>(ovl->score), 100.0 * identity,
			static_cast<int32_t>(ovl->a_rev), ovl->a_start, ovl->a_end, ovl->a_len,
			static_cast<int32_t>(t_is_rev), t_start, t_end, t_len);

}

}  // namespace sova
}  // namespace raptor
