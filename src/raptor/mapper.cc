/*
 * mapper.cc
 *
 *  Created on: Sep 2, 2017
 *      Author: Ivan Sovic
 */

#include <raptor/mapper.h>
#include <raptor/mapper_tools.h>

#include <limits>
#include <iostream>
#include <functional>
#include <unordered_map>
#include <tuple>
#include <thread>
#include <cmath>

#include <log/log_tools.h>
#include <algorithm/lis_is.hpp>
#include <writer/output_formatter.h>
#include <index/kxsort.h>
#include <raptor/dp_chain.h>
#include <raptor/interval_tree_builder.h>
#include <debug_tools/write_seed_hit_1.h>
#include <utility/tictoc.h>

// #define USE_LIS_FILTER

namespace raptor {

std::unique_ptr<raptor::Mapper> createMapper(const mindex::IndexPtr index,
                                         const std::shared_ptr<raptor::ParamsMapper> params) {
    return std::unique_ptr<raptor::Mapper>(new raptor::Mapper(index, params));
}

Mapper::Mapper(const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsMapper> params)
    : index_(index), params_(params) {}

Mapper::~Mapper() {}

template <class T>
std::vector<std::shared_ptr<raptor::TargetHits<T>>> raptor::Mapper::FilterForOverlapping_(
    const std::vector<std::shared_ptr<raptor::TargetHits<T>>>& target_hits, const mindex::SequencePtr& qseq) {
    std::vector<std::shared_ptr<raptor::TargetHits<T>>> ret;

    for (size_t i = 0; i < target_hits.size(); i++) {
        // Hashes of strings could have collisions, depending on the STL implementation.
        // auto target_val = index_->seqs()->GetSeqByID(target_hits[i]->env()->t_id)->header_hash();
        // auto query_val = qseq->header_hash();
        auto target_internal_id = target_hits[i]->env()->t_id;
        auto target_val = index_->seqs()->GetSeqByID(target_internal_id)->abs_id();
        auto query_val = qseq->abs_id();
        if ((params_->overlap_skip_self_hits && target_val == query_val) ||
            (params_->overlap_single_arc && target_val > query_val)) {
            continue;
        }
        ret.emplace_back(target_hits[i]);
    }

    return ret;
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
std::shared_ptr<raptor::LinearMappingResult> raptor::Mapper::Map(const mindex::SequencePtr& qseq) {
    // std::thread::id thread_id = std::this_thread::get_id();

    TicToc tt_total;
    tt_total.start();

    std::shared_ptr<raptor::LinearMappingResult> result = raptor::createMappingResult(
        qseq->abs_id(), qseq->data().size(), qseq->header(), index_);

    if (result->qseq_len() == 0) {
        result->return_value(raptor::MapperReturnValueBase::QlenIsZero);
        tt_total.stop();
        result->timings()["map_total"] = tt_total.get_msecs();
        return result;
    }

    if (result->qseq_len() < params_->min_qlen) {
        result->return_value(raptor::MapperReturnValueBase::QseqTooShort);
        tt_total.stop();
        result->timings()["map_total"] = tt_total.get_msecs();
        return result;
    }

    TicToc tt_collect_hits;
    tt_collect_hits.start();
    std::vector<mindex::SeedHitPacked> seed_hits = index_->CollectHits(
        (const int8_t*) &(qseq->data()[0]), qseq->data().size(), qseq->abs_id());
    tt_collect_hits.stop();

    // kx::radix_sort(seed_hits.begin(), seed_hits.end());

    TicToc tt_sort, tt_chain;

#ifndef USE_LIS_FILTER
    /////////////////////////////
    /// This is good for now. ///
    /////////////////////////////
    /*
    Current results with default options:
        $ bash run-simulated-celegans-default.sh 98
        sim-tests$ export PATH=$(pwd)/sim-tests/minimap2/:$PATH
        sim-tests$ minimap2/misc/paftools.js mapeval -Q60 out-celegans/raptor-98-aligned-no_HP_k15_w5-default-noalign.paf

        $ paftools.js mapeval -Q60 out-celegans/raptor-92-aligned-no_HP_k15_w5-default-noalign.paf
        Q	60	63676	0	0.000000000	63676
        Q	3	869	245	0.003795801	64545
        Q	2	103	31	0.004269274	64648
        Q	1	188	54	0.005089765	64836
        Q	0	79	24	0.005453285	64915
    */

    tt_sort.start();
    std::sort(seed_hits.begin(), seed_hits.end());  // Sort by TargetId() and TargetPos() first.
    tt_sort.stop();

    tt_chain.start();
    auto filtered_hits = raptor::ChainHits(
                                seed_hits,
                                index_,
                                qseq->abs_id(), qseq->len(),
                                params_->chain_max_skip,
                                params_->chain_max_predecessors,
                                params_->seed_join_dist,
                                params_->diag_margin,
                                params_->min_num_seeds,
                                params_->min_cov_bases,
                                params_->min_dp_score,
                                index_->params()->k);
    tt_chain.stop();
#else
    ////////////////////////////////////
    /// Testing with LIS + dp_chain. ///
    ////////////////////////////////////
    TicToc tt_filter_hits;
    tt_filter_hits.start();
    auto filtered_hits_1 = LISFilterAndGroupByTarget_(
        seed_hits, qseq->abs_id(), qseq->data().size(),
        params_->diag_margin, params_->min_cov_bases, params_->min_num_seeds);
    std::vector<mindex::SeedHitPacked> new_seed_hits;
    for (auto& th : filtered_hits_1) {
        new_seed_hits.insert(new_seed_hits.end(), th->hits().begin(), th->hits().end());
    }
    tt_filter_hits.stop();

    tt_sort.start();
    std::sort(new_seed_hits.begin(), new_seed_hits.end());  // Sort by TargetId() and TargetPos() first.
    tt_sort.stop();

    tt_chain.start();
    auto filtered_hits = raptor::ChainHits(
                                new_seed_hits,
                                index_,
                                qseq->abs_id(), qseq->len(),
                                params_->chain_max_skip,
                                params_->chain_max_predecessors,
                                params_->seed_join_dist,
                                params_->diag_margin,
                                params_->min_num_seeds,
                                params_->min_cov_bases,
                                params_->min_dp_score,
                                index_->params()->k);
    tt_chain.stop();
    ///////////////////////////
#endif

    TicToc tt_filter_overlap;
    tt_filter_overlap.start();
    if (params_->ref_and_reads_path_same && (params_->overlap_skip_self_hits || params_->overlap_single_arc)) {
        filtered_hits = FilterForOverlapping_<mindex::SeedHitPacked>(
            filtered_hits, qseq);
    }
    tt_filter_overlap.stop();

    TicToc tt_make_anchors;
    tt_make_anchors.start();
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> anchors = raptor::mapper::MakeAnchors(filtered_hits);
    tt_make_anchors.stop();

    TicToc tt_sort_anchors;
    tt_sort_anchors.start();
    for (auto& target_anchors : anchors) {
        std::sort(target_anchors->hits().begin(), target_anchors->hits().end(),
                    [](const raptor::AnchorPtr& a, const raptor::AnchorPtr& b) {
                            return ((a->TargetStart() < b->TargetStart()) ||
                                    (a->TargetStart() == b->TargetStart() &&
                                            a->TargetEnd() < b->TargetEnd())); } );
    }
    tt_sort_anchors.stop();

    TicToc tt_score_anchors;
    tt_score_anchors.start();
    // If requested, each anchor will be aligned and the
    // number of match bases used as the score.
    if (params_->score_anchors == true) {
        for (auto& target_anchors : anchors) {
            for (auto& anchor : target_anchors->hits()) {
                int32_t matches = raptor::mapper::CalcMatchRate(qseq, index_, anchor);
                anchor->score(matches);
            }
        }
    }
    tt_score_anchors.stop();

    // Store the results.
    result->target_anchors(anchors);
    result->target_hits(filtered_hits);
    result->return_value(raptor::MapperReturnValueBase::OK);

    TicToc tt_debug_out;
    tt_debug_out.start();

#ifdef RAPTOR_TESTING_MODE
    if (params_->debug_qid == qseq->abs_id() ||
        params_->debug_qname == std::string(qseq->header())) {
        LOG_ALL("Showing debug info for read %d: %s\n", qseq->id(), qseq->header().c_str());

        raptor::WriteSeedHits("temp/debug/mapping-0-seed_hits.csv", seed_hits, index_->params()->k,
                      qseq->header(), qseq->data().size(), std::string("ref"), 0);
        // raptor::WriteSeedHits("temp/debug/mapping-1-filtered_hits_by_target.csv", new_seed_hits, index_->params()->k,
        //               qseq.get_header(), qseq.get_sequence_length(), std::string("ref"), 0);
        // raptor::WriteTargetHits("temp/debug/mapping-1-filtered_hits_by_target.csv", filtered_hits_1,
        //                 index_->params()->k, qseq.get_header(), qseq.get_sequence_length(),
        //                 std::string("ref"), 0);
        raptor::WriteTargetHits("temp/debug/mapping-2-chained_hits.csv", filtered_hits, index_->params()->k,
                        qseq->header(), qseq->data().size(), std::string("ref"), 0);

        LOG_ALL("Collected seed hits: %ld\n", seed_hits.size());
        LOG_ALL("Grouped and filtered hits before anchoring: %ld\n", filtered_hits.size());
        for (size_t i = 0; i < filtered_hits.size(); i++) {
            LOG_ALL("   [i = %ld] filtered_hits[i].hits().size() = %ld\n", i,
                    filtered_hits[i]->hits().size());
        }
        LOG_NEWLINE;
        LOG_ALL("Anchors grouped on %ld targets:\n", anchors.size());
        for (size_t i = 0; i < anchors.size(); i++) {
            LOG_ALL("   [anchor i = %ld] %s\n", i, anchors[i]->VerbosePointers().c_str());
        }
        LOG_NEWLINE;
    }
#endif
    tt_debug_out.stop();

    tt_total.stop();
    result->timings()["map_total"] = tt_total.get_msecs();
    result->timings()["map_hits_collect"] = tt_collect_hits.get_msecs();
    result->timings()["map_hits_sort"] = tt_sort.get_msecs();
    result->timings()["map_hits_chain"] = tt_chain.get_msecs();
    result->timings()["map_hits_filter_overlap"] = tt_filter_overlap.get_msecs();
    result->timings()["map_anchors_make"] = tt_make_anchors.get_msecs();
    result->timings()["map_anchors_sort"] = tt_sort_anchors.get_msecs();
    result->timings()["map_anchors_score"] = tt_score_anchors.get_msecs();
    result->timings()["map_debug_out"] = tt_debug_out.get_msecs();

#ifdef USE_LIS_FILTER
    result->timings()["map_hits_filter_lis"] = tt_filter_hits.get_msecs();
#endif

    return result;
}



#ifdef USE_LIS_FILTER
// #define DETAILED_FILTER_DEBUG_
std::vector<raptor::ChainPtr>
raptor::Mapper::LISFilterAndGroupByTarget_(std::vector<mindex::SeedHitPacked>& seed_hits, indid_t q_id,
                                    ind_t q_len, ind_t diag_margin, ind_t min_cov_bases,
                                    int32_t min_num_hits) {
    auto& index = index_;
    int32_t seed_len = index->params()->k;

    std::sort(seed_hits.begin(), seed_hits.end());

    // Each vector element corresponds to one target ID.
    std::vector<raptor::ChainPtr> ret_target_hits;

    // Relation from target ID to the ordinal number of the vector element where the target is.
    std::unordered_map<indid_t, size_t> t_id_to_ret;

    if (seed_hits.size() == 0) {
        return ret_target_hits;
    }

    size_t n_hits = seed_hits.size();

    // Comparison function to sort the seed hits for LIS.
    std::function<bool(const mindex::SeedHitPacked& a, const mindex::SeedHitPacked& b)>
        comp_packed_sort = [](const mindex::SeedHitPacked& a, const mindex::SeedHitPacked& b) {
            auto a_qpos = a.QueryPos();
            auto a_tpos = a.TargetPos();
            auto b_qpos = b.QueryPos();
            auto b_tpos = b.TargetPos();
            return (a_qpos < b_qpos || (a_qpos == b_qpos && a_tpos < b_tpos));
        };
    std::function<bool(const mindex::SeedHitPacked& a, const mindex::SeedHitPacked& b)>
        comp_packed_lis = [](const mindex::SeedHitPacked& a, const mindex::SeedHitPacked& b) {
            auto a_qpos = a.QueryPos();
            auto a_tpos = a.TargetPos();
            auto b_qpos = b.QueryPos();
            auto b_tpos = b.TargetPos();
            return (a_qpos < b_qpos && a_tpos < b_tpos);
        };

    // Keep track of the first seed in the LIS streak.
    int32_t first = 0;

    // Values of the previous seed in the streak, needed to break
    // the elongation of the streak.
    indid_t prev_t_id = seed_hits[0].TargetId();
    bool prev_t_rev = seed_hits[0].TargetRev();
    // ind_t prev_diag = seed_hits[0].Diag();
    ind_t prev_q_pos = seed_hits[0].QueryPos();
    ind_t prev_t_pos = seed_hits[0].TargetPos();
    ind_t prev_diag = prev_t_pos - prev_q_pos;

    for (size_t i = 1; i <= n_hits; i++) {
        // Set the default values for the current seed hit,
        // in case i == n_hits this will ensure everything works fine.
        indid_t curr_t_id = prev_t_id;
        indid_t curr_t_rev = prev_t_rev;
        ind_t curr_diag = prev_diag;
        ind_t curr_q_pos = prev_q_pos;
        ind_t curr_t_pos = prev_t_pos;

        // Get the actual values.
        if (i < n_hits) {
            curr_t_id = seed_hits[i].TargetId();
            // curr_diag = seed_hits[i].Diag();
            curr_q_pos = seed_hits[i].QueryPos();
            curr_t_pos = seed_hits[i].TargetPos();
            curr_diag = curr_t_pos - curr_q_pos;
        }

        // Calculates the distance between two seeds.
        ind_t seed_dist = std::max(abs(curr_q_pos - prev_q_pos), abs(curr_t_pos - prev_t_pos));

#ifdef DETAILED_FILTER_DEBUG_
            fprintf (stderr, "[i = %llu, n_hits = %ld] first = %d, curr_q_pos = %d, curr_t_pos = %d, curr_diag = %d, curr_t_id = %ld, seed_dist = %d\n",
                i, n_hits, first, curr_q_pos, curr_t_pos, curr_diag, curr_t_id, seed_dist);
        // if (curr_q_pos == 49 || prev_q_pos == 49) {
        //     fprintf (stderr, "[curr_q_pos = %d, first = %d] i = %llu, n_hits = %ld, curr_q_pos = %d, curr_t_pos = %d, curr_diag = %d, curr_t_id = %ld, seed_dist = %d\n",
        //         first, curr_q_pos, i, n_hits, curr_q_pos, curr_t_pos, curr_diag, curr_t_id, seed_dist);
        //     fprintf (stderr, "[prev_q_pos = %d] i = %llu, n_hits = %ld, prev_q_pos = %d, prev_t_pos = %d, prev_diag = %d, prev_t_id = %ld\n",
        //         prev_q_pos, i, n_hits, prev_q_pos, prev_t_pos, prev_diag, prev_t_id);
        //     fflush(stderr);
        // }
#endif

        // Break the streak if required and calculate LIS.
        if (i == n_hits || curr_t_id != prev_t_id || curr_t_rev != prev_t_rev || abs(curr_diag - prev_diag) >= diag_margin ||
            seed_dist > params_->seed_join_dist) {
            int32_t last = i;

#ifdef DETAILED_FILTER_DEBUG_
            fprintf (stderr, "    -> Rounding: first = %d, last = %d\n", first, last);
            if (i == n_hits) { fprintf (stderr, "  i == n_hits\n"); }
            if (curr_t_id != prev_t_id) { fprintf (stderr, "  curr_t_id != prev_t_id\n"); }
            if (curr_t_rev != prev_t_rev) { fprintf (stderr, "  curr_t_rev != prev_t_rev\n"); }
            if (abs(curr_diag - prev_diag) >= diag_margin ) { fprintf (stderr, "  abs(curr_diag - prev_diag) >= diag_margin\n"); }
            if (seed_dist > params_->seed_join_dist) { fprintf (stderr, "  seed_dist > params_->seed_join_dist\n"); }
#endif

            if ((last - first) >= min_num_hits) {
#ifdef DETAILED_FILTER_DEBUG_
                fprintf (stderr, "        -> Running LIS: first = %d, last = %d\n", first, last);
#endif

                // Calculate the LIS of the hits.
                std::sort(seed_hits.begin() + first, seed_hits.begin() + last, comp_packed_sort);
                auto lis_hits = raptor::LIS(seed_hits, first, last, comp_packed_lis);

                // This just fetches the target details.
                indid_t index_t_id = prev_t_id;
                bool t_rev = prev_t_rev;
                int32_t index_t_start = index->starts()[index_t_id];
                int32_t t_len = index->lens()[index_t_id];

                int32_t cov_bases_q = 0, cov_bases_t = 0;
                raptor::mapper::CalcHitCoverage(lis_hits, seed_len, 0, lis_hits.size(), cov_bases_q, cov_bases_t);

#ifdef DETAILED_FILTER_DEBUG_
                fprintf (stderr, "        -> After LIS: lis_hits.size() = %llu, cov_bases_q = %d, cov_bases_t = %d\n", lis_hits.size(), cov_bases_q, cov_bases_t);
#endif

                if (lis_hits.size() >= min_num_hits) {
                // if (cov_bases_q >= min_cov_bases && cov_bases_t >= min_cov_bases) {

                    indid_t prev_t_id_rev_packed = (prev_t_id << 1) | ((prev_t_rev) ? 0x01 : 0x00);

                    auto it_loc = t_id_to_ret.find(prev_t_id_rev_packed);

                    // If the target ID does not already exist in the output, add it.
                    if (it_loc == t_id_to_ret.end()) {
                        t_id_to_ret[prev_t_id_rev_packed] = ret_target_hits.size();
                        auto new_env = raptor::createMappingEnv(index_t_id, index_t_start, t_len, t_rev,
                                                            q_id, q_len, false);
                        auto new_target_hits =
                            raptor::ChainPtr(
                                new raptor::TargetHits<mindex::SeedHitPacked>(new_env));
                        ret_target_hits.emplace_back(new_target_hits);
                        it_loc = t_id_to_ret.find(prev_t_id_rev_packed);
                    }

                    // Add the new hits. The cov_bases values are only an estimate from now on,
                    // because neighboring hits from different streaks may overlap.
                    auto& th = ret_target_hits[it_loc->second];
                    th->AppendHits(lis_hits);
                    th->cov_bases_q(th->cov_bases_q() + cov_bases_q);
                    th->cov_bases_t(th->cov_bases_t() + cov_bases_t);

#ifdef DETAILED_FILTER_DEBUG_
                    fprintf (stderr, "  ================ Appended.\n");
#endif

                }
            }

            prev_t_id = curr_t_id;
            prev_t_rev = curr_t_rev;
            prev_diag = curr_diag;
            prev_q_pos = curr_q_pos;
            prev_t_pos = curr_t_pos;
            first = last;
        }
// Do ovdje!

    }

    // exit(1);

    return ret_target_hits;
}
#endif

}  // namespace raptor
