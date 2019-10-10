/*
 * dp_chain.cc
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#include <raptor/dp_chain.h>
#include <raptor/mapper_tools.h>
#include <cmath>
#include <algorithm>
#include <utility/math.hpp>

namespace raptor {

constexpr int32_t PlusInf = std::numeric_limits<int32_t>::max() - 10000;  // Leave a margin.

// #define DEBUG_DP_VERBOSE_

std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>> ChainHits(
    const mindex::IndexPtr index,       // Optional. If nullptr, then the target sequence length won't be initialized in the MappingEnv.
    int64_t qseq_abs_id,
    int64_t qseq_len,
    int32_t chain_max_skip,
    int32_t chain_max_predecessors,
    int32_t seed_join_dist,
    int32_t diag_margin,
    int32_t min_num_seeds,
    int32_t min_cov_bases, int32_t min_dp_score, int32_t k,
    const std::vector<mindex::SeedHitPacked>& hits) {

    std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>> chains;

    if (hits.size() == 0) {
        return chains;
    }

    if (chain_max_skip <= 0) {
        return chains;
    }

    int32_t n_hits = (int32_t) hits.size();

    // Zeroth element will be the "Null" state.
    std::vector<int32_t> dp(n_hits + 1, 0);         // Initial dp score is 0, for local alignment.
    std::vector<int32_t> pred(n_hits + 1, 0);
    std::vector<int32_t> chain_id(n_hits + 1, -1);  // For each node, it's chain ID is the same as of it's predecessor.
    int32_t num_chains = 0;

    auto hits_unpacked = mindex::UnpackMinimizerHitVector(hits);

    const double lin_factor = 0.01 * k;

    for (int32_t i = 1; i < (n_hits + 1); i++) {
        int32_t x_i_start = hits_unpacked[i - 1].QueryPos();
        int32_t y_i_start = hits_unpacked[i - 1].TargetPos();
        int32_t x_i_end = x_i_start + k;  // TODO: This is not accurate for homopolymer suppression.
        int32_t y_i_end = y_i_start + k;
        int32_t l_i = y_i_start - x_i_start;
        int32_t t_id_i = hits_unpacked[i - 1].TargetId();
        int32_t t_rev_i = hits_unpacked[i - 1].TargetRev();

        // Add the initial gap open penalty.
        int32_t x_i_score = k;
        int32_t new_dp_val = x_i_score;
        int32_t new_dp_pred = 0;
        int32_t new_dp_chain = num_chains;

        int32_t num_skipped_predecessors = 0;

        int32_t num_processed = 0;

        int32_t min_j = 0;
        // int32_t min_j = std::max(0, i - 1000);
        min_j = (chain_max_predecessors <= 0) ? 0 : std::max(static_cast<int32_t>(0), (i - 1 - chain_max_predecessors));
        for (int32_t j = (i - 1); j > min_j; j--) {
            // if (num_processed > chain_max_predecessors) {
            //     // std::cerr << "Bump 1! num_processed = " << num_processed << ", chain_max_predecessors = " << chain_max_predecessors << std::endl;
            //     break;
            // }

#ifdef EXPERIMENTAL_QUERY_MASK
            bool is_tandem = hits_unpacked[j - 1].QueryMask() & MINIMIZER_HIT_TANDEM_FLAG;
#endif

            int32_t x_j_start = hits_unpacked[j - 1].QueryPos();
            int32_t y_j_start = hits_unpacked[j - 1].TargetPos();
            int32_t t_id_j = hits_unpacked[j - 1].TargetId();
            int32_t t_rev_j = hits_unpacked[j - 1].TargetRev();
            // int32_t x_j_end = x_j_start + k;  // TODO: This is not accurate for homopolymer suppression.
            // int32_t y_j_end = y_j_start + k;

            int32_t dist_x = x_i_start - x_j_start;   // If < 0, it's not a predecessor.
            int32_t dist_y = y_i_start - y_j_start;
            // int32_t l_j = y_j_start - x_j_start;

            int32_t gap_dist = (dist_x < dist_y) ? (dist_y - dist_x) : (dist_x - dist_y);

            if (t_id_j != t_id_i || t_rev_j != t_rev_i) {
                break;
            }
            if (dist_y > seed_join_dist) {
                break;
            }

#ifdef EXPERIMENTAL_QUERY_MASK
            if (is_tandem) {
                continue;
            }
#endif

            if (x_i_start <= x_j_start || y_i_start <= y_j_start) {
                continue;
            }
            if (gap_dist > diag_margin) {
                continue;
            }
            if (dist_x > seed_join_dist) {
                continue;
            }

            num_processed += 1;

            int32_t lin_part = (int32_t) (gap_dist * lin_factor);
            int32_t log_part = (int32_t) ((gap_dist == 0) ? 0 : raptor::utility::ilog2_32(gap_dist));
            int32_t edge_score = lin_part + (log_part >> 1);

            int32_t x_j_score = std::min(k, (int32_t) std::min(abs(dist_x), abs(dist_y)));
            int32_t score_ij = dp[j] + x_j_score - edge_score;

            if (score_ij >= new_dp_val) {
                new_dp_pred = j;
                new_dp_val = score_ij;
                new_dp_chain = chain_id[j];

                // This is the main difference to how I previously calculated the scan_depth.
                num_skipped_predecessors -= 1;
                num_skipped_predecessors = std::max(0, num_skipped_predecessors);

            } else {
                num_skipped_predecessors += 1;
                if (num_skipped_predecessors > chain_max_skip) {
                    // std::cerr << "Bump! num_skipped_predecessors = " << num_skipped_predecessors << ", chain_max_skip = " << chain_max_skip << std::endl;
                    break;
                }
            }
        }

        dp[i] = new_dp_val;
        pred[i] = new_dp_pred;
        chain_id[i] = new_dp_chain;
        if (new_dp_chain == num_chains) {
            num_chains += 1;
        }

    }

    // Find the maximum of every chain for backtracking.
    std::vector<int32_t> chain_maxima(num_chains, -PlusInf);
    for (int32_t i = 1; i < (n_hits + 1); i++) {
        if (chain_maxima[chain_id[i]] == -PlusInf || dp[i] >= dp[chain_maxima[chain_id[i]]]) {
            chain_maxima[chain_id[i]] = i;
        }
    }
    // Testing with taking the last seed as the start for backtrack, instead
    // of the maximum scoring one. Mismatch distances play a significant role
    // in the scoring.
    // std::vector<int32_t> chain_last_score(num_chains, -PlusInf);
    // for (int32_t i = 1; i < (n_hits + 1); i++) {
    //     chain_last_score[chain_id[i]] = i;
    // }

    std::shared_ptr<raptor::MappingEnv> new_env = nullptr;

    // Backtrack.
    for (int32_t i = 0; i < (int32_t)chain_maxima.size(); i++) {

        // Trace back from the maxima.
        int32_t node_id = chain_maxima[i];
        int32_t score = dp[node_id];

        if (score < min_dp_score) {
            continue;
        }

        std::vector<int32_t> nodes;
        while (node_id > 0) {
            nodes.emplace_back(node_id - 1);  // The "- 1" is because of the DP offset.
            node_id = pred[node_id];
        }
        // Reverse the backtracked nodes.
        std::reverse(nodes.begin(), nodes.end());
        if (nodes.size() == 0 || nodes.size() < min_num_seeds) {
            continue;
        }

        /////////////////////////
        /// Create the chain. ///
        /////////////////////////

        // Construct the MappingEnv.

        indid_t index_t_id = hits_unpacked[nodes.front()].TargetId();
        bool t_rev = hits_unpacked[nodes.front()].TargetRev();
        if (new_env == nullptr || new_env->t_id != index_t_id || new_env->t_rev != t_rev) {
            int32_t index_t_start = 0; // index->starts()[index_t_id];
            int32_t t_len = (index == nullptr) ? 0 : index->len(index_t_id);
            new_env = raptor::createMappingEnv(index_t_id, index_t_start, t_len, t_rev, qseq_abs_id, qseq_len, false);
        }

        auto chain = std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>(new raptor::TargetHits<mindex::SeedHitPacked>(new_env));

        for (auto& node : nodes) {
            chain->hits().emplace_back(hits[node]);
        }

        // Penalize the distance from the end of the query.
        // Otherwise, shorted chains near the beginning would
        // prevail longer ones in some cases.
        // int32_t chain_dist_to_end = qseq.get_sequence_length() - chain->hits().back().QueryPos();
        // int32_t chain_score = score - chain_dist_to_end * params->chain_penalty_gap;
        // chain->score(chain_score);
        chain->score(score);
        int32_t cov_bases_q = 0;
        int32_t cov_bases_t = 0;

        raptor::mapper::CalcHitCoverage(chain->hits(), k, 0, chain->hits().size(), cov_bases_q, cov_bases_t);
        chain->cov_bases_q(cov_bases_q);
        chain->cov_bases_t(cov_bases_t);

        int32_t qspan = chain->hits().back().QueryPos() - chain->hits().front().QueryPos();
        double frac = (qspan == 0) ? 0 : ((double) cov_bases_q) / ((double) qspan);

        // Add the new chain.
        if (cov_bases_q >= min_cov_bases && cov_bases_t >= min_cov_bases) {
            chains.emplace_back(chain);
        }
        /////////////////////////

    }

    #ifdef DEBUG_DP_VERBOSE_
        printf("The DP:\n");
        for (int32_t i = 0; i < dp.size(); i++) {
            printf("[%d] dp[i] = %d, pred[i] = %d, chain_id[i] = %d\n", i, dp[i], pred[i],
                   chain_id[i]);
        }
    #endif

    return chains;
}

}  // namespace raptor
