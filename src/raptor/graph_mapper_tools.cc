/*
 * graph_mapper_tools.cc
 *
 *  Created on: Jan 26, 2019
 *      Author: Ivan Sovic
 */

#include <raptor/graph_mapper_tools.h>
#include <raptor/mapper_tools.h>

#include <algorithm>
#include <cmath>
#include <log/log_tools.h>
#include <raptor/mapper_tools.h>

namespace raptor {
namespace graphmapper {

constexpr int64_t PlusInf64 = std::numeric_limits<int64_t>::max() - 10000; // Leave a margin.

GraphEdgePenalties CalcAnchorGraphEdgePenalty(const raptor::AnchorPtr& node_v, const raptor::AnchorPtr& node_w,
                    const std::shared_ptr<raptor::AnchorGraphEdge>& local_graph_edge,
                    int32_t node_v_score,       // Source node (predecessor) score. This is _not_ the DP value, but the score of that single node. Required only for score compensation.
                    int32_t node_w_score,        // Sink node score.
                    bool verbose_debug
                    ){
                    // const std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>>& e) {

    // int32_t node_w_score = (node_w->Score() > 0) ? node_w->Score() : node_w->QuerySpan(); // CoveredBasesQuery(); // x_w_span;
    // int32_t node_v_score = (node_v->Score() > 0) ? node_v->Score() : node_v->QuerySpan(); // CoveredBasesQuery();

    GraphEdgePenalties ret;

    // const double lin_factor = 0.01 * 15;
    const double lin_factor = 1.0;

    // Sink data. Easier to configure.
    ret.w_qstart = node_w->QueryStart();
    ret.w_qend = node_w->QueryEnd();
    ret.w_tstart = node_w->TargetStart();
    ret.w_tend = node_w->TargetEnd();
    ret.w_l = ret.w_tstart - ret.w_qstart;

    ret.v_qstart = node_v->QueryStart();
    ret.v_qend = node_v->QueryEnd();
    ret.v_tstart = node_v->TargetStart();
    ret.v_tend = node_v->TargetEnd();
    ret.v_l = ret.v_tstart - ret.v_qstart;
    ret.v_virtual_tstart = node_v->TargetStart();       // Valid for implicit edges.
    ret.v_virtual_tend = node_v->TargetEnd();           // Valid for implicit edges.
    ret.v_virtual_l = ret.v_virtual_tend - ret.v_qend;

    ret.v_pivot_tpos = (local_graph_edge->IsImplicit() == false) ?
                            std::max((int32_t) node_v->TargetStart(), std::min((int32_t) node_v->TargetEnd(), (int32_t) local_graph_edge->source_start())) :
                            node_v->TargetEnd();

    ret.v_score_compensation = 0;

    if (local_graph_edge->IsImplicit() == false) {
        // Implicit edges are between nodes on the same segment. There
        // is no actual SegmentEdge object between them.
        // For implicit edges, virtual target position is the same
        // as the actual anchor target position.
        // Everything is already calculated above.

        // In case there are multiple consecutive SegmentEdge objects
        // without an anchor in between any of them, we still need
        // to shift the target coordinate.
        int64_t jump_length = local_graph_edge->CalcJumpLength();

        // This shifts the predecessors (node_v) coordinate to node_w's space.
        // Using local_graph_edge->source_end() should allow for edge overlaps.
        // This is the legacy formulation:
        ret.v_virtual_tstart = local_graph_edge->sink_start() - (local_graph_edge->source_start() - node_v->TargetStart()) - jump_length;
        ret.v_virtual_tend = ret.v_virtual_tstart + ret.v_tspan();
        ret.v_pivot_tpos = (local_graph_edge->sink_start() - local_graph_edge->source_start() - jump_length) + ret.v_pivot_tpos;
        ret.v_virtual_l = ret.v_virtual_tend - ret.v_qend;

        ret.v_score_compensation = (local_graph_edge->source_start() >= node_v->TargetStart() && local_graph_edge->source_start() < node_v->TargetEnd() && ret.v_tspan() != 0) ?
                                ((int32_t) (node_v_score * (static_cast<double>(node_v->TargetEnd() - local_graph_edge->source_start())) / (static_cast<double>(ret.v_tspan())))) : 0;
    }

    // int64_t dist_min = (int64_t) abs(std::min(abs(dist_y), abs(dist_x)));
    // int64_t mismatch_dist = dist_min - score_compensation;

    int64_t gap_dist = ret.gap_dist();
    // std::cerr << "ret.v_virtual_tstart = " << ret.v_virtual_tstart << ", ret.v_virtual_tend = " << ret.v_virtual_tend << "\n";

    int64_t lin_part = (int32_t) (gap_dist * lin_factor); // * 0.01 * 15); // x_w_span);
    int64_t log_part = (int32_t) ((gap_dist == 0) ? 0 : log2(gap_dist));    // Before changing the log to log2 and moving the "(dist_y > seed_voin_dist) break" above, there was only 1 off.
    ret.edge_score = -(lin_part + (log_part >> 1));
    ret.score_vw = node_w_score + ret.edge_score - ret.v_score_compensation;

    ret.status_flag = 0;

    // Sanity check that the coordinates make sense.
    if (ret.w_tend <= ret.v_pivot_tpos) {
        ret.status_flag |= (1 << 0);
    }
    if (ret.w_qstart <= ret.v_qstart) {
        ret.status_flag |= (1 << 1);
    }
    if (ret.w_tstart <= ret.v_virtual_tstart) {
        ret.status_flag |= (1 << 2);
    }


    #ifdef RAPTOR_TESTING_MODE
        if (verbose_debug) {
            LOG_NOHEADER("v_qstart = %ld, v_qend = %ld, v_virtual_tstart = %ld, v_virtual_tend = %ld, gap_dist = %ld, edge_score = %lf, v_score_compensation = %ld, score_vw = %lf\n",
                ret.v_qstart, ret.v_qend, ret.v_virtual_tstart, ret.v_virtual_tend, gap_dist, ret.edge_score, ret.v_score_compensation, ret.score_vw);
        }
    #endif

    return ret;


            // Consider this for overlapping anchors (this is taken from dp_chain.cc): int32_t x_j_score = std::min(k, (int32_t) std::min(abs(x_i_start - x_j_start), abs(y_i_start - y_j_start)));

            // #ifdef RAPTOR_TESTING_MODE
            //     if (verbose_debug_qid) {
            //         LOG_NOHEADER("  [node_j_name = %ld -> node_i_name = %ld]\tvirtual_y_j_start = %ld, virtual_y_j_end = %ld, pivot_pos_y = %d, x_j_start = %ld, x_j_end = %ld, l_j = %ld, mismatch_dist = %ld, gap_dist = %ld, chain_max_bandwidth = %l\n",
            //                 node_j_name, node_i_name, virtual_y_j_start, virtual_y_j_end, pivot_pos_y, x_j_start, x_j_end, l_j, mismatch_dist, gap_dist, chain_max_bandwidth);
            //         LOG_NOHEADER("                                   \ty_i_start = %ld, y_i_end = %ld, x_i_start = %ld, x_i_end = %ld, l_i = %ld\n",
            //                 y_i_start, y_i_end, x_i_start, x_i_end, l_i);
            //     }
            // #endif

            // #ifdef RAPTOR_TESTING_MODE
            //         if (verbose_debug_qid) {
            //         LOG_NOHEADER("                                   \t");
            //         LOG_NOHEADER("lin_part = %ld, log_part = %ld, x_i_span = %ld, gap_dist = %ld, edge_score = %ld, score_ij = %ld\n", lin_part, log_part, x_i_span, gap_dist, edge_score, score_ij);
            //         LOG_NOHEADER("                                   \t");
            //         LOG_NOHEADER("score_ij = %ld, predecessor_score = %ld, edge_score = %ld, new_dp_val = %ld, new_dp_path = %ld, node_i_score = %d, node_j_score = %d, gap_dist = %ld, mismatch_dist = %ld, score_compensation = %d\n",
            //                 score_ij, dp[node_j_name + 1], edge_score, new_dp_val, new_dp_path, node_i_score, node_j_score, gap_dist, mismatch_dist, score_compensation);
            //     }
            // #endif

}

std::shared_ptr<raptor::LinearMappingResult> BreakAnchors(
                    const raptor::GraphPtr graph,
                    const std::shared_ptr<raptor::LinearMappingResult>& in_linear_result,
                    const int32_t k
                    ) {

    std::shared_ptr<raptor::LinearMappingResult> out_linear_result = raptor::createMappingResult(
            in_linear_result->QueryId(), in_linear_result->QueryLen(),
            in_linear_result->QueryHeader(), in_linear_result->Index());

    std::vector<raptor::ChainPtr> filtered_hits;

    /*
     * For each target anchor, look if there are any out-/in- edges forking out of that
     * sequence between it's start and end coordinates.
     * If there are, we need to split the anchor for better graph chaining afterwards.
    */
    for (auto& single_target_anchors: in_linear_result->target_anchors()) {

        for (size_t anchor_ordinal_id = 0; anchor_ordinal_id < single_target_anchors->hits().size(); anchor_ordinal_id++) {
            auto& anchor = single_target_anchors->hits()[anchor_ordinal_id];
            int32_t lookup_t_id = anchor->TargetID();
            int32_t lookup_start = anchor->TargetStart();
            int32_t lookup_end = anchor->TargetEnd();

            std::vector<int32_t> breakpoints;

            // Find all input edge junctions, and accumulate their coordinates.
            std::vector<std::shared_ptr<raptor::SegmentEdge>> avail_in_edges = graph->FindSegmentInputEdges(lookup_t_id, lookup_start, lookup_end);
            for (auto& edge: avail_in_edges) {
                // Take care of the strand.
                if (anchor->TargetRev() != edge->sink_is_rev()) {
                    continue;
                }

                breakpoints.emplace_back(edge->sink_start());

                if (edge->sink_end() != edge->sink_start()) {
                    breakpoints.emplace_back(edge->sink_end());
                }
            }

            // Find all output edge junctions, and accumulate their coordinates.
            std::vector<std::shared_ptr<raptor::SegmentEdge>> avail_out_edges = graph->FindSegmentOutputEdges(lookup_t_id, lookup_start, lookup_end);
            for (auto& edge: avail_out_edges) {
                // Take care of the strand.
                if (anchor->TargetRev() != edge->source_is_rev()) {
                    continue;
                }

                breakpoints.emplace_back(edge->source_start());

                if (edge->source_end() != edge->source_start()) {
                    breakpoints.emplace_back(edge->source_end());
                }
            }

            std::sort(breakpoints.begin(), breakpoints.end());

            auto anchor_th_id = anchor->target_hits_id();
            auto& th = in_linear_result->target_hits()[anchor_th_id];

            // #ifdef RAPTOR_TESTING_MODE
            //     LOG_ALL("breakpoints.size() = %d\n", breakpoints.size());
            //     for (int32_t bp_id = 0; bp_id < breakpoints.size(); bp_id++) {
            //         LOG_ALL("Breakpoint %d: %d\n", bp_id, breakpoints[bp_id]);
            //     }
            // #endif

            // Copy the seeds into new anchors, broken by breakpoints.
            size_t curr_breakpoint = 0;
            auto new_chain = raptor::ChainPtr(new raptor::ChainType(single_target_anchors->env()));
            for (size_t hit_id = 0; hit_id < th->hits().size(); ++hit_id) {
                auto& hit = th->hits()[hit_id];

                size_t prev_breakpoint = curr_breakpoint;

                // The "<=" in the while loop definition below is intentional and necessary.
                // The breakpoint should always be larger than the hit.TargetPos() because
                // we're binning everything before it on the target.
                while (curr_breakpoint < breakpoints.size() && breakpoints[curr_breakpoint] <= hit.TargetPos()) {
                    ++curr_breakpoint;
                }

                // The breakpoint ID changed. Summarize the old group and create an empty one.
                if (curr_breakpoint != prev_breakpoint && new_chain->hits().size() > 0) {
                    int32_t cov_bases_q = 0, cov_bases_t = 0;
                    raptor::mapper::CalcHitCoverage(new_chain->hits(), k, 0, new_chain->hits().size(), cov_bases_q, cov_bases_t);
                    new_chain->cov_bases_q(cov_bases_q);
                    new_chain->cov_bases_t(cov_bases_t);
                    new_chain->score(std::min(cov_bases_q, cov_bases_t));
                    if (new_chain->hits().size() > 0) {
                        filtered_hits.emplace_back(new_chain);
                        new_chain = raptor::ChainPtr(new raptor::ChainType(single_target_anchors->env()));
                    }
                }

                new_chain->hits().emplace_back(hit);
            }

            if (new_chain->hits().size() > 0) {
                new_chain->score(th->score());
                int32_t cov_bases_q = 0, cov_bases_t = 0;
                raptor::mapper::CalcHitCoverage(new_chain->hits(), k, 0, new_chain->hits().size(), cov_bases_q, cov_bases_t);
                new_chain->cov_bases_q(cov_bases_q);
                new_chain->cov_bases_t(cov_bases_t);
                new_chain->score(std::min(cov_bases_q, cov_bases_t));
                filtered_hits.emplace_back(new_chain);
            }
        }

    }

    auto broken_anchors = raptor::mapper::MakeAnchors(filtered_hits);

    for (auto& ba : broken_anchors) {
        for (auto& anchor : ba->hits()) {
            anchor->score(filtered_hits[anchor->target_hits_id()]->score());
        }
    }

#ifdef RAPTOR_TESTING_MODE
    // if (params->debug_qid == qseq->abs_id() ||
    //     params->debug_qname == std::string(qseq->header())) {
    //     LOG_ALL("Showing debug info for read %d: %s\n", qseq->id(), qseq->header().c_str());

    //     raptor::WriteTargetHits("temp/debug/mapping-3-broken_anchors.csv", filtered_hits, index_->params()->k,
    //                     qseq.get_header(), qseq.get_sequence_length(), std::string("ref"), 0);
    // }
#endif

    out_linear_result->target_hits(filtered_hits);
    out_linear_result->target_anchors(broken_anchors);

    return out_linear_result;
}

/*
 * Finds leaf nodes for every labeled path. It doesn't return _all_ leaf nodes
 * but only those which have the score >= (allowed_score_diff_frac * max_leaf_node_score).
 * Of course, if allowed_score_diff_frac == 0.0, all leaf nodes will be reported.
 *
 * Returns a vector of vectors. Inner vector is the list of leaf nodes for a path.
 * The outter vector is one element per path. The ID of the element corresponds to the
 * path ID.
*/
std::vector<std::vector<int64_t>> FindBestLeafNodes(
                    const std::shared_ptr<raptor::AnchorGraph>& graph,
                    double allowed_score_diff_frac,
                    bool verbose_debug_qid
                    ) {

    // Find out the maximum number of paths in the graph. Needed to preallocate space for leafs.
    int32_t num_paths = 0;
    for (const auto& v_item: graph->nodes()) {
        if (v_item->is_removed()) {
            continue;
        }
        num_paths = std::max(num_paths, v_item->data()->path_id() + 1);
    }

    if (num_paths == 0) {
        return {};
    }

    std::vector<int64_t> max_pid_leaf_scores(num_paths, -PlusInf64);
    std::vector<std::vector<double>> pid_leaf_scores(num_paths, std::vector<double>());
    std::vector<std::vector<int64_t>> pid_leafs(num_paths, std::vector<int64_t>());

    // Find all the leaf scores.
    for (const auto& v_item: graph->nodes()) {
        int64_t v_name = v_item->name();
        if (v_item->is_removed()) {
            continue;
        }
        const auto& v_data = v_item->data();
        int32_t pid = v_data->path_id();

        // All nodes should have a non-negative path ID, but let's be safe.
        if (pid < 0) {
            LOG_ALL("The pid is < 0! pid = %ld. Skipping.\n", pid);
            continue;
        }

        // Count the number of valid out-edges.
        // This means: non-deleted edges, which belong to the same path.
        // Edges which do not belong to any paths should be discarded.
        auto edges = graph->GetOutEdges(v_name);
        int64_t num_valid_out_edges = 0;
        for (const auto& e_item: edges) {
            if (e_item->data()->path_id() < 0) {
                continue;
            }
            if (e_item->data()->path_id() != pid) {
                continue;
            }
            ++num_valid_out_edges;
        }

        // Take only the sinks.
        if (num_valid_out_edges == 0) {
            pid_leafs[pid].emplace_back(v_name);
            pid_leaf_scores[pid].emplace_back(v_data->score());
            if (v_data->score() > max_pid_leaf_scores[pid]) {
                max_pid_leaf_scores[pid] = v_data->score();
            }
        }
    }

    // For each path, remove poor scoring leafs.
    for (size_t pid = 0; pid < pid_leafs.size(); pid++) {
        // Using the abs because if the score is below zero, our minimum threshold needs to be
        // even below that.
        double min_allowed_leaf_score = max_pid_leaf_scores[pid] - std::abs(max_pid_leaf_scores[pid]) * (1.0 - allowed_score_diff_frac);
        std::vector<int64_t> filtered_pid_leafs;
        for (size_t leaf_id = 0; leaf_id < pid_leafs[pid].size(); ++leaf_id) {
            if (pid_leaf_scores[pid][leaf_id] >= min_allowed_leaf_score) {
                filtered_pid_leafs.emplace_back(pid_leafs[pid][leaf_id]);
            }
        }
        std::swap(pid_leafs[pid], filtered_pid_leafs);
    }

    #ifdef RAPTOR_TESTING_MODE
        if (verbose_debug_qid) {
            LOG_NEWLINE;
            LOG_ALL("Max scoring leafs for each path:\n");
            for (size_t pid = 0; pid < pid_leafs.size(); pid++) {
                LOG_NOHEADER("  - [path_id = %ld]: ", pid);
                auto& path_dp_ids = pid_leafs[pid];
                for (size_t j = 0; j < path_dp_ids.size(); j++) {
                    LOG_NOHEADER("%ld ", path_dp_ids[j]);
                }
                LOG_NOHEADER("\n");
                // LOG_NOHEADER("      - max_pid_leaf_node = %ld, max_pid_leaf_score = %ld\n", max_pid_leaf_nodes[pid], max_pid_leaf_scores[pid]);
                LOG_NOHEADER("      - max_pid_leaf_score = %ld\n", max_pid_leaf_scores[pid]);
            }
            LOG_NOHEADER("\n");

        }
    #endif

    return pid_leafs;
}

}
}
