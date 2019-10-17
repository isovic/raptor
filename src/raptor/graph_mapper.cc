/*
 * graph_mapper.cc
 *
 *  Created on: Dec 18, 2017
 *      Author: Ivan Sovic
 */

#include <raptor/graph_mapper.h>
#include <raptor/graph_mapper_tools.h>
#include <raptor/mapper_tools.h>

#include <cmath>
#include <algorithm/graph.hpp>
#include <graph/segment_edge.h>
#include <algorithm>
#include <memory>
#include <set>
#include <tuple>
#include <unordered_map>
#include <raptor/interval_tree_builder.h>
#include <log/log_tools.h>
#include <raptor/mapper.h>
#include <raptor/dp_chain.h>
#include <debug_tools/write_seed_hit_1.h>
#include <utility/tictoc.h>

namespace raptor {

// A const used in DP.
static const int64_t PlusInf64 = std::numeric_limits<int64_t>::max() - 10000; // Leave a margin.

typedef std::tuple<raptor::AnchorPtr, int64_t, int64_t> AnchorTuple;
typedef std::vector<AnchorTuple> AnchorTupleVector;
typedef std::vector<std::shared_ptr<raptor::SegmentEdge>> PathType;

std::unique_ptr<raptor::GraphMapper> createGraphMapper(
                                        const mindex::IndexPtr index,
                                        const raptor::GraphPtr graph,
                                        const raptor::SplitSegmentGraphPtr& ssg,
                                        const std::shared_ptr<raptor::ParamsMapper> params) {
  return std::unique_ptr<raptor::GraphMapper>(new raptor::GraphMapper(index, graph, ssg, params));
}

GraphMapper::GraphMapper(
                    const mindex::IndexPtr index,
                    const raptor::GraphPtr graph,
                    const raptor::SplitSegmentGraphPtr& ssg,
                    const std::shared_ptr<raptor::ParamsMapper> params)
                    :   index_(index),
                        graph_(graph),
                        ssg_(ssg),
                        params_(params)
                    {

}

GraphMapper::~GraphMapper() {

}

std::shared_ptr<raptor::GraphMappingResult> GraphMapper::Map(const mindex::SequencePtr& qseq,
                                                         const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result) {

    bool verbose_debug_qid = false;
    DEBUG_QSEQ(params_, qseq, verbose_debug_qid = true);

    auto result = raptor::createGraphMappingResult(qseq->abs_id(),
                                               qseq->data().size(),
                                               qseq->header(), index_);

    TicToc tt_total;
    tt_total.start();

    // auto& target_anchors = input_mapping_result->target_anchors();
    TicToc tt_break_anchors;
    tt_break_anchors.start();
    std::shared_ptr<raptor::LinearMappingResult> final_linear_mapping = raptor::graphmapper::BreakAnchors(
                                    graph_,
                                    input_mapping_result,
                                    index_->params()->k);
    tt_break_anchors.stop();

    TicToc tt_align_anchors;
    tt_align_anchors.start();
    if (params_->graph_score_anchors) {
        AlignAnchors_(qseq, final_linear_mapping->target_anchors());
    }
    tt_align_anchors.stop();

    TicToc tt_construct_anchor_graph;
    tt_construct_anchor_graph.start();
    auto anchor_graph = ConstructAnchorGraph(
                                            graph_,
                                            final_linear_mapping->target_anchors(),
                                            params_->chain_max_dist,
                                            params_->chain_max_dist,
                                            params_->graph_max_path_edges,
                                            params_->chain_max_skip,
                                            verbose_debug_qid);
    tt_construct_anchor_graph.stop();

    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
            LOG_NOHEADER("%s\n", anchor_graph->Verbose().c_str());
            // LOG_ALL("graph_:\n%s\n", graph_->Verbose().c_str());
            LOG_ALL("GraphDP.\n");
        }
    #endif

    TicToc tt_graph_dp;
    tt_graph_dp.start();

    raptor::AnchorGraphPtr mapped_anchor_graph = GraphDP2(anchor_graph, params_->chain_max_bandwidth, verbose_debug_qid);

    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
            LOG_NOHEADER("After GraphDP2:\n%s\n", mapped_anchor_graph->Verbose().c_str());
        }
    #endif

    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs = CreateAnchorGraphConnectedComponents(
                                            mapped_anchor_graph,
                                            params_->graph_allowed_score_diff_frac,
                                            verbose_debug_qid);

    std::vector<std::shared_ptr<raptor::LocalPath>> paths;

    for (size_t pgraph_id = 0; pgraph_id < pgraphs.size(); ++pgraph_id) {
        auto& vals = pgraphs[pgraph_id];
        raptor::AnchorGraphPtr pid_graph = std::get<0>(vals);
        int64_t pid = std::get<1>(vals);
        double pid_score = std::get<2>(vals);

        #ifdef RAPTOR_TESTING_MODE
            if (params_->debug_qid == qseq->abs_id() ||
                params_->debug_qname == std::string(qseq->header())) {
                LOG_ALL("[pid = %ld] Graph component before backtrack:\n%s\n", pid, pid_graph->Verbose().c_str());
            }
        #endif

        // The following will create linear paths from a given graph.
        std::unordered_map<int64_t, int64_t> node_to_path_score;
        std::vector<raptor::AnchorGraphPtr> pid_path_graphs = BacktrackReducedMappedAnchorGraph(
                                                                    pid_graph,
                                                                    params_->graph_allowed_score_diff_frac,
                                                                    node_to_path_score,
                                                                    verbose_debug_qid);

        for (size_t pid_graph_id = 0; pid_graph_id < pid_path_graphs.size(); ++pid_graph_id) {
            auto new_paths = GraphToPaths2_(pid_path_graphs[pid_graph_id], node_to_path_score, verbose_debug_qid);
            paths.insert(paths.end(), new_paths.begin(), new_paths.end());
        }

        #ifdef RAPTOR_TESTING_MODE
            if (params_->debug_qid == qseq->abs_id() ||
                params_->debug_qname == std::string(qseq->header())) {
                LOG_ALL("[after BacktrackReducedMappedAnchorGraph] pid = %ld, paths.size() = %ld\n", pid, paths.size());
        }
        #endif
    }

    tt_graph_dp.stop();

    result->paths(paths);
    result->SetReturnValue(raptor::MapperReturnValueBase::OK);

    tt_total.stop();

    result->timings()["gmap_total"] = tt_total.get_millisecs();
    result->timings()["gmap_break_anchors"] = tt_break_anchors.get_millisecs();
    result->timings()["gmap_align_anchors"] = tt_align_anchors.get_millisecs();
    result->timings()["gmap_constr_anchor_graph"] = tt_construct_anchor_graph.get_millisecs();
    result->timings()["gmap_graph_dp"] = tt_graph_dp.get_millisecs();

    return result;
}

std::shared_ptr<raptor::GraphMappingResult> GraphMapper::MapBetter(const mindex::SequencePtr& qseq,
                                                         const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result) {

    bool verbose_debug_qid = false;
    DEBUG_QSEQ(params_, qseq, verbose_debug_qid = true);

    auto result = raptor::createGraphMappingResult(qseq->abs_id(),
                                               qseq->data().size(),
                                               qseq->header(), index_);

    TicToc tt_total;
    tt_total.start();

    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
                LOG_ALL("Debug verbose of input_mapping_result->target_anchors():\n");

                LOG_NOHEADER("%s\n", input_mapping_result->WriteAsCSV(',').c_str());

            // // LOG_NOHEADER("TargetAnchors:\n%s\n", input_mapping_result->target_anchors()->Verbose().c_str());
            // for (size_t i = 0; i < input_mapping_result->target_anchors().size(); ++i) {
            //     const auto& single_target_anchors = input_mapping_result->target_anchors()[i];
            // // for (auto& single_target_anchors: input_mapping_result->target_anchors()) {
            //     LOG_NOHEADER("TargetAnchors [%ld] for TargetID = %ld:\n", i, single_target_anchors->env()->t_id);
            //     LOG_NOHEADER("%s\n", single_target_anchors->VerbosePointers().c_str());

            //     // for (size_t anchor_ordinal_id = 0; anchor_ordinal_id < single_target_anchors->hits().size(); anchor_ordinal_id++) {
            //     //     auto& anchor = single_target_anchors->hits()[anchor_ordinal_id];
            //     // }
            // }

        }
    #endif

    // auto& target_anchors = input_mapping_result->target_anchors();
    TicToc tt_break_anchors;
    tt_break_anchors.start();
    std::shared_ptr<raptor::LinearMappingResult> final_linear_mapping = raptor::graphmapper::BreakAnchors(
                                    graph_,
                                    input_mapping_result,
                                    index_->params()->k);
    tt_break_anchors.stop();

    TicToc tt_align_anchors;
    tt_align_anchors.start();
    if (params_->graph_score_anchors) {
        AlignAnchors_(qseq, final_linear_mapping->target_anchors());
    }
    tt_align_anchors.stop();

    TicToc tt_construct_anchor_graph;
    tt_construct_anchor_graph.start();

    auto anchor_graph = ConstructAnchorGraph(
                                            graph_,
                                            final_linear_mapping->target_anchors(),
                                            params_->chain_max_dist,
                                            params_->chain_max_dist,
                                            params_->graph_max_path_edges,
                                            params_->chain_max_skip,
                                            verbose_debug_qid);
    tt_construct_anchor_graph.stop();

    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
            LOG_NOHEADER("%s\n", anchor_graph->Verbose().c_str());
            // LOG_ALL("graph_:\n%s\n", graph_->Verbose().c_str());
            LOG_ALL("GraphDP.\n");
        }
    #endif

    TicToc tt_graph_dp;
    tt_graph_dp.start();

    // auto paths = GraphDP_(anchor_graph, params_->chain_max_bandwidth, verbose_debug_qid);

    raptor::AnchorGraphPtr mapped_anchor_graph = GraphDP2(anchor_graph, params_->chain_max_bandwidth, verbose_debug_qid);

    #ifdef RAPTOR_TESTING_MODE
        if (params_->debug_qid == qseq->abs_id() ||
            params_->debug_qname == std::string(qseq->header())) {
            LOG_NOHEADER("After GraphDP2:\n%s\n", mapped_anchor_graph->Verbose().c_str());
        }
    #endif

    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs = CreateAnchorGraphConnectedComponents(
                                            mapped_anchor_graph,
                                            params_->graph_allowed_score_diff_frac,
                                            verbose_debug_qid);

    /*
     * If use_approx_aln == false, then:
     *  - backtracking will be applied to each AnchorGraph, and the graph broken into
     *    one or more paths.
     *
    */
    bool use_approx_aln = false;

    if (use_approx_aln) {
        /*
         * In this case, break-up every graph component into linear paths.
        */
        std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> new_pgraphs;

        for (size_t pgraph_id = 0; pgraph_id < pgraphs.size(); ++pgraph_id) {
            auto& vals = pgraphs[pgraph_id];
            raptor::AnchorGraphPtr pid_graph = std::get<0>(vals);
            int64_t pid = std::get<1>(vals);
            double pid_score = std::get<2>(vals);

            LOG_ALL("[pid = %ld] Graph component before backtrack:\n%s\n", pid, pid_graph->Verbose().c_str());

            // The BacktrackReducedMappedAnchorGraph will create linear paths from a given graph.
            std::unordered_map<int64_t, int64_t> node_to_path_score;
            std::vector<raptor::AnchorGraphPtr> pid_path_graphs = BacktrackReducedMappedAnchorGraph(
                                                                        pid_graph,
                                                                        params_->graph_allowed_score_diff_frac,
                                                                        node_to_path_score,
                                                                        verbose_debug_qid);

            for (size_t pid_graph_id = 0; pid_graph_id < pid_path_graphs.size(); ++pid_graph_id) {
                new_pgraphs.emplace_back(std::make_tuple(pid_path_graphs[pid_graph_id], pid, pid_score));
            }

            LOG_ALL("[after BacktrackReducedMappedAnchorGraph] pid = %ld, pid_path_graphs.size() = %ld\n", pid, pid_path_graphs.size());
        }

        std::swap(pgraphs, new_pgraphs);
    }

    DEBUG_RUN(verbose_debug_qid, LOG_ALL("About to run ConstructAlignmentGraph on %ld pgraph components.\n", pgraphs.size()));

    // Convert each graph component into the alignment graph.
    std::vector<std::pair<raptor::SplitSegmentGraphPtr, double>> aln_graphs;
    for (size_t pgraph_id = 0; pgraph_id < pgraphs.size(); ++pgraph_id) {
        auto& vals = pgraphs[pgraph_id];
        raptor::AnchorGraphPtr pid_graph = std::get<0>(vals);
        int64_t pid = std::get<1>(vals);
        double pid_score = std::get<2>(vals);

        raptor::SplitSegmentGraphPtr aln_graph = GraphMapper::ConstructAlignmentGraph(ssg_, pid_graph, qseq, use_approx_aln, verbose_debug_qid);
        aln_graphs.emplace_back(std::make_pair(aln_graph, pid_score));

        #ifdef RAPTOR_TESTING_MODE
            if (params_->debug_qid == qseq->abs_id() ||
                params_->debug_qname == std::string(qseq->header())) {
                LOG_ALL("Writing graphs to temp/debug-graph/*.\n");

                std::ofstream ofs_debug_graph_anchor("temp/debug-graph/pid_" + std::to_string(pid) + ".graph.anchor.html");
                if (ofs_debug_graph_anchor.is_open()) {
                    pid_graph->WriteAsHTML(ofs_debug_graph_anchor);
                }

                std::ofstream ofs_debug_graph_aln("temp/debug-graph/pid_" + std::to_string(pid) + ".graph.aln.html");
                if (ofs_debug_graph_aln.is_open()) {
                    aln_graph->WriteAsHTML(ofs_debug_graph_aln);
                }
            }
        #endif
    }

    tt_graph_dp.stop();

    // result->paths(paths);
    result->SetReturnValue(raptor::MapperReturnValueBase::OK);

    tt_total.stop();

    result->timings()["gmap_total"] = tt_total.get_millisecs();
    result->timings()["gmap_break_anchors"] = tt_break_anchors.get_millisecs();
    result->timings()["gmap_align_anchors"] = tt_align_anchors.get_millisecs();
    result->timings()["gmap_constr_anchor_graph"] = tt_construct_anchor_graph.get_millisecs();
    result->timings()["gmap_graph_dp"] = tt_graph_dp.get_millisecs();

    return result;
}

std::shared_ptr<raptor::GraphMappingResult> GraphMapper::DummyMap(const mindex::SequencePtr& qseq,
                                                              const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result) {
    TicToc tt_total;
    tt_total.start();

    auto result = raptor::createGraphMappingResult(qseq->abs_id(),
                                               qseq->data().size(),
                                               qseq->header(), index_);

    std::vector<std::shared_ptr<LocalPath>> paths;

    for (auto& single_target_anchors: input_mapping_result->target_anchors()) {
        for (size_t anchor_ordinal_id = 0; anchor_ordinal_id < single_target_anchors->hits().size(); anchor_ordinal_id++) {
            auto& anchor = single_target_anchors->hits()[anchor_ordinal_id];
            auto path = LocalPath::createPath(static_cast<int64_t>(paths.size()), anchor);
            path->score(anchor->Score());
            paths.emplace_back(path);
        }
    }

    result->paths(paths);
    result->SetReturnValue(raptor::MapperReturnValueBase::OK);

    tt_total.stop();

    result->timings()["gmap_total"] = tt_total.get_millisecs();
    result->timings()["gmap_break_anchors"] = 0.0;
    result->timings()["gmap_align_anchors"] = 0.0;
    result->timings()["gmap_constr_anchor_graph"] = 0.0;
    result->timings()["gmap_graph_dp"] = 0.0;

    return result;
}

/*
Algorithm for constructing the AnchorGraph:
1. Add implicit edges for every anchor. Implicit edges are edges within the same target, and are represented as LocalEdges composed of an empty vector of SegmentEdges.
2. BFS from every node (sink) to find the local edges in the graph:
    2.1 Discover all the nearby edges and add them to the processing queue.
    2.2 Perform a BFS starting from each SegmentEdge in the processing queue.
        2.2.1 For each edge in the processing queue, discover nodes in front of it.
        2.2.2 For every discovered node (source), create a LocalEdge object which holds a list of all SegmentEdges traversed between the two nodes.
        2.2.3 Extend the queue with all reachable SegmentEdges from the current starting SegmentEdge. This will ensure continuation of the BFS.
        2.2.4 Break the BFS when a depth limit is reached (e.g. distance in bp from the initial sink node).
*/

void AddImplicitEdges_(const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors,
                       const std::unordered_map<raptor::AnchorPtr, int64_t>& anchor_to_oid,
                       int32_t max_allowed_dist,
                       int32_t chain_max_skip,
                       std::shared_ptr<raptor::AnchorGraph> local_graph,
                       bool verbose_debug_qid) {

    // Add implicit edges.
    for (auto& single_target_anchors: target_anchors) {
        for (int32_t i = 1; i < (int32_t) single_target_anchors->hits().size(); i++) {
            auto& anchor = single_target_anchors->hits()[i];

            for (int32_t j = (i - 1); j >= 0 && j > (i - chain_max_skip); j--) {
                auto& pred_anchor = single_target_anchors->hits()[j];
                bool is_rev = pred_anchor->TargetRev() != anchor->TargetRev();
                // Skip hits on the wrong strand.
                if (is_rev) {
                    continue;
                }
                // Avoid meaningless cycles.
                if (pred_anchor->QueryEnd() > anchor->QueryStart()) {
                    continue;
                }
                // Unlike the discovery of nodes via edges (as bellow), here
                // we need to check the target coordinate as well because the
                // nodes are on the same segment.
                if (pred_anchor->TargetEnd() > anchor->TargetStart()) {
                    continue;
                }
                // Skip very distant edges.
                if ((anchor->QueryStart() - pred_anchor->QueryEnd()) > max_allowed_dist) {
                    continue;
                }
                // Skip very distant edges.
                if ((anchor->TargetStart() - pred_anchor->TargetEnd()) > max_allowed_dist) {
                    continue;
                }

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid,
                        LOG_ALL("Implicit edges between:\n  - %s\n  - %s\n", anchor->Verbose().c_str(), pred_anchor->Verbose().c_str())
                    );
                #endif

                auto new_implicit_edge = raptor::createImplicitAnchorGraphEdge(is_rev);
                auto edge_scores = raptor::graphmapper::CalcAnchorGraphEdgePenalty(pred_anchor, anchor, new_implicit_edge, 0, 0, verbose_debug_qid);
                new_implicit_edge->score(edge_scores.edge_score);
                local_graph->AddEdge(pred_anchor->id(), anchor->id(), new_implicit_edge);

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid,
                        LOG_NOHEADER("\n");
                    );
                #endif
            }
        }
    }

}

void DiscoverPredecessorNodes_(
                               const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors,
                               const std::unordered_map<raptor::AnchorPtr, int64_t>& anchor_to_oid,
                               const std::unordered_map<int64_t, int64_t>& tid_to_taid,             // Target ID to target_anchors ID.
                               const PathType& path,
                               int32_t source_t_id,
                               int32_t source_lookup_start,
                               int32_t source_lookup_end,
                               bool source_t_is_rev,
                               const raptor::AnchorPtr& anchor,
                               int32_t max_allowed_dist,
                               int32_t chain_max_skip,
                               std::shared_ptr<raptor::AnchorGraph> local_graph,
                               bool verbose_debug_qid) {

    auto edge = path.back();

    int64_t t_id_with_rev = source_t_id << 1;
    t_id_with_rev |= (source_t_is_rev) ? 1 : 0;

    auto it_taid = tid_to_taid.find(t_id_with_rev);

    if (it_taid == tid_to_taid.end()) {
        return;
    }

    int64_t taid = it_taid->second;
    auto& single_target_anchors = target_anchors[taid];

    #ifdef RAPTOR_TESTING_MODE
        DEBUG_RUN(verbose_debug_qid,
            LOG_ALL("  Running a binary search on taid = %ld, source_lookup_end = %ld.\n", taid, source_lookup_end)
        );

        DEBUG_RUN(verbose_debug_qid,
            LOG_ALL("  Hits used for binary search:\n%s\n", single_target_anchors->VerbosePointers().c_str())
        );

    #endif

    // Find the first anchor with the start coordinate _larger_ than the source_lookup_end.
    // If all values are smaller, it returns .end(). If all values are larger, then it returns .begin().
    // Then, we'll iterate in reverse from the last anchor until we find an anchor out of frame.
    auto it_anchor = std::upper_bound(single_target_anchors->hits().begin(),
                                    single_target_anchors->hits().end(),
                                    source_lookup_end,
                                    [](int32_t val, const std::shared_ptr<raptor::RegionMapped>& a) { return val < a->TargetStart(); });

    int64_t num_predecessors = 0;
    int64_t num_nodes_added = 0;
    int64_t count = it_anchor - single_target_anchors->hits().begin();

    for (int64_t ii = count; ii > 0; --ii) {
        --it_anchor;
        auto& pred_anchor = *(it_anchor);
        num_predecessors += 1;

        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid,
                LOG_ALL("        [pred: %ld/%ld] Candidate predecessor anchor (pred_anchor->id() = %ld, anchor->id() = %ld): %s\n",
                        num_predecessors, count, pred_anchor->id(), anchor->id(), pred_anchor->Verbose().c_str())
            );
        #endif

        // Safety heuristic.
        if (num_predecessors > chain_max_skip) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          Breaking because num_predecessors > chain_max_skip (%ld > %ld).\n",
                            num_predecessors, chain_max_skip)
                );
            #endif

            break;
        }

        // This condition is sketchy! Because sorting is on start coordinate and not on end coordinate, there can be a valid end
        // coordinate which comes after this condition turns true!
        // Perhaps the best option would be to revert back to using the interval tree.
        if ((source_lookup_start - pred_anchor->TargetEnd()) > max_allowed_dist && (source_lookup_start - pred_anchor->TargetStart()) > max_allowed_dist) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          Breaking because (source_lookup_end - pred_anchor->TargetEnd()) > max_allowed_dist (%ld > %ld).\n",
                            (source_lookup_end - pred_anchor->TargetEnd()), max_allowed_dist)
                );
            #endif

            break;
        }

        if (pred_anchor->id() == anchor->id()) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          Candidate (pred_anchor->id() = %ld, anchor->id() = %ld) not added because "
                            "they are the same. Anchor: %s\n",
                            pred_anchor->id(), anchor->id(), pred_anchor->Verbose().c_str())
                );

            #endif

            continue;

        }

        // Skip hits on the wrong strand of an edge source.
        if (pred_anchor->TargetRev() != edge->sink_is_rev()) {
            #ifdef RAPTOR_TESTING_MODE
                if (verbose_debug_qid) {
                    LOG_ALL("          Candidate (pred_anchor->id() = %ld, anchor->id() = %ld) not added because "
                            "pred_anchor->TargetRev() != edge->sink_is_rev();"
                            "%d != %d\n",
                            pred_anchor->id(), anchor->id(), (int32_t) pred_anchor->TargetRev(), (int32_t) edge->sink_is_rev());
                }
            #endif

            continue;
        }

        // Avoid meaningless cycles.
        if (pred_anchor->QueryStart() > anchor->QueryStart()) {
            // This used to be (pred_anchor->QueryStart() >= anchor->QueryStart()) but that is bad if there is a perfect circular
            // alignment.
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          Candidate (pred_anchor->id() = %ld, anchor->id() = %ld) not added because "
                            "pred_anchor->QueryStart() > anchor->QueryStart();"
                            "%ld >= %ld\n",
                            pred_anchor->id(), anchor->id(), pred_anchor->QueryStart(), anchor->QueryStart())
                );

            #endif

            continue;
        }
        // Skip very distant edges.
        if ((anchor->QueryStart() - pred_anchor->QueryEnd()) > max_allowed_dist) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          Candidate (pred_anchor->id() = %ld, anchor->id() = %ld) not added because "
                            "(pred_anchor->QueryEnd() - anchor->QueryStart()) > max_allowed_dist;"
                            "(%ld - %ld) = %ld > %ld\n",
                            pred_anchor->id(), anchor->id(), pred_anchor->QueryEnd(),anchor->QueryStart(), pred_anchor->QueryEnd() - anchor->QueryStart(), max_allowed_dist)
                );
            #endif

            continue;
        }
        auto path_copy = path;
        std::reverse(path_copy.begin(), path_copy.end());
        auto new_edge = raptor::createAnchorGraphEdge(path_copy);

        // Deduplicate multiedges.
        auto edge_scores = raptor::graphmapper::CalcAnchorGraphEdgePenalty(pred_anchor, anchor, new_edge, 0, 0, verbose_debug_qid);
        new_edge->score(edge_scores.edge_score);

        auto existing_edges = local_graph->GetEdges(pred_anchor->id(), anchor->id());
        double penalty_max = -1e200;
        for (const auto& existing_edge: existing_edges) {
            penalty_max = std::max(penalty_max, existing_edge->data()->score());
        }

        bool edge_added = false;

        // The  && edge_scores.status_flag == 0 means that all coords check-up went fine.
        if (existing_edges.size() == 0 && edge_scores.status_flag == 0) {
            local_graph->AddEdge(pred_anchor->id(), anchor->id(), new_edge);
            edge_added = true;

        } else if (edge_scores.edge_score > penalty_max && edge_scores.status_flag == 0) {
            // If all existing edges are poorer, mark them as removed.
            local_graph->RemoveEdges(pred_anchor->id(), anchor->id());
            local_graph->AddEdge(pred_anchor->id(), anchor->id(), new_edge);
            edge_added = true;
        }

        #ifdef RAPTOR_TESTING_MODE
            if (edge_added == true) {
                num_nodes_added += 1;
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("          penalty_new_edge = %lf, min_existing_penalty = %lf\n", edge_scores.edge_score, penalty_max);
                    LOG_ALL("          Added new AnchorGraph edge between: %ld -> %ld, %s\n", pred_anchor->id(), anchor->id(), new_edge->Verbose().c_str());
                );
            }
        #endif
    }

    #ifdef RAPTOR_TESTING_MODE
        DEBUG_RUN(verbose_debug_qid,
            LOG_ALL("          Total added nodes via the popped edge: %ld\n", num_nodes_added)
        );
    #endif
}

void AddLocalEdges_(
                    const raptor::GraphPtr graph,
                    const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors,
                    const std::vector<raptor::AnchorPtr>& flat_anchors,
                    const std::unordered_map<raptor::AnchorPtr, int64_t>& anchor_to_oid,
                    const std::unordered_map<int64_t, int64_t>& tid_to_taid,
                    int32_t max_allowed_dist,
                    int32_t predecessor_lookup_window,
                    int32_t max_path_edges,
                    int32_t chain_max_skip,
                    std::shared_ptr<raptor::AnchorGraph> local_graph,
                    bool verbose_debug_qid) {
    // From every anchor, start a BFS to the furthes reach any walk can have,
    // constrained by the query length and other parameters.
    size_t num_anchors = flat_anchors.size();

    for (size_t anchor_ordinal_id = 0; anchor_ordinal_id < num_anchors; anchor_ordinal_id++) {

        auto& anchor = flat_anchors[anchor_ordinal_id];

        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid,
                LOG_ALL("  Running a BFS from anchor->id() = %ld (anchor_ordinal_id = %ld / %ld).\n", anchor->id(), anchor_ordinal_id, num_anchors)
                LOG_ALL("    Anchor: %s\n", anchor->Verbose().c_str());
            );
        #endif

        // To avoid cycles.
        std::unordered_map<std::shared_ptr<raptor::SegmentEdge>, bool> visited_edges;

        // Keep track of the BFS path inefficiently.
        PathType path;

        // Recursion traversion. Tuple: current_path, jump_length
        std::deque<std::pair<PathType, int64_t>> edge_queue;

        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid,
                LOG_ALL("  Initializing the queue.\n")
            );
        #endif

        // Initialize the queue for the current node.
        // Discovers all the nearby edges and adds them to the queue.
        {
            // I changed the lookup_end and lookup_start here to allow for non 0-length edges.
            int32_t lookup_t_id = anchor->TargetID();
            int32_t lookup_end = anchor->TargetEnd();
            int32_t lookup_start = std::max(0, anchor->TargetStart() - predecessor_lookup_window);

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("    Discovering nearby edges: lookup_t_id = %ld, lookup_start = %ld, lookup_end = %ld\n", lookup_t_id, lookup_start, lookup_end)
                );
            #endif

            std::vector<std::shared_ptr<raptor::SegmentEdge>> avail_edges = graph->FindSegmentInputEdges(lookup_t_id, lookup_start, lookup_end);

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("    avail_edges.size() = %ld\n", avail_edges.size())
                );
            #endif

            for (auto& in_edge: avail_edges) {
                // Take care of the strand.
                if (anchor->TargetRev() != in_edge->sink_is_rev()) {
                    continue;
                }

                // Tuple: (path, jump_len). The `jump_len` is the distance covered
                // by hopping through SegmentEdges (i.e. sum of the length of the
                // sequences in between any two SegmentEdges).
                edge_queue.emplace_back(std::make_pair(PathType{in_edge}, 0));

                #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                        LOG_ALL("      - Adding edge '%s' to q.\n", in_edge->Verbose().c_str())
                    );
                #endif
            }
        }

        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid,
                LOG_ALL("    edge_queue.size() = %ld\n", edge_queue.size());
                LOG_ALL("    Running BFS.\n")
            );
        #endif

        /////////////////////////////////////////////////
        /// Everything up to here sounds good.
        /// Work needed below.
        /////////////////////////////////////////////////

        // BFS for the node. Using a BFS to limit the search depth if user specified.
        while (edge_queue.size() > 0) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("    while(%ld > 0)\n", edge_queue.size())
                );
            #endif

            // Get the path and the jump_len from the queue, and pop it.
            auto& edge_path_pair = edge_queue.front();
            auto path = std::get<0>(edge_path_pair);
            int64_t jump_len = std::get<1>(edge_path_pair);
            edge_queue.pop_front();           // Imitate stack.

            auto edge = path.back();
            visited_edges[edge] = true;

            // Define some values here for easier usage.
            int32_t source_t_id = edge->source_id();
            int32_t source_lookup_end = edge->source_end() + 0;
            int32_t source_lookup_start = std::max((int64_t) 0, (int64_t) (edge->source_start() - predecessor_lookup_window));
            bool source_t_is_rev = edge->source_is_rev();

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid,
                    LOG_ALL("      Popped edge: '%s'\n", edge->Verbose().c_str());
                    LOG_ALL("        Discovering nodes in range: source_t_id = %ld, source_lookup_start = %ld, source_lookup_end = %ld\n", source_t_id, source_lookup_start, source_lookup_end)
                );
            #endif

            // Discover all nodes (anchors) withing edge's range and add edges to them.
            DiscoverPredecessorNodes_(target_anchors, anchor_to_oid, tid_to_taid,
                                      path, source_t_id, source_lookup_start,
                                      source_lookup_end, source_t_is_rev, anchor,
                                      max_allowed_dist, chain_max_skip, local_graph,
                                      verbose_debug_qid);

            // Limit the number of successive edges if specified by the user, to ensure linear runtime.
            if (max_path_edges > 0 && path.size() >= max_path_edges) {
                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid,
                        LOG_ALL("    Stopping to extend this particular path further because path.size() = %ld >= max_path_edges = %ld. "
                                "This means, graph->FindSegmentInputEdges will not be run anymore to add predecessor edges.\n", path.size(), max_path_edges)
                    );
                #endif
                // continue;

            } else {
                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid,
                        LOG_ALL("        Discovering input segment edges in range: source_t_id = %ld, source_lookup_start = %ld, source_lookup_end = %ld\n",
                                            source_t_id, source_lookup_start, source_lookup_end)
                    );

                    DEBUG_RUN(verbose_debug_qid,
                            LOG_ALL("        Visited edges: ");
                            for (const auto& ve: visited_edges) {
                                LOG_NOHEADER("'%s', ", ve.first->symbolic_edge_name().c_str());
                            }
                            LOG_NOHEADER("\n");
                    );
                #endif

                // Add all edge's predecessor edges to the queue.
                auto avail_in_edges = graph->FindSegmentInputEdges(source_t_id, source_lookup_start, source_lookup_end);
                for (auto& in_edge: avail_in_edges) {
                    // Skip visited edges.
                    auto it_visited = visited_edges.find(in_edge);
                    if (it_visited != visited_edges.end() && it_visited->second == true) {
                        #ifdef RAPTOR_TESTING_MODE
                            DEBUG_RUN(verbose_debug_qid,
                                LOG_ALL("          - Continue (1) because edge '%s' is already on the visited list: visited_edges.find(in_edge) != visited_edges.end().\n", in_edge->symbolic_edge_name().c_str())
                            );
                        #endif
                        continue;
                    }
                    // Keep in mind the strand.
                    if (edge->source_is_rev() != in_edge->sink_is_rev()) {
                        #ifdef RAPTOR_TESTING_MODE
                            DEBUG_RUN(verbose_debug_qid,
                                LOG_ALL("          - Continue (2) because edge '%s' is on the wrong strand of the segment edge->source_is_rev() != in_edge->sink_is_rev().\n", in_edge->symbolic_edge_name().c_str())
                            );
                        #endif
                        continue;
                    }

                    // Calculate distance between edges. If they overlap, mark the distance as 0.
                    int64_t edge_dist = (in_edge->sink_end() < edge->source_start()) ? (edge->source_start() - in_edge->sink_end()) : 0;
                    int64_t new_jump_len = jump_len + edge_dist;

                    // Don't extend the queue if we are already too far.
                    if (new_jump_len >= max_allowed_dist) {
                        #ifdef RAPTOR_TESTING_MODE
                            DEBUG_RUN(verbose_debug_qid,
                                LOG_ALL("          - Continue (3) before adding edge '%s' because new_jump_len = %ld >= max_allowed_dist = %ld.\n", in_edge->symbolic_edge_name().c_str(), new_jump_len, max_allowed_dist)
                            );
                        #endif
                        continue;
                    }

                    #ifdef RAPTOR_TESTING_MODE
                        DEBUG_RUN(verbose_debug_qid,
                            LOG_ALL("          - Extending the path with the new edge: %s\n", in_edge->Verbose().c_str());
                        );
                    #endif

                    // Extend the path. Expensive operation.
                    PathType new_path(path);
                    new_path.emplace_back(in_edge);
                    // Add to queue.
                    edge_queue.emplace_back(std::make_pair(new_path, new_jump_len));
                }
            }

            visited_edges[edge] = false;
        }

    }

}

std::shared_ptr<raptor::AnchorGraph> GraphMapper::ConstructAnchorGraph(
                                                    const raptor::GraphPtr graph,
                                                    const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors,
                                                    int32_t max_allowed_dist,
                                                    int32_t predecessor_lookup_window,
                                                    int32_t max_path_edges,
                                                    int32_t chain_max_skip,
                                                    bool verbose_debug_qid) {

    // Create an empty graph.
    auto local_graph = raptor::AnchorGraph::createGraph();

    // Moram napraviti kopiju svih target anchora, da ih mogu sortirati. Edit: stavio sam sortiranje u Mapper, tako da ovdje mogu samo provjeriti da li je sve sortirano, i failati
    // ako nije.
    // Ako je sortirano, onda mozemo usparati puno vremena na kopiranje, i na binary search i sve.

    // Za svaki anchor zelim uzeti prethodnih 50. To je jednostavno ako imam implicitne edgeve, ali sto ako imam explicitne?
    // Umjesto IntervalTree-ja za koji ne znam tocno kako je order kojim mi vraca, da koristim obican sort i binary search?


    // std::vector<std::shared_ptr<raptor::TargetAnchorType>> target_anchors;
    // auto target_anchors = input_mapping_result->target_anchors();

    // Create a flat representation of anchors, for easier traversal.
    std::vector<raptor::AnchorPtr> flat_anchors;
    for (auto& single_target_anchors: target_anchors) {
        for (size_t i = 0; i < single_target_anchors->hits().size(); i++) {
            auto& anchor = single_target_anchors->hits()[i];
            flat_anchors.emplace_back(anchor);

            #ifdef RAPTOR_TESTING_MODE
                if (verbose_debug_qid) {
                    LOG_ALL("i = %ld, id = %ld, anchor = %s\n", i, anchor->id(), anchor->Verbose().c_str());
                }
            #endif
        }
    }

    std::sort(flat_anchors.begin(), flat_anchors.end(),
        [](const raptor::AnchorPtr& a, const raptor::AnchorPtr& b) {
                return (a->TargetID() < b->TargetID()) || (a->TargetID() == b->TargetID() && (a->TargetRev() < b->TargetRev())) ||
                (a->TargetID() == b->TargetID() && a->TargetRev() == b->TargetRev() && (a->TargetStart() < b->TargetStart())); } );

    // Hash anchor pointers to "ordinal ID" of the anchor in the vector.
    // This will be the ID used in DP.
    std::unordered_map<raptor::AnchorPtr, int64_t> anchor_to_oid;
    for (int64_t i = 0; i < (int64_t) flat_anchors.size(); i++) {
        anchor_to_oid[flat_anchors[i]] = i;
    }

    // Create a lookup of target ID to the correct element of the target_anchors
    // vector, for easier access when discovering the predecessors.
    std::unordered_map<int64_t, int64_t> tid_to_taid;
    for (int64_t i = 0; i < (int64_t) target_anchors.size(); i++) {
        int64_t key = target_anchors[i]->env()->t_id << 1;
        key |= ((int64_t) (target_anchors[i]->env()->t_rev) ? 1 : 0);
        tid_to_taid[key] = i;
        // printf ("target_anchors", target_anchors[i]->env()->
        // printf ("key = %ld, target_anchors[i]->env()->t_id = %ld, target_anchors[i]->env()->t_rev = %d, i = %ld\n", key, target_anchors[i]->env()->t_id, (int) target_anchors[i]->env()->t_rev, i);
    }

    // Initialize graph nodes.
    for (size_t anchor_ordinal_id = 0; anchor_ordinal_id < flat_anchors.size(); anchor_ordinal_id++) {
        auto& anchor = flat_anchors[anchor_ordinal_id];
        local_graph->AddNode(anchor->id(), raptor::createAnchorGraphNode(anchor));

        // LOG_ALL("anchor_ordinal_id = %ld, id = %ld, anchor = %s\n", anchor_ordinal_id, anchor->id(), anchor->Verbose().c_str());
    }

    // As the name says, add implicit edges.
    AddImplicitEdges_(target_anchors, anchor_to_oid, max_allowed_dist, chain_max_skip, local_graph, verbose_debug_qid);

    // Runs a BFS over all reachable SegmentEdges, and for each edge in the path
    // discovers nodes. A series of SegmentEdges required to traverse from node
    // N1 to node N2 is called a LocalEdge.
    AddLocalEdges_(graph, target_anchors, flat_anchors, anchor_to_oid, tid_to_taid, max_allowed_dist, predecessor_lookup_window, max_path_edges, chain_max_skip, local_graph, verbose_debug_qid);

    return local_graph;
}

void GraphMapper::AlignAnchors_(const mindex::SequencePtr& qseq,
                                std::vector<raptor::TargetAnchorPtr>& target_anchors) const {
    for (auto& ta : target_anchors) {
        for (auto& anchor : ta->hits()) {
            int32_t matches = raptor::mapper::CalcMatchRate(qseq, index_, anchor);
            anchor->score(matches);
        }
    }
}

std::shared_ptr<raptor::AnchorGraph> GraphMapper::GraphDP2(
                std::shared_ptr<raptor::AnchorGraph>& anchor_graph,
                int32_t chain_max_bandwidth,
                bool verbose_debug_qid) {

    // Stores the final results, a new graph with scores on nodes and edges.
    // Edges can be removed here, marking the unlikely branching.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = anchor_graph;

    // Sort the graph.
    std::vector<int64_t> sorted_nodes;
    bool is_dag = mapped_graph->TopologicalSort(sorted_nodes);
    if (is_dag == false) {
        LOG_ALL("Warning: Graph is cyclic, sorted_nodes will be 0 in size.\n");
    }
    size_t num_nodes = sorted_nodes.size();
    size_t num_edges = mapped_graph->CountEdges();

    #ifdef RAPTOR_TESTING_MODE
        if (verbose_debug_qid) {
            LOG_ALL("Topologically sorted nodes:\n");
            for (auto& node: sorted_nodes) {
                LOG_NOHEADER("%ld ", node);
            }
            LOG_NOHEADER("\n");
        }
    #endif

    // I need to reset the scores here.
    for (auto& v_item: mapped_graph->nodes()) {
        v_item->data()->score(-PlusInf64);
        v_item->data()->path_id(-1);
        v_item->data()->predecessor(-1);
    }
    for (auto& e_item: mapped_graph->edges()) {
        e_item->data()->score(-PlusInf64);
        e_item->data()->path_id(-1);
    }

    // Nothing to do after this.
    if (sorted_nodes.size() == 0) {
        return mapped_graph;
    }

    int64_t num_paths = 0;

    // Fill out the DP matrix.
    for (int64_t i = 0; i < (int64_t) sorted_nodes.size(); i++) {
        const int64_t w_name = sorted_nodes[i];
        auto w_item = mapped_graph->GetNode(w_name);
        const auto& w_anchor = w_item->anchor();
        const int32_t w_score = (w_anchor->Score() > 0) ? w_anchor->Score() : w_anchor->QuerySpan();

        int64_t new_dp_val = w_score;
        int64_t new_dp_pred = -1;
        int64_t new_dp_path = num_paths;
        std::shared_ptr<raptor::AnchorGraphEdge> new_dp_pred_edge = nullptr;

        #ifdef RAPTOR_TESTING_MODE
            if (verbose_debug_qid) {
                    LOG_NOHEADER("\n");
                    LOG_NOHEADER("[i = %ld, w_name = %ld]\n", i, w_name);
                    LOG_NOHEADER("w = '%s'\n", w_anchor->Verbose().c_str());
                    LOG_NOHEADER("-> Initial dp value for i = %ld: new_dp_val = %ld, new_dp_pred = %ld, new_dp_path = %ld\n\n", i, new_dp_val, new_dp_pred, new_dp_path);
            }
        #endif

        for (auto& e_item: mapped_graph->GetInEdges(w_name)) {
            const int64_t v_name = e_item->source_name();
            const auto& v_item = mapped_graph->GetNode(v_name);
            const auto& v_anchor = v_item->anchor();
            const int32_t v_score = (v_anchor->Score() > 0) ? v_anchor->Score() : v_anchor->QuerySpan();

            #ifdef RAPTOR_TESTING_MODE
                if (verbose_debug_qid) {
                    LOG_NOHEADER("\t[i = %ld, v_name = %ld, w_name = %ld, v = '%s']\n", i, v_name, w_name, v_anchor->Verbose().c_str());
                    LOG_NOHEADER("\t\t");
                }
            #endif

            auto es_vals = raptor::graphmapper::CalcAnchorGraphEdgePenalty(v_anchor, w_anchor, e_item->data(), v_score, w_score, verbose_debug_qid);

            // Sanity check. Verify that the edge satisfies
            // rational coordinate restrictions.
            if (es_vals.status_flag != 0) {
                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Continue. es_vals.status_flag = %ld\n", es_vals.status_flag));
                #endif
                continue;
            }
            // This is not a sanity check, but actually necessary to compensate for the wormhole.
            if (es_vals.gap_dist() > chain_max_bandwidth) {
                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Continue. es_vals.gap_dist() = %ld, chain_max_bandwidth = %ld\n", es_vals.gap_dist(), chain_max_bandwidth));
                #endif
                continue;
            }

            int64_t score_vw = v_item->score() + es_vals.score_vw;
            e_item->data()->score(score_vw);

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- v_score = %ld, w_score = %ld, score_vw = %ld\n", v_score, w_score, score_vw););
            #endif

            if (score_vw > new_dp_val) {
                new_dp_pred = v_name;
                new_dp_val = score_vw;
                new_dp_path = v_item->path_id();
                new_dp_pred_edge = e_item->data();

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t=> New score is larger! score_vw = %ld, new_dp_val = %ld, new_dp_path = %ld, new_dp_pred = %ld\n", score_vw, new_dp_val, new_dp_path, new_dp_pred));
                #endif
            }

            #ifdef RAPTOR_TESTING_MODE
                if (verbose_debug_qid) {
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\n"));
                }
            #endif

        }           // for (auto& e: local_graph->GetInEdges(node_i_name))
        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\n"));
        #endif

        w_item->score(new_dp_val);
        w_item->predecessor(new_dp_pred);
        w_item->path_id(new_dp_path);
        if (new_dp_pred_edge != nullptr) {
            new_dp_pred_edge->path_id(new_dp_path);
        }

        // Counts the number of chains that were found.
        if (new_dp_path == num_paths) {
            num_paths += 1;
        }
    }

    for (auto& e_item: mapped_graph->edges()) {
        if (e_item->is_removed()) {
            continue;
        }
        auto v_data = mapped_graph->GetNode(e_item->source_name());
        auto w_data = mapped_graph->GetNode(e_item->sink_name());
        if (v_data->path_id() == w_data->path_id()) {
            e_item->data()->path_id(v_data->path_id());
        } else {
            e_item->data()->path_id(-1);
        }
    }

    return mapped_graph;
}

// void GraphMapper::ReduceMappedAnchorGraph(
//                                             std::shared_ptr<raptor::AnchorGraph>& graph,
//                                             double edge_retain_frac,    // For a node with multiple out edges, pick all those with score >= (max_score * edge_retain_frac).
//                                             bool verbose_debug_qid) {

//     std::set<std::pair<int64_t, int64_t>> edges_for_removal;

//     // Find all edges which are not part of the maximum scoring path, and
//     // which do not have scores _close_ to the maximum (within the edge_retain_frac from max).
//     for (const auto& v_item: graph->nodes()) {

//         int64_t v = v_item->name();
//         if (v_item->is_removed()) {
//             continue;
//         }
//         const auto& v_data = v_item->data();

//         auto out_edges = graph->GetOutEdges(v);
//         auto in_edges = graph->GetInEdges(v);

//         // Filter out edges. Retain only the ones within allowed fraction
//         // from maximum.
//         {
//             // Find the maximum score.
//             double max_out_score = -static_cast<double>(PlusInf64);
//             for (const auto& e_item: out_edges) {
//                 if (e_item->is_removed()) {
//                     continue;
//                 }
//                 // int64_t w = e_item->sink_name();
//                 // const auto& w_data = graph->GetNode(w);
//                 // double score = e_item->data()->score() + v_data->score() + w_data->score();
//                 double score = e_item->data()->score();
//                 max_out_score = std::max(max_out_score, score);
//             }
//             // Take care of both positive and negative scores.
//             double min_out_score = max_out_score - (std::abs(max_out_score) * edge_retain_frac);
//             // Do the filtering.
//             for (const auto& e_item: out_edges) {
//                 if (e_item->is_removed()) {
//                     continue;
//                 }
//                 int64_t w = e_item->sink_name();
//                 // const auto& w_data = graph->GetNode(w);
//                 // double score = e_item->data()->score() + v_data->score() + w_data->score();
//                 double score = e_item->data()->score();
//                 if (score < min_out_score) {
//                     edges_for_removal.emplace(std::make_pair(v, w));
//                 }
//             }
//         }

//         {
//             // Find the maximum score.
//             double max_in_score = -static_cast<double>(PlusInf64);
//             for (const auto& e_item: in_edges) {
//                 if (e_item->is_removed()) {
//                     continue;
//                 }
//                 // int64_t vv = e_item->source_name();
//                 // const auto& vv_data = graph->GetNode(vv);
//                 // double score = e_item->data()->score() + vv_data->score() + v_data->score();
//                 double score = e_item->data()->score();
//                 max_in_score = std::max(max_in_score, score);
//             }
//             // Take care of both positive and negative scores.
//             double min_in_score = max_in_score - (std::abs(max_in_score) * edge_retain_frac);
//             // Do the filtering.
//             for (const auto& e_item: in_edges) {
//                 if (e_item->is_removed()) {
//                     continue;
//                 }
//                 double score = e_item->data()->score();
//                 int64_t vv = e_item->source_name();
//                 // const auto& vv_data = graph->GetNode(vv);
//                 // double score = e_item->data()->score() + vv_data->score() + v_data->score();
//                 if (score < min_in_score) {
//                     edges_for_removal.emplace(std::make_pair(vv, v));
//                 }
//             }
//         }
//     }

//     // Clean up the edges which do not belong to any path.
//     for (const auto& e_item: graph->edges()) {
//         if (e_item->is_removed()) {
//             continue;
//         }
//         const auto& e_data = e_item->data();

//         // Remove edges that don't belong to any path.
//         if (e_data->path_id() < 0) {
//             edges_for_removal.emplace(std::make_pair(e_item->source_name(), e_item->sink_name()));
//             continue;
//         }

//         // Safety guard, this should not happen, but still.
//         // Remove edges between different paths.
//         const auto& v_data = graph->GetNode(e_item->source_name());
//         const auto& w_data = graph->GetNode(e_item->sink_name());
//         if (v_data->path_id() != w_data->path_id()) {
//             edges_for_removal.emplace(std::make_pair(e_item->source_name(), e_item->sink_name()));
//             continue;
//         }
//     }

//     // Remove all the edges from the graph.
//     for (const auto& e_pair: edges_for_removal) {
//         graph->RemoveEdges(std::get<0>(e_pair), std::get<1>(e_pair));
//     }
// }

std::vector<std::shared_ptr<raptor::LocalPath>> GraphMapper::GraphToPaths2_(
                                        const std::shared_ptr<raptor::AnchorGraph>& graph,
                                        const std::unordered_map<int64_t, int64_t>& node_to_path_score,
                                        const bool verbose_debug_qid) {

    std::vector<std::shared_ptr<raptor::LocalPath>> ret;

    #ifdef RAPTOR_TESTING_MODE
        DEBUG_RUN(verbose_debug_qid, LOG_ALL("Backtrack graph:\n%s\n\n", graph->Verbose().c_str()));
    #endif

    // Find source nodes. There can be many paths in a graph.
    std::vector<int64_t> sources;
    for (auto& node_item: graph->nodes()) {
        if (node_item->is_removed()) {
            continue;
        }
        int64_t node_name = node_item->name();
        auto edge_ptr_vector = graph->GetInEdges(node_name);
        if (edge_ptr_vector.size() == 0) {
            sources.push_back(node_name);
        }
    }

    for (const auto& source_name: sources) {
        int64_t v_name = source_name;
        auto v_data = graph->GetNode(v_name);
        auto v_anchor = v_data->anchor();
        auto out_edges = graph->GetOutEdges(v_name);

        int64_t leaf_score = v_data->score();
        int64_t max_path_score = leaf_score;
        int64_t max_path_node = v_name;

        // Find the maximum node in the path.
        while (out_edges.size() == 1) {
            const auto& edge_item = out_edges[0];
            auto w_name = edge_item->sink_name();
            auto w_data = graph->GetNode(w_name);
            auto w_anchor = w_data->anchor();
            int64_t w_score = w_data->score();
            if (w_score > max_path_score) {
                max_path_score = w_score;
                max_path_node = w_name;
            }
            out_edges = graph->GetOutEdges(w_name);
        }

        // Extract the path up to the maximum node.
        int64_t w_name = v_name;
        auto w_data = v_data;
        auto new_path = raptor::LocalPath::createPath(w_name, w_data->anchor());
        out_edges = graph->GetOutEdges(w_name);
        while (out_edges.size() == 1 && w_name != max_path_node) {
            auto& edge_item = out_edges[0];
            w_name = edge_item->sink_name();
            w_data = graph->GetNode(w_name);
            auto w_anchor = w_data->anchor();
            // path_score = w_data->score();
            new_path->Add(w_name, w_anchor, edge_item->data());
            out_edges = graph->GetOutEdges(w_name);
        }

        auto it_score = node_to_path_score.find(v_name);
        int64_t path_score = 0;
        if (it_score == node_to_path_score.end()) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_ALL("Node not found! v_name = %ld\n", v_name));
            #endif
        } else {
            path_score = it_score->second;
        }
        new_path->score(path_score);

        ret.emplace_back(new_path);
    }

    return ret;
}

std::vector<raptor::AnchorGraphPtr> GraphMapper::BacktrackReducedMappedAnchorGraph(
                                            const std::shared_ptr<raptor::AnchorGraph>& graph,
                                            const double edge_retain_frac,
                                            std::unordered_map<int64_t, int64_t>& node_to_path_score,    // For a node with multiple out edges, pick all those with score >= (max_score * edge_retain_frac).
                                            const bool verbose_debug_qid) {

    std::vector<raptor::AnchorGraphPtr> ret_graphs;

    // Find out the maximum number of paths in the graph. Needed to preallocate space for leafs.
    int32_t num_paths = 0;
    for (const auto& v_item: graph->nodes()) {
        if (v_item->is_removed()) {
            continue;
        }
        num_paths = std::max(num_paths, v_item->data()->path_id() + 1);
    }

    // std::unordered_map<int64_t, int64_t> node_to_path_score;
    node_to_path_score.clear();

    // Find all leafs, and find the maximum scoring leaf.
    std::vector<std::vector<int64_t>> pid_leafs(num_paths, std::vector<int64_t>());
    std::vector<int64_t> max_pid_leaf_nodes(num_paths, -PlusInf64);
    std::vector<int64_t> max_pid_leaf_scores(num_paths, -PlusInf64);
    for (const auto& v_item: graph->nodes()) {
        int64_t v = v_item->name();
        if (v_item->is_removed()) {
            continue;
        }
        const auto& v_data = v_item->data();
        int32_t pid = v_data->path_id();

        // All nodes should have a non-negative path ID, but let's be safe.
        if (pid < 0) {
            continue;
        }

        auto out_edges = graph->GetOutEdges(v);
        int64_t num_valid_out_edges = 0;
        for (const auto& e_item: out_edges) {
            if (e_item->data()->path_id() == -1) {
                continue;
            }
            if (e_item->data()->path_id() != pid) {
                continue;
            }
            ++num_valid_out_edges;
        }

        // Take only sources.
        if (num_valid_out_edges == 0) {
            pid_leafs[pid].emplace_back(v);
            node_to_path_score[v] = v_data->score();
            if (v_data->score() > max_pid_leaf_scores[pid]) {
                max_pid_leaf_scores[pid] = v_data->score();
                max_pid_leaf_nodes[pid] = v;
            }
        }
    }

    #ifdef RAPTOR_TESTING_MODE
        if (verbose_debug_qid) {
            LOG_NEWLINE;
            LOG_ALL("Leafs for each path:\n");
            for (size_t pid = 0; pid < pid_leafs.size(); pid++) {
                LOG_NOHEADER("  - [path_id = %ld]: ", pid);
                auto& path_dp_ids = pid_leafs[pid];
                for (size_t j = 0; j < path_dp_ids.size(); j++) {
                    LOG_NOHEADER("%ld ", path_dp_ids[j]);
                }
                LOG_NOHEADER("\n");
                LOG_NOHEADER("      - max_pid_leaf_node = %ld, max_pid_leaf_score = %ld\n", max_pid_leaf_nodes[pid], max_pid_leaf_scores[pid]);
            }
            LOG_NOHEADER("\n");

        }
    #endif

    for (size_t pid = 0; pid < pid_leafs.size(); pid++) {
        std::unordered_map<int64_t, int8_t> visited_nodes;

        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid, LOG_ALL("Processing path %ld\n", pid));
        #endif

        // Skip empty components.
        if (pid_leafs[pid].size() == 0) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_ALL("No leaf nodes for path pid = %ld. Skipping.\n", pid));
            #endif
            continue;
        }

        // Add the maximum leaf node to the queue. Backtrack will start from there.
        // The path score will be propagated with the queue, so that alternative bubbles
        // will appear in the output.
        std::deque<std::pair<int64_t, int64_t>> fork_queue;

        // Emplace all leafs. We expect that the filtering was applied on the outside.
        double min_leaf_score = max_pid_leaf_scores[pid] - std::abs(max_pid_leaf_scores[pid] * edge_retain_frac);
        for (const auto& v: pid_leafs[pid]) {
            const auto& v_data = graph->GetNode(v);
            if (v_data->score() < min_leaf_score) {
                continue;
            }
            fork_queue.emplace_back(std::make_pair(v, v_data->score()));
            node_to_path_score[v] = v_data->score();
        }

        // Create a graph, where supplementary mappings will be individual
        // connected components, where each component is a linear path.
        auto backtrack_graph = raptor::AnchorGraph::createGraph();

        // The queue contains the beginnings of all similarly good subpaths.
        while (fork_queue.size() > 0) {

            // Pop the queue to get the first node.
            std::pair<int64_t, int64_t> name_score_pair = fork_queue.front();
            int64_t v_name = std::get<0>(name_score_pair);
            int64_t leaf_score = std::get<1>(name_score_pair);  // This score will be propagated to all successive nodes recursively.
            auto v_data = graph->GetNode(v_name);
            int64_t score = v_data->score();
            fork_queue.pop_front();

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\tPopped v_name = %ld, v_score = %ld\n", v_name, score));
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\tFetched the first node data.\n"));
            #endif

            if (visited_nodes.find(v_name) != visited_nodes.end()) {
                continue;
            }

            int64_t w_name = v_name;
            auto w_data = v_data;
            backtrack_graph->AddNode(w_name, w_data);
            node_to_path_score[w_name] = leaf_score;

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\tStarting with: curr_node_name = %ld\n", w_name));
            #endif

            // Backtrack the best path, and add the nodes and edges to the backtrack_graph.
            while (true) {
                // Mark the node as visited.
                visited_nodes[w_name] = 1;

                // Find the best incomming edges.
                std::vector<raptor::EdgeItemAnchorGraphPtr>in_edges = graph->GetInEdges(w_name);

                // Find the maximum in-edge.
                raptor::EdgeItemAnchorGraphPtr e_in_max = nullptr;
                double max_in_score = -static_cast<double>(PlusInf64);
                for (const auto& e_item: in_edges) {
                    if (e_item->is_removed()) {
                        // Skip if edge is marked as removed.
                        continue;
                    }
                    if (e_item->data()->path_id() != pid) {
                        // Skip if the edge is not classified as on this path.
                        continue;
                    }
                    if (graph->GetNode(e_item->source_name())->path_id() != pid) {
                        // Skip if the source is not of the same path.
                        continue;
                    }
                    if (visited_nodes.find(e_item->source_name()) != visited_nodes.end()) {
                        // Skip if visited.
                        continue;
                    }
                    // Finally, check max.
                    double score = e_item->data()->score();
                    if (e_in_max == nullptr || e_item->data()->score() > max_in_score) {
                        e_in_max = e_item;
                        max_in_score = e_item->data()->score();
                    }
                }

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - max_in_score = %lf\n", max_in_score));
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - in_edges.size() = %ld\n", in_edges.size()));
                #endif

                if (e_in_max == nullptr) {
                    // There are either no input edges, or all are visited.
                    #ifdef RAPTOR_TESTING_MODE
                        DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - Break: e_in_max == nullptr\n"));
                    #endif
                    break;
                }

                double min_allowed_in_edge_score = max_in_score * edge_retain_frac;

                // Reiterate over all edges, and add all alternate paths to the queue.
                // This expects that all edges which are not of interest are already
                // filtered prior to entering this function.
                for (const auto& e_item: in_edges) {
                    if (e_item->is_removed()) {
                        // Skip if edge is marked as removed.
                        continue;
                    }
                    if (e_item->data()->path_id() != pid) {
                        // Skip if the edge is not classified as on this path.
                        continue;
                    }
                    if (graph->GetNode(e_item->source_name())->path_id() != pid) {
                        // Skip if the source is not of the same path.
                        continue;
                    }
                    if (visited_nodes.find(e_item->source_name()) != visited_nodes.end()) {
                        // Skip if visited.
                        continue;
                    }
                    if (e_item == e_in_max) {
                        // Skip the max edge.
                        continue;
                    }
                    if (e_item->data()->score() < min_allowed_in_edge_score) {
                        // Skip poor edges.
                        continue;
                    }
                    // Finally, add to the queue.
                    fork_queue.push_back(std::make_pair(e_item->source_name(), leaf_score));
                }

                w_name = e_in_max->source_name();
                w_data = graph->GetNode(w_name);

                backtrack_graph->AddNode(w_name, w_data);
                backtrack_graph->AddEdge(e_in_max->source_name(), e_in_max->sink_name(), e_in_max->data());

                // Set the new node to have a score of the leaf.
                // backtrack_graph->GetNode(w_name)->score(leaf_score);
                node_to_path_score[w_name] = leaf_score;

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\tMoving to: w_name = %ld\n", w_name));
                #endif
            }
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\tbacktrack_graph.nodes().size() = %d\n", backtrack_graph->nodes().size()));
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\n"));
            #endif
        }

        ret_graphs.emplace_back(backtrack_graph);
    }

    return ret_graphs;
}

raptor::EdgeItemAnchorGraphPtr FindMaxInEdgeForPathId(
                const std::shared_ptr<raptor::AnchorGraph>& graph,
                const std::vector<raptor::EdgeItemAnchorGraphPtr>& in_edges,
                const std::unordered_map<int64_t, int8_t>& visited_edges,
                int64_t pid
                ) {
    raptor::EdgeItemAnchorGraphPtr e_in_max = nullptr;
    double max_in_score = -static_cast<double>(PlusInf64);
    for (const auto& e_item: in_edges) {
        if (e_item->is_removed()) {
            // Skip if edge is marked as removed.
            continue;
        }
        if (e_item->data()->path_id() != pid) {
            // Skip if the edge is not classified as on this path.
            continue;
        }
        if (graph->GetNode(e_item->source_name())->path_id() != pid) {
            // Skip if the source is not of the same path.
            continue;
        }
        if (visited_edges.find(e_item->name()) != visited_edges.end()) {
            // Skip if visited.
            continue;
        }
        // Finally, check max.
        double score = e_item->data()->score();
        if (e_in_max == nullptr || e_item->data()->score() > max_in_score) {
            e_in_max = e_item;
            max_in_score = e_item->data()->score();
        }
    }
    return e_in_max;
}

/*
 * The CreateAnchorGraphConnectedComponents method separates connected components in the graph.
 * It prunes the graph by selecting the highest scoring paths, and filtering the rest.
 *
 * For each path ID component, this method first collects the leaf nodes. Leafs aren't in any
 * particular sorted order. The maximum scoring leaf is found, and used as a representative
 * score for this graph component.
 * For every leaf node in the queue, we begin backtracking along the highest scoring path,
 * while also keeping edges and traversing alternate paths to within a margin from the
 * top score (allowed_score_diff_frac).
 * All such reachable nodes and edges are added to the pid graph.
 *
 * At most one implicit edge is allowed, the best scoring one.
*/
std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> GraphMapper::CreateAnchorGraphConnectedComponents(
                                            const std::shared_ptr<raptor::AnchorGraph>& graph,
                                            const double allowed_score_diff_frac,    // For a node with multiple out edges, pick all those with score >= (max_score * edge_retain_frac).
                                            const bool verbose_debug_qid) {

    // Find leafs for all paths.
    std::vector<std::vector<int64_t>> pid_leafs = raptor::graphmapper::FindBestLeafNodes(graph, allowed_score_diff_frac, verbose_debug_qid);

    // Initialize the node scores for the selected leafs.
    std::unordered_map<int64_t, int64_t> node_to_path_score;

    /*
     * Template args: <anchor_graph, path_id, path_score>
    */
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> component_graphs;

    for (size_t pid = 0; pid < pid_leafs.size(); pid++) {
        #ifdef RAPTOR_TESTING_MODE
            DEBUG_RUN(verbose_debug_qid, LOG_ALL("Processing the path pid = %ld\n", pid));
        #endif

        // Collect all leaf nodes and scores. Instead of relying on them being
        // filtered already, prune them here just in case.
        std::vector<std::pair<int64_t, double>> candidate_leafs;
        double max_pid_score = -static_cast<double>(PlusInf64);
        for (const auto& v_name: pid_leafs[pid]) {
            auto v_data = graph->GetNode(v_name);
            candidate_leafs.emplace_back(std::make_pair(v_name, v_data->score()));
            max_pid_score = std::max(max_pid_score, v_data->score());
        }
        // Sort the candidate leafs by scorel
        std::sort(candidate_leafs.begin(), candidate_leafs.end(),
                [](const std::pair<int64_t, double>& a, const std::pair<int64_t, double>& b) {
                    return std::get<1>(a) > std::get<1>(b);
                }
        );
        // Initialize the processing queue with top scoring leafs which are above
        // a heuristic threshold.
        double min_leaf_score = max_pid_score * allowed_score_diff_frac;
        std::deque<std::pair<int64_t, double>> proc_queue;
        for (const auto& leaf_pair: candidate_leafs) {
            if (std::get<1>(leaf_pair) < min_leaf_score) {
                break;
            }
            proc_queue.emplace_back(leaf_pair);
        }

        // Skip empty components.
        if (proc_queue.size() == 0) {
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_ALL("No leaf nodes for path pid = %ld. Skipping.\n", pid));
            #endif
            continue;
        }

        // // Create a graph, where supplementary mappings will be individual
        // // connected components, where each component is a linear path.
        // auto backtrack_graph = raptor::AnchorGraph::createGraph();

        std::unordered_map<int64_t, int8_t> visited_nodes;
        std::unordered_map<int64_t, int8_t> visited_edges;

        // Create an empty graph.
        raptor::AnchorGraphPtr pid_graph = raptor::AnchorGraph::createGraph();

        // All of the nodes in the processing queue will be added to the subgraph component.
        while (proc_queue.size() > 0) {
            // Pop the queue to get the first node.
            std::pair<int64_t, double> name_score_pair = proc_queue.front();
            int64_t v_name = std::get<0>(name_score_pair);
            double leaf_score = std::get<1>(name_score_pair);  // This score will be propagated to all successive nodes recursively.
            auto v_data = graph->GetNode(v_name);
            double score = v_data->score();
            proc_queue.pop_front();

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\tPopped v_name = %ld, v_score = %ld\n", v_name, score));
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\tFetched the first node data.\n"));
            #endif

            // Need to think if this is an ok break condition. I think it is, because
            // otherwise we might make alternate paths which reuse nodes.
            if (visited_nodes.find(v_name) != visited_nodes.end()) {
                continue;
            }

            int64_t w_name = v_name;
            auto w_data = v_data;
            node_to_path_score[w_name] = leaf_score;

            pid_graph->AddNode(w_name, w_data);

            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Starting with: w_name = %ld\n", w_name));
            #endif

            // Backtrack the best path FIRST. The best path follows the maximum scoring edges.
            // Along the way, look for similarly scoring alternate edges, and add the source
            // nodes to the processing queue.
            // This while loop iterates over incoming best scoring path through w_name.
            // Alternative paths are processed through proc_queue, and edges are added here.
            while (true) {
                // Mark the node as visited.
                visited_nodes[w_name] = 1;

                // Find the best incomming edges.
                std::vector<raptor::EdgeItemAnchorGraphPtr> in_edges = graph->GetInEdges(w_name);

                // Find the maximum scoring in-edge.
                raptor::EdgeItemAnchorGraphPtr e_in_max = FindMaxInEdgeForPathId(graph, in_edges, visited_edges, pid);

                // If that one doesn't exist, break. No way we are entering a cycle,
                // because we mark nodes as visited.
                if (e_in_max == nullptr) {
                    // There are either no input edges, or all are visited.
                    #ifdef RAPTOR_TESTING_MODE
                        DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - Break: e_in_max == nullptr\n"));
                    #endif
                    break;
                }

                double max_in_score = e_in_max->data()->score();
                double min_allowed_in_edge_score = max_in_score * allowed_score_diff_frac;

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - max_in_score = %lf\n", max_in_score));
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - in_edges.size() = %ld\n", in_edges.size()));
                #endif

                // Allow only one implicit edge (the best scoring one, within the threshold);
                int64_t selected_implicit_edge = -1;
                std::vector<int64_t> selected_edges;

                // This loop finds all candidate edges, and stores them into the
                // selected_edges vector and the selected_implicit_edge variable
                // which will be appended to the vector.
                for (int64_t e_item_id = 0; e_item_id < static_cast<int64_t>(in_edges.size()); ++e_item_id) {
                    const auto& e_item = in_edges[e_item_id];

                    if (e_item->is_removed()) {
                        // Skip if edge is marked as removed.
                        continue;
                    }
                    if (e_item->data()->path_id() != pid) {
                        // Skip if the edge is not classified as on this path.
                        continue;
                    }
                    if (graph->GetNode(e_item->source_name())->path_id() != pid) {
                        // Skip if the source is not of the same path.
                        continue;
                    }
                    if (visited_edges.find(e_item->name()) != visited_edges.end()) {
                        continue;
                    }
                    if (e_item == e_in_max) {
                        if (e_item->data()->IsImplicit()) {
                            selected_implicit_edge = e_item_id;
                        }
                        // Skip the max edge. We will traverse the max edge
                        // right away.
                        continue;
                    }
                    if (e_item->data()->score() < min_allowed_in_edge_score) {
                        // Skip poor edges.
                        continue;
                    }

                    if (e_item->data()->IsImplicit() == false) {
                        // If the edge is explicit, then just add it.
                        selected_edges.emplace_back(e_item_id);

                    } else if (selected_implicit_edge < 0 ||
                                (selected_implicit_edge >= 0 && e_item->data()->score() > in_edges[selected_implicit_edge]->data()->score())) {
                        // This ensures that only one implicit edge is added.
                        selected_implicit_edge = e_item_id;

                    }
                }

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - selected_implicit_edge = %d\n", selected_implicit_edge));
                    if (selected_implicit_edge >= 0) {
                        DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t  - selected_implicit_edge_data = '%s'\n", in_edges[selected_implicit_edge]->data()->ToJSON().c_str()));
                    }
                    DEBUG_RUN(verbose_debug_qid,
                        LOG_NOHEADER("\t\t  - selected_edges = [");
                        for (size_t ii = 0; ii < selected_edges.size(); ++ii) {
                            if (ii > 0) {
                                LOG_NOHEADER(", ");
                            }
                            LOG_NOHEADER("%ld", selected_edges[ii]);
                        }
                        LOG_NOHEADER("]\n");
                    );
                #endif

                // If there was an implicit edge selected, append it. Otherwise, of course, no.
                // Don't add implicit max edge, because max edges are handled below.
                if (selected_implicit_edge >= 0 && in_edges[selected_implicit_edge] != e_in_max) {
                    selected_edges.emplace_back(selected_implicit_edge);
                }

                // Add all alternate paths to the queue.
                for (const auto& e_item_id: selected_edges) {
                    const auto& e_item = in_edges[e_item_id];

                    // Add the alternate node to the new graph.
                    if (visited_nodes.find(e_item->source_name()) == visited_nodes.end()) {
                        pid_graph->AddNode(e_item->source_name(), graph->GetNode(e_item->source_name()));
                    }
                    pid_graph->AddEdge(e_item->source_name(), e_item->sink_name(), e_item->data());

                    // Finally, add to the queue.
                    proc_queue.push_back(std::make_pair(e_item->source_name(), leaf_score));

                    visited_edges[e_item->name()] = 1;

                    node_to_path_score[e_item->source_name()] = leaf_score;

                    #ifdef RAPTOR_TESTING_MODE
                        DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Added alternate path node: %ld\n", e_item->source_name()));
                        DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Added alternate path edge: %ld -> %ld\n", e_item->source_name(), e_item->sink_name()));
                    #endif
                }

                w_name = e_in_max->source_name();
                w_data = graph->GetNode(w_name);

                pid_graph->AddNode(w_name, w_data);
                pid_graph->AddEdge(e_in_max->source_name(), e_in_max->sink_name(), e_in_max->data());
                visited_edges[e_in_max->name()] = 1;

                // Set the new node to have a score of the leaf.
                node_to_path_score[w_name] = leaf_score;

                #ifdef RAPTOR_TESTING_MODE
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Added main path node: %ld\n", w_name));
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Added main path edge: %ld -> %ld\n", e_in_max->source_name(), e_in_max->sink_name()));
                    DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\t- Moving to: w_name = %ld\n", w_name));
                #endif
            }
            #ifdef RAPTOR_TESTING_MODE
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\t\tbacktrack_graph.nodes().size() = %d\n", pid_graph->nodes().size()));
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\n"));
            #endif
        }

        if (pid_graph->CheckIsEmpty() == false) {
            component_graphs.emplace_back(std::make_tuple(pid_graph, pid, max_pid_score));
        }

        // auto new_paths = GraphToPaths2_(pid_graph, node_to_path_score, verbose_debug_qid);

        // paths.insert(paths.end(), new_paths.begin(), new_paths.end());

    }

    return component_graphs;
}

bool GraphMapper::FindLinearSplitSegmentNodesAndEdges_(const raptor::SplitSegmentGraphPtr& ssg,
                    int64_t t_id, bool t_rev, int64_t start, int64_t end,
                    std::vector<raptor::SplitSegmentGraphNameType>& ret_nodes,
                    std::vector<raptor::SplitSegmentGraphNameType>& ret_edges
                    ) {

    std::vector<raptor::SplitSegmentGraphNameType> ss_node_names =
                        ssg->FindSplitNodesByTargetId(t_id, t_rev, start, end);

    std::vector<std::pair<raptor::SplitSegmentGraphNameType, int64_t>> node_start_pair;
    for (const auto& v_name: ss_node_names) {
        auto v_data = ssg->GetNode(v_name);
        if (v_data == nullptr) {
            return false;
        }
        int64_t v_start = v_data->start();
        node_start_pair.emplace_back(std::make_pair(v_name, v_start));
    }

    std::sort(node_start_pair.begin(), node_start_pair.end(),
                [](const std::pair<raptor::SplitSegmentGraphNameType, int64_t>& a, const std::pair<raptor::SplitSegmentGraphNameType, int64_t>& b) {
                    return std::get<1>(a) < std::get<1>(b);
                });

    ret_nodes.clear();
    ret_edges.clear();

    if (node_start_pair.size() > 0) {
        ret_nodes.emplace_back(std::get<0>(node_start_pair[0]));
    }

    for (size_t i = 1; i < node_start_pair.size(); ++i) {
        auto prev_v_name = std::get<0>(node_start_pair[i-1]);
        auto v_name = std::get<0>(node_start_pair[i]);

        ret_nodes.emplace_back(v_name);

        const auto& edge_items = ssg->GetEdges(prev_v_name, v_name);
        if (edge_items.empty()) {
            return false;
        }
        for (const auto& e_item: edge_items) {
            ret_edges.emplace_back(e_item->name());
        }
    }

    return true;
}

raptor::SplitSegmentGraphPtr GraphMapper::ConstructAlignmentGraph(const raptor::SplitSegmentGraphPtr& ssg,
                                                 const raptor::AnchorGraphPtr& ag,
                                                 const mindex::SequencePtr& qseq,
                                                 bool use_approx_aln,         // TODO: This could require a blacklist of paths which were discarded during GraphDP, e.g. some have mappings but of poorer quality.
                                                 bool verbose_debug_qid) {

    raptor::SplitSegmentGraphPtr aln_graph = raptor::createSplitSegmentGraph();

    const int8_t COLOR_BLANK = 0;       // Uninitialized.
    const int8_t COLOR_TRAVERSED = 1;   // The node/edge was reached, but not in the subgraph.
    const int8_t COLOR_SELECTED = 2;    // In-play, to be used as subgraph.

    LOG_ALL("SSG:\n%s\n", ssg->ToJSON().c_str());
    LOG_ALL("AG:\n%s\n", ag->ToJSON().c_str());

    // std::cerr << "AG:\n" << local_graph->ToJSON() << "\n";
    // exit(1);

    // Make a lookup between AnchorNode -> SplitSegmentNode.
    // There should be 1:1 relation, because GraphMapper breaks anchors on edges.
    std::unordered_map<int64_t, int64_t> agname_to_ssgname;
    for (const auto& node_item: ag->nodes()) {
        if (node_item->is_removed()) {
            LOG_NOHEADER("Skipping node_item for node_name = '%d' when building agname_to_ssgname, it's removed.\n", node_item->name());
            continue;
        }
        const auto& node_data = node_item->data();
        const auto& anchor = node_data->anchor();
        int64_t aname = node_item->name();

        // Take only the start coordinate to detect the SplitSegmentGraph node, because k-mers might span edges even though anchors were broken?
        // TODO: Reconsider this in the future.
        std::vector<raptor::SplitSegmentGraphNameType> ss_node_names =
                            ssg->FindSplitNodesByTargetId(
                                        anchor->TargetID(), anchor->TargetRev(),
                                        anchor->TargetStart(), anchor->TargetStart() + 1);
        // std::vector<raptor::SplitSegmentGraphNameType> ss_node_names = ssg->FindSplitNodesByTargetId(anchor->TargetID(), anchor->TargetRev(), anchor->TargetStart(), anchor->TargetEnd());

        // There should be only 1 SplitSegment node per anchor, otherwise the BreakAnchors function is buggy.
        // Warn in that case.
        if (ss_node_names.size() != 1) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "An anchor matches != 1 SplitSegmentGraph nodes. "
                "This should not happen, likely a bug in BreakAnchors. Returning nullptr for AlignmentGraph.\n"
                "qseq = '%s', anchor_graph_node_name = %ld, ss_node_names.size() = %ld, "
                "anchor_graph_node = '%s'\n",
                qseq->header().c_str(), aname, ss_node_names.size(), node_data->Verbose().c_str());

            #ifdef RAPTOR_TESTING_MODE
                // if (verbose_debug_qid) {
                    LOG_NOHEADER("  - ag_name = %ld, ag_node = '%s'\n",
                        aname, node_data->Verbose().c_str());
                    for (const auto &v_name: ss_node_names) {
                        LOG_NOHEADER("  - ssg_node_name = %ld, ssg_node = '%s'\n",
                            v_name, ssg->GetNode(v_name)->Verbose().c_str());

                    }
                // }
            #endif

            return nullptr;
        }
        // Add the relation.
        agname_to_ssgname[aname] = ss_node_names[0];
        DEBUG_RUN(verbose_debug_qid, LOG_ALL("Adding agname_to_ssgname[%ld] = %ld\n  - ag_node = '%s'\n  - ssg_node = '%s'\n\n",
            aname, ss_node_names[0], node_data->Verbose().c_str(), ssg->GetNode(ss_node_names[0])->Verbose().c_str()));

    }

    for (auto it: agname_to_ssgname) {
        LOG_NOHEADER("agname_to_ssgname[%ld] = %ld\n", it.first, it.second);
    }


    // Color the graph to find a subgraph of interest.
    std::unordered_map<int64_t, int8_t> node_colors;
    std::unordered_map<int64_t, int8_t> edge_colors;

    // Initialize the node and edge colors to unused.
    for (const auto& node_item: ssg->nodes()) {
        if (node_item->is_removed()) {
            continue;
        }
        node_colors[node_item->name()] = COLOR_BLANK;
    }
    for (const auto& edge_item: ssg->edges()) {
        if (edge_item->is_removed()) {
            continue;
        }
        edge_colors[edge_item->name()] = COLOR_BLANK;
    }

    // Set the colors of every SplitSegmentGraph node which contains an anchor to SELECTED.
    for (const auto& it: agname_to_ssgname) {
        node_colors[it.second] = COLOR_SELECTED;
    }

    // Find the exact traversal path of every AnchorGraph edge,
    // and resolve it across underlying split segments and segment edges.
    // Mark all of those nodes and edges as COLOR_SELECTED.
    for (size_t eid = 0; eid < ag->edges().size(); ++eid) {
        if (ag->edges()[eid] == nullptr) {
            continue;
        }

        const raptor::EdgeItemAnchorGraphPtr& edge_item = ag->edges()[eid];

        if (edge_item->is_removed()) {
            continue;
        }

        const raptor::AnchorGraphEdgePtr& edge = edge_item->data();

        int64_t ssg_v = agname_to_ssgname.find(edge_item->source_name())->second;
        int64_t ssg_w = agname_to_ssgname.find(edge_item->sink_name())->second;
        auto ssg_v_data = ssg->GetNode(ssg_v);
        auto ssg_w_data = ssg->GetNode(ssg_w);

        if (ssg_v_data == nullptr || ssg_w_data == nullptr) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Either the ssg_v_data or ssg_w_data is nullptr! ssg_v = %ld, ssg_w = %ld. Returning nullptr for AlignmentGraph.\n", ssg_v, ssg_w);
            return nullptr;
        }

        node_colors[ssg_v] = COLOR_SELECTED;
        node_colors[ssg_w] = COLOR_SELECTED;

        // If an edge is implicit, that's fine, that means a SegmentNode was broken.
        if (edge->IsImplicit()) {
            // An implicit edge in SegmentGraph can span multiple nodes in the
            // SplitSegment graph because the path might have been broken.
            // Here we identify the linear chain of SplitSegmentNodes between
            // the source and sink node of the implicit edge.

            std::vector<raptor::SplitSegmentGraphNameType> nodes_to_color;
            std::vector<raptor::SplitSegmentGraphNameType> edges_to_color;

            bool rv1 = FindLinearSplitSegmentNodesAndEdges_(ssg,
                    ssg_v_data->t_id(), ssg_v_data->is_rev(), ssg_v_data->start(), ssg_w_data->start() + 1, nodes_to_color, edges_to_color);

            for (const auto& node_to_color: nodes_to_color) {
                node_colors[node_to_color] = COLOR_SELECTED;
            }

            for (const auto& edge_to_color: edges_to_color) {
                edge_colors[edge_to_color] = COLOR_SELECTED;
            }

        } else {
            // The following is simple but bloated because of all of the checks.
            // One AnchorGraphEdge can contain 1 or more SegmentEdges, which
            // make a jump path. We need to add all SegmentEdges in the vector,
            // and also discover all SplitSegmentNodes between every pair.

            // Color the SegmentEdges first. Perform checks just in case.
            for (size_t i = 0; i < edge->segment_edges().size(); i++) {
                auto seg_edge = edge->segment_edges()[i];

                if (seg_edge == nullptr) {
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "seg_edge == nullptr! Returning nullptr alignment graph.\n");
                    return nullptr;
                }

                int64_t ssg_source = -1, ssg_sink = -1;
                bool rv_find = ssg->FindSegmentEdgeSourceAndSink(seg_edge, ssg_source, ssg_sink);
                if (rv_find == false) {
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Skipping edge: '%s'. Returning nullptr alignment graph.\n", seg_edge->ToJSON().c_str());
                    return nullptr;
                }

                node_colors[ssg_source] = COLOR_SELECTED;
                node_colors[ssg_sink] = COLOR_SELECTED;

                auto ssg_edges = ssg->GetEdges(ssg_source, ssg_sink);
                if (ssg_edges.size() != 1) {
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "There are %ld edges between SplitSegmentGraph source and sink, even though there should be only 1! Edge: '%s'. Returning nullptr alignment graph.\n", ssg_edges.size(), seg_edge->Verbose().c_str());
                    return nullptr;
                }

                edge_colors[ssg_edges[0]->name()] = COLOR_SELECTED;
            }

            // Add all SplitSegmentNodes in between the source and sink on the same
            // SegmentNode sequence. There could have been edges in between which
            // broke one SegmentNode into multiple SplitSegmentNodes, but none
            // of those edges are part of the final AnchorGraph.
            for (size_t i = 1; i < edge->segment_edges().size(); i++) {
                // We just checked everything, so no need to do it again.
                auto prev_seg_edge = edge->segment_edges()[i-1];
                auto seg_edge = edge->segment_edges()[i];

                int32_t t_id = seg_edge->source_id();
                bool t_rev = seg_edge->source_is_rev();
                auto t_start = prev_seg_edge->sink_start();
                auto t_end = seg_edge->source_start() + 1;

                std::vector<raptor::SplitSegmentGraphNameType> nodes_to_color;
                std::vector<raptor::SplitSegmentGraphNameType> edges_to_color;

                bool rv1 = FindLinearSplitSegmentNodesAndEdges_(ssg,
                        t_id, t_rev, t_start, t_end, nodes_to_color, edges_to_color);

                for (const auto& node_to_color: nodes_to_color) {
                    node_colors[node_to_color] = COLOR_SELECTED;
                }

                for (const auto& edge_to_color: edges_to_color) {
                    edge_colors[edge_to_color] = COLOR_SELECTED;
                }
            }
        }
    }

    if (use_approx_aln == false) {
        // Reuse the same node and edge color map to mark the visited edges and nodes.
        // Reset the non-selected nodes, and collect all selected ones as sources for DFS.
        std::vector<int64_t> source_nodes;
        for (const auto it: node_colors) {
            if (it.second != COLOR_SELECTED) {
                node_colors[it.first] = COLOR_BLANK;
            } else {
                // proc_node_queue.emplace_back(it.first);
                source_nodes.emplace_back(it.first);
            }
        }
        // Reset the non-selected edges.
        for (const auto it: edge_colors) {
            if (it.second != COLOR_SELECTED) {
                edge_colors[it.first] = COLOR_BLANK;
            }
        }

        // For each source node, perform a DFS over non-visited edges.
        // If we hit any selected node, we mark the currently traversed
        // path as COLOR_SELECTED.
        for (const auto& source_name: source_nodes) {
            DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("Going down from source_name = %ld\n", source_name));

            // <sink_node_name, which_edge_was_used_to_reach, dfs_depth>
            std::deque<std::tuple<int64_t, int64_t, int64_t>> dfs_queue = {{source_name, -1, 0}};
            // std::deque<int64_t> dfs_path = {source_name};
            std::deque<std::tuple<int64_t, int64_t, int64_t>> dfs_path = {{source_name, -1, 0}};

            // Initialize the DFS up to depth 1. Otherwise, the break condition below which checks if
            // the node is COLOR_SELECTED will kick us right out.
            auto source_out_e_items = ssg->GetOutEdges(source_name);
            DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("Initializing the DFS. Candidates size = %ld.\n", source_out_e_items.size()));
            for (const auto& e_item: source_out_e_items) {
                auto it_e_color = edge_colors.find(e_item->name());
                // Perform filtering first. If all filters passed, then
                // mark the edge and extend DFS.
                if (it_e_color == edge_colors.end()) {
                    // Sanity check.
                    // All edges should be in the unordered_map, but check just in case.
                    // Something weird happened, and the edge is not in the map.
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Edge name %d was not found in the edge_colors map. Skipping path.", e_item->name());
                    continue;
                } else if (it_e_color->second != COLOR_BLANK) {
                    // Even though the edge color can be "selected", if we
                    // arrived at a node which wasn't colored, but the out edge
                    // is colored, then something is wrong probably.
                    // Somehow, this edge was processed already, and it could
                    // likely be some repeat cycle.
                    continue;
                }
                edge_colors[e_item->name()] = COLOR_TRAVERSED;
                dfs_queue.emplace_back(std::make_tuple(e_item->sink_name(), e_item->name(), 1));
            }

            while (dfs_queue.size()) {
                // Make a copy to be safe after popping.
                const auto back_tuple = dfs_queue.back();
                int64_t v_name = std::get<0>(back_tuple);
                int64_t v_via_edge = std::get<1>(back_tuple);
                int64_t v_depth = std::get<2>(back_tuple);
                dfs_queue.pop_back();

                // Remove all distant visited edges from the path. (Backtrack.)
                while (dfs_path.size() && dfs_path.size() >= v_depth) {
                    dfs_path.pop_back();
                }
                // Add the currently popped node to the visited list to track the path.
                dfs_path.emplace_back(back_tuple);

                DEBUG_RUN(verbose_debug_qid,
                    LOG_NOHEADER("popped: (v = %ld, via_e = %ld, depth = %ld)\n", v_name, v_via_edge, v_depth);
                    LOG_NOHEADER("  - dfs_queue (size = %ld):\n", dfs_queue.size());
                    size_t temp_i1 = 0;
                    for (const auto& vals: dfs_queue) {
                        LOG_NOHEADER("      [%ld] (v = %ld, via_e = %ld, depth = %ld)\n", temp_i1, std::get<0>(vals), std::get<1>(vals), std::get<2>(vals));
                        ++temp_i1;
                    }
                    LOG_NOHEADER("  - dfs_path (size = %ld):\n", dfs_path.size());
                    size_t temp_i2 = 0;
                    for (const auto& vals: dfs_path) {
                        LOG_NOHEADER("      [%ld] (v = %ld, via_e = %ld, depth = %ld)\n", temp_i2, std::get<0>(vals), std::get<1>(vals), std::get<2>(vals));
                        ++temp_i2;
                    }
                );

                auto it_v_color = node_colors.find(v_name);
                // Sanity check.
                if (it_v_color == node_colors.end()) {
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Node v_name = %d was not found in the node_colors map. Skipping path.", v_name);
                    node_colors[v_name] = COLOR_TRAVERSED;
                    continue;
                } else if (it_v_color->second == COLOR_SELECTED) {
                    // The sink node was already visited, but it's part of a
                    // selected set of nodes and edges.
                    // This means all nodes on the current DFS path need to
                    // be marked as COLOR_SELECTED, because they are reachable.
                    // Also, do not extend the DFS down this line any further.
                    for (const auto& dfs_path_tuple: dfs_path) {
                        int64_t w_name = std::get<0>(dfs_path_tuple);
                        int64_t w_via_edge = std::get<1>(dfs_path_tuple);
                        int64_t w_depth = std::get<2>(dfs_path_tuple);
                        node_colors[w_name] = COLOR_SELECTED;
                        if (w_via_edge >= 0) {
                            edge_colors[w_via_edge] = COLOR_SELECTED;
                        }
                    }
                    continue;
                } else if (it_v_color->second != COLOR_BLANK) {
                    // The sink node was already visited, don't continue
                    // with DFS extension down this path again.
                    continue;
                }

                // If we're here, the node is COLOR_BLANK and hasn't been visited so far.
                // Mark it as traversed.
                node_colors[v_name] = COLOR_TRAVERSED;

                // If we're here, then the DFS needs to be extended further.
                // Add all the valid out edges.
                auto v_out_e_items = ssg->GetOutEdges(v_name);
                DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("  - Looking at next depth to extend DFS. Candidates size = %ld.\n", v_out_e_items.size()));
                for (const auto& e_item: v_out_e_items) {
                    auto it_e_color = edge_colors.find(e_item->name());

                    // Perform filtering first. If all filters passed, then
                    // mark the edge and extend DFS.
                    if (it_e_color == edge_colors.end()) {
                        // Sanity check.
                        // All edges should be in the unordered_map, but check just in case.
                        // Something weird happened, and the edge is not in the map.
                        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Edge name %d was not found in the edge_colors map. Skipping path.", e_item->name());
                        continue;
                    } else if (it_e_color->second != COLOR_BLANK) {
                        // Even though the edge color can be "selected", if we
                        // arrived at a node which wasn't colored, but the out edge
                        // is colored, then something is wrong probably.
                        // Somehow, this edge was processed already, and it could
                        // likely be some repeat cycle.
                        continue;
                    }

                    // Extend the DFS with the sink node. We don't care here
                    // if the sink node was visited, because we have a
                    // recursion-exit condition before this loop.

                    edge_colors[e_item->name()] = COLOR_TRAVERSED;
                    dfs_queue.emplace_back(std::make_tuple(e_item->sink_name(), e_item->name(), v_depth + 1));
                }
            }
            DEBUG_RUN(verbose_debug_qid, LOG_NOHEADER("\n"));
        }
    }

    // Filter out the selected nodes and edges.
    std::unordered_map<int64_t, int8_t> selected_nodes;
    for (const auto it: node_colors) {
        if (it.second == COLOR_SELECTED) {
            selected_nodes[it.first] = it.second;
        }
    }
    std::unordered_map<int64_t, int8_t> selected_edges;
    for (const auto it: edge_colors) {
        if (it.second == COLOR_SELECTED) {
            selected_edges[it.first] = it.second;
        }
    }

    // Extract the subgraph as a new graph.
    // The nodes of the aln_graph will have the same name as in the original graph.
    aln_graph = ssg->ExtractSubgraph(selected_nodes, selected_edges);

    // For every source, we need to modify the corresponding SplitSegmentNode in the
    // aln_graph, to start at the coordinate of the anchor.
    // This has to be done for every sink in the AnchorGraph too.
    // NOTE: We need to implement a check beforehand to avoid multiple sources on
    // the same SplitSegment node. This could mess things up greatly for this
    // augmentation procedure, and is meaningless as far as I can think of.
    // NOTE2: Multiple sources on the same SplitSegmentNode can be realistic. I should
    // pick the minimum coordinate one for clipping.

    LOG_NOHEADER("Selected nodes:\n");
    for (auto it: selected_nodes) {
        LOG_NOHEADER("  - Node name: %ld, color: %d\n", it.first, it.second);
    }

    LOG_NOHEADER("Selected edges:\n");
    for (auto it: selected_edges) {
        LOG_NOHEADER("  - Edge name: %ld, color: %d\n", it.first, it.second);
    }

    return aln_graph;
}

}
