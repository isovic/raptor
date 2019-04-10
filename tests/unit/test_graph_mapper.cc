#include <gtest/gtest.h>
#include <raptor/graph_mapper.h>
#include <graph/segment_graph_parser.h>

#include <cstdint>
#include <algorithm>
#include <set>
#include <vector>
#include <log/log_tools.h>

#include <tests/unit/test_graph_mapper.h>

// mindex::IndexPtr CreateSampleIndex() {
//     int32_t k = 5;
//     int32_t w = 1;
//     bool hp_supp = false;
//     int32_t max_hp_len = 10;
//     double freq_percentile = 0.002;
//     int32_t min_occ_cutoff = -1;
//     bool index_only_fwd_strand = false;

//     auto index_params = mindex::createIndexParams(k, w, hp_supp,max_hp_len, freq_percentile, min_occ_cutoff, index_only_fwd_strand);

//     auto index = mindex::createMinimizerIndex(index_params);

//     // Define the reference sequences and their headers.
//     std::vector<std::string> seqs = {
//                                     "AAAAAAAAAA",
//                                     };
//     std::vector<std::string> headers = {
//                                     "fake_seq1",
//                                     };

//     // Add them to the index. This only adds the data, but does not build the index.
//     index->AddSequences(seqs, headers);

//     // Build the index from the added sequences.
//     index->Build();

//     return index;
// }

namespace raptor {
namespace unit {

// struct AnchorM4Plus {
//     int32_t qid = 0, tid = 0,
// };

        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
raptor::AnchorPtr ConvertM4PlusToRegionMapped(
        int32_t qid, int32_t tid, double score, double idt,
        bool qrev, int32_t qstart, int32_t qend, int32_t qlen,
        bool trev, int32_t tstart, int32_t tend, int32_t tlen,
        int32_t cov_bases_q, int32_t cov_bases_t, int32_t num_seeds,
        int32_t target_hits_id,
        int32_t anchor_id
        ) {

    auto map_env = raptor::createMappingEnv(tid, 0, tlen, trev, qid, qlen, qrev);
    auto new_anchor = raptor::createRegionMapped(anchor_id, target_hits_id, map_env,
                                                    qstart, qend, tstart, tend,
                                                    cov_bases_q, cov_bases_t, num_seeds,
                                                    -1, -1, -1, -1, -1);                    // Path IDs. Not used here.
    return new_anchor;

}

raptor::AnchorPtr ConvertAnchorM4PlusToRegionMapped(const std::vector<int32_t>& anchors_m4_plus_line) {
    const auto& anchor = anchors_m4_plus_line;

    int32_t qid = anchor[0];
    bool qrev = anchor[4];
    int32_t qstart = anchor[5];
    int32_t qend = anchor[6];
    int32_t qlen = anchor[7];

    int32_t tid = anchor[1];
    bool trev = anchor[8];
    int32_t tstart = anchor[9];
    int32_t tend = anchor[10];
    int32_t tlen = anchor[11];

    int32_t cov_bases_q = anchor[12];
    int32_t cov_bases_t = anchor[13];
    int32_t num_seeds = anchor[14];

    int32_t anchor_id = anchor[15];
    int32_t target_hits_id = -1;    // We aren't generating a list of k-mer hits here.

    auto map_env = raptor::createMappingEnv(tid, 0, tlen, trev, qid, qlen, qrev);

    auto new_anchor = raptor::createRegionMapped(anchor_id, target_hits_id, map_env,
                                                    qstart, qend, tstart, tend,
                                                    cov_bases_q, cov_bases_t, num_seeds,
                                                    -1, -1, -1, -1, -1);                    // Path IDs. Not used here.

    return new_anchor;
}

void AddNodesToGraph(raptor::GraphPtr& graph, const std::vector<std::vector<std::string>>& nodes_str) {
    for (const auto& node_str: nodes_str) {
        // Add nodes.
        int64_t seq_len = std::stol(node_str[1]);
        int64_t seq_id = std::stol(node_str[2]);
        graph->AddNode(seq_id, raptor::createSegmentNode(node_str[0], seq_len, false, seq_id));
    }
}

void AddEdgesToGraph(raptor::GraphPtr& graph, const std::vector<std::vector<std::string>>& edges_gfa2) {
    for (const auto& edge_gfa2: edges_gfa2) {
        raptor::SegmentEdgePtr edge = raptor::GraphLoader::ParseGFA2Edge(graph, edge_gfa2);
        EXPECT_NE(edge, nullptr);
        graph->AddEdge(edge->source_id(), edge->sink_id(), edge);
    }
}

raptor::TargetAnchorPtrVector ConstructTargetAnchors(std::vector<std::vector<int32_t>> anchors_m4_plus) { // Copies for sorting.

    raptor::TargetAnchorPtrVector target_anchors;

    // Sorting is required to construct the local graph properly.
    std::sort(anchors_m4_plus.begin(), anchors_m4_plus.end(),
                    [](const std::vector<int32_t>& a, const std::vector<int32_t>& b) {
                    return (a[1] < b[1] ||  // Different target.
                                (a[1] == b[1] &&    // Same target.
                                    (a[8] < b[8] || // Different strand.
                                        (a[8] == b[8] && a[9] < b[9]))  // Same strand, sort by start coord.
                                )
                            );
                    });

    std::shared_ptr<raptor::MappingEnv> map_env = nullptr;
    std::shared_ptr<raptor::TargetAnchorType> target_anchor_container = nullptr;

    // This line is only used to verify the input. Not actually testing the accuracy of output.
    if (anchors_m4_plus.size() == static_cast<size_t>(16)) {
        std::cerr << "ERROR: Test is wrong, not all fields are specified! Exiting.\n";
        exit(1);
    }

    for (size_t i = 0; i < anchors_m4_plus.size(); ++i) {

        auto new_anchor = ConvertAnchorM4PlusToRegionMapped(anchors_m4_plus[i]);

        if (i == 0 || (i > 0 && anchors_m4_plus[i][1] != anchors_m4_plus[i-1][1])) {
            map_env = new_anchor->env();
            target_anchor_container = std::shared_ptr<raptor::TargetAnchorType>(
                                            new raptor::TargetAnchorType(map_env));
            target_anchors.emplace_back(target_anchor_container);
        }

        new_anchor->env(map_env);

        target_anchors.back()->hits().emplace_back(new_anchor);

        // const auto& anchor = anchors_m4_plus[i];

        // int32_t qid = anchor[0];
        // bool qrev = anchor[4];
        // int32_t qstart = anchor[5];
        // int32_t qend = anchor[6];
        // int32_t qlen = anchor[7];

        // int32_t tid = anchor[1];
        // bool trev = anchor[8];
        // int32_t tstart = anchor[9];
        // int32_t tend = anchor[10];
        // int32_t tlen = anchor[11];

        // int32_t cov_bases_q = anchor[12];
        // int32_t cov_bases_t = anchor[13];
        // int32_t num_seeds = anchor[14];

        // int32_t anchor_id = anchor[15];
        // int32_t target_hits_id = -1;    // We aren't generating a list of k-mer hits here.

        // if (i == 0 || (i > 0 && anchor[1] != anchors_m4_plus[i-1][1])) {
        //     map_env = raptor::createMappingEnv(tid, 0, tlen, trev, qid, qlen, qrev);

        //     target_anchor_container = std::shared_ptr<raptor::TargetAnchorType>(
        //                                     new raptor::TargetAnchorType(map_env));
        //     target_anchors.emplace_back(target_anchor_container);
        // }

        // auto new_anchor = raptor::createRegionMapped(anchor_id, target_hits_id, map_env,
        //                                                 qstart, qend, tstart, tend,
        //                                                 cov_bases_q, cov_bases_t, num_seeds,
        //                                                 -1, -1, -1, -1, -1);                    // Path IDs. Not used here.
        // target_anchors.back()->hits().emplace_back(new_anchor);
    }

    return target_anchors;
}

void CompareEdges(const raptor::AnchorGraphPtr& local_graph, const std::vector<std::tuple<int32_t, int32_t, std::string>>& expected_edges) {
    std::set<std::tuple<int32_t, int32_t, std::string>> edge_set;
    for (const auto& edge_item: local_graph->edges()) {
        if (edge_item->is_removed()) {
            continue;
        }
        const auto& edge = edge_item->data();
        const auto& source = local_graph->GetNode(edge_item->source_name())->anchor();
        const auto& sink = local_graph->GetNode(edge_item->sink_name())->anchor();
        edge_set.emplace(std::make_tuple(source->id(), sink->id(), edge->SummarizeSegmentEdgesAsString()));

        std::cerr << "source->id() = " << source->id() << ", sink->id() = " << sink->id() << "\n";
        std::cerr << edge->SummarizeSegmentEdgesAsString() << "\n";
    }
    EXPECT_EQ(edge_set.size(), expected_edges.size());
    for (const auto exp_edge: expected_edges) {
        EXPECT_NE(edge_set.find(exp_edge), edge_set.end());
    }
}

raptor::SegmentGraphPtr WrapConstructSegmentGraph(const std::string& test_name,
                                                const std::vector<std::vector<std::string>>& nodes_str,
                                                const std::vector<std::vector<std::string>>& edges_gfa2,
                                                bool verbose_debug_qid
                                                ) {
    // Create an empty graph.
    raptor::SegmentGraphPtr graph = raptor::createSegmentGraph();

    // Add nodes.
    AddNodesToGraph(graph, nodes_str);

    // Add edges.
    AddEdgesToGraph(graph, edges_gfa2);

    // Construct the lookup interval trees. This should in the future be done automatically when
    // new node or edge are added.
    graph->BuildSegmentTrees();

    if (verbose_debug_qid) {
        LOG_ALL("SegmentGraph:\n%s\n", graph->Verbose().c_str());
    }

    return graph;
}

raptor::AnchorGraphPtr WrapConstructAnchorGraph(const std::string& test_name,
                                                const raptor::SegmentGraphPtr& graph,
                                                const std::vector<std::vector<int32_t>>& anchors_m4_plus,
                                                int32_t max_allowed_dist,
                                                int32_t predecessor_lookup_window,
                                                int32_t max_path_edges,
                                                int32_t chain_max_skip,
                                                bool verbose_debug_qid
                                                ) {

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Create the anchors from the M4-like formatted vector above.
    auto target_anchors = ConstructTargetAnchors(anchors_m4_plus);

    std::shared_ptr<raptor::AnchorGraph> anchor_graph =
                raptor::GraphMapper::ConstructAnchorGraph(
                                            graph,
                                            target_anchors,
                                            max_allowed_dist,
                                            predecessor_lookup_window,
                                            max_path_edges,
                                            chain_max_skip,
                                            verbose_debug_qid);

    if (verbose_debug_qid) {
        LOG_ALL("Local graph:\n%s\n", anchor_graph->Verbose().c_str());
    }
	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    // std::ofstream ofs_graph("test.graph." + test_name + ".html");
    // local_graph->WriteAsHTML(ofs_graph);

    return anchor_graph;
}



TEST(GraphMapper, ConstructAnchorGraph1) {
    std::string test_name("ConstructAnchorGraph1");
    /*
    Simple example for a circular genome.
    */
    int32_t max_allowed_dist = 10000;
    int32_t predecessor_lookup_window = 10000;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "edge1", "ref1+", "ref1+", "100000$", "100000$", "0", "0", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        {0, 0, -1, -1,    0,    0, 1000, 2000,     0, 99000, 100000, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0, 1000, 2000, 2000,     0,     0,   1000, 100000,   1000, 1000, 10, 1},
    };

    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        {0, 1, "['edge1']"},
    };

    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);


}

TEST(GraphMapper, ConstructAnchorGraph2) {
    std::string test_name("ConstructAnchorGraph2");
    /*
    Testing simple detection of implicit edges, without explicit segmental edges.
    */
    int32_t max_allowed_dist = 10000;
    int32_t predecessor_lookup_window = 10000;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}

        // Test 1. These should be linked.
        {0, 0, -1, -1,    0,    0, 1000, 2000,     0, 1000, 2000, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0, 1000, 2000, 2000,     0, 2100, 3000, 100000,   1000, 1000, 10, 1},

        // Test 2. Separate test, mapping onto different target. Both qid and tid are different. This should still link the anchors.
        {1, 1, -1, -1,    0,    0, 1000, 4000,     0, 1000, 2000, 100000,   1000, 1000, 10, 2},
        {1, 1, -1, -1,    0, 3000, 4000, 4000,     0, 3000, 4000, 100000,   1000, 1000, 10, 3},

        // Test 3. Separate test, mapping onto different target. Both qid and tid are different. This should still link the anchors.
        {2, 2, -1, -1,    0,    0, 1000, 15000,     0, 1000, 2000, 100000,   1000, 1000, 10, 4},
        {2, 2, -1, -1,    0, 9000, 10000, 15000,     0, 9000, 10000, 100000,   1000, 1000, 10, 5},

        // Test 4. These are too far apart (max_allowed_dist).
        {3, 3, -1, -1,    0,    0, 1000, 15000,     0, 1000, 2000, 100000,   1000, 1000, 10, 6},
        {3, 3, -1, -1,    0, 11001, 12000, 15000,     0, 12001, 13000, 100000,   1000, 1000, 10, 7},

        // Test 5. Another example which should produce no edges because the anchors are too far apart.
        {4, 4, -1, -1,    0,    0, 1000, 2000,     0, 1000, 2000, 100000,   1000, 1000, 10, 8},
        {4, 4, -1, -1,    0, 1000, 2000, 2000,     0, 20100, 30000, 100000,   1000, 1000, 10, 9},
    };

    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        // Test 1.
        {0, 1, "[]"},
        // Test 2.
        {2, 3, "[]"},
        // Test 3.
        {4, 5, "[]"},
        // Test 4 has no edges.
        // Test 5 has no edges.
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}

TEST(GraphMapper, ConstructAnchorGraph3) {
    std::string test_name("ConstructAnchorGraph3");
    /*
    Simple transcriptome graph.
    */
    int32_t max_allowed_dist = 10000;
    int32_t predecessor_lookup_window = 10000;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
        {"ref2", "100000", "1"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "edge1.1", "ref1+", "ref1+", "11000", "11000", "20000", "20000", "*"},
        {"E", "edge1.2", "ref1+", "ref1+", "20400", "20400", "80000", "80000", "*"},
        {"E", "edge1.3", "ref1+", "ref1+", "80100", "80100", "89990", "89990", "*"},

        {"E", "edge2.1", "ref2+", "ref2+", "11000", "11000", "20000", "20000", "*"},
        {"E", "edge2.2", "ref2+", "ref2+", "20400", "20400", "80000", "80000", "*"},
        {"E", "edge2.3", "ref2+", "ref2+", "80100", "80100", "89990", "89990", "*"},
        {"E", "edge2.4", "ref2+", "ref2+", "90505", "90505", "10000", "10000", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // Test 1. Transcripts that are far apart.
        {0, 0, -1, -1,    0,    0, 1000, 3000,     0, 10000, 11000, 100000,   1000, 1000, 10, 0}, // 0
        {0, 0, -1, -1,    0, 2000, 2400, 3000,     0, 20000, 20400, 100000,   1000, 1000, 10, 1}, // 1
        {0, 0, -1, -1,    0, 2400, 2500, 3000,     0, 80000, 80100, 100000,   1000, 1000, 10, 2}, // 2
        {0, 0, -1, -1,    0, 2573, 3000, 3000,     0, 90000, 90500, 100000,   1000, 1000, 10, 3}, // 3

        // Test 2, circRNA.
        {1, 1, -1, -1,    0, 1000, 1400, 3000,     0, 20000, 20400, 100000,   1000, 1000, 10, 4}, // 4
        {1, 1, -1, -1,    0, 1400, 1500, 3000,     0, 80000, 80100, 100000,   1000, 1000, 10, 5}, // 5
        {1, 1, -1, -1,    0, 1573, 2000, 3000,     0, 90000, 90500, 100000,   1000, 1000, 10, 6}, // 6
        {1, 1, -1, -1,    0, 2001, 3000, 3000,     0, 10001, 11000, 100000,   1000, 1000, 10, 7}, // 7
    };

    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        {0, 1, "['edge1.1']"},
        {0, 2, "['edge1.1','edge1.2']"},
        {0, 3, "['edge1.1','edge1.2','edge1.3']"},
        {1, 2, "['edge1.2']"},
        {1, 3, "['edge1.2','edge1.3']"},
        {2, 3, "['edge1.3']"},

        {4, 5, "['edge2.2']"},
        {4, 6, "['edge2.2','edge2.3']"},
        {4, 7, "['edge2.2','edge2.3','edge2.4']"},

        {5, 6, "['edge2.3']"},
        {5, 7, "['edge2.3','edge2.4']"},

        {6, 7, "['edge2.4']"},

        // {7, 5, "['edge2.4','edge2.1']"},
        // {7, 6, "['edge2.1','edge2.2','edge2.3']"},

        // // Explicit edges.
        // {0, 1},
        // {0, 2},             // This one is weird, but 11000 + 10000 == 21000 and the next edge is at 20400, so we actually pick edge2
        // {1, 2}, {1, 3},
        // {2, 3}, // The {1, 3} is just by a margin.
        // // Multi-hop edges.
        // {0, 2}, {0, 3},
        // {0, 3},             // This one is the same weirdness as the one before, we pick up anchor1 -> edge2 -> anchor2 -> edge3 -> anchor3
        // // Implicit edges.
        // {0, 1}, {2, 3}
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}

TEST(GraphMapper, ConstructAnchorGraph4) {
    std::string test_name("ConstructAnchorGraph4");
    /*
    Graph with a simple cycle.
    */
    int32_t max_allowed_dist = 1000;
    int32_t predecessor_lookup_window = 1000;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
        {"ref2", "100000", "1"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "edge1", "ref1+", "ref1+", "45000", "45000", "40000", "40000", "*"},    // Circular edge
        {"E", "edge2", "ref2+", "ref2+", "45000", "45000", "40000", "40000", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // Test 1, a simple cycle.
        {0, 0, -1, -1,    0,    0, 5000, 25000,     0, 35000, 40000, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0,    5000, 10000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10, 1},
        {0, 0, -1, -1,    0,    10000, 15000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10, 2},
        {0, 0, -1, -1,    0,    15000, 20000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10, 3},
        {0, 0, -1, -1,    0,    20000, 25000, 25000,     0, 45000, 50000, 100000,   1000, 1000, 10, 4},

        // Test 2, a simple cycle, but one anchor is missing.
        {1, 1, -1, -1,    0,    0, 5000, 25000,     0, 35000, 40000, 100000,   1000, 1000, 10, 5},
        {1, 1, -1, -1,    0,    5000, 10000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10, 6},
        // {1, 1, -1, -1,    0,    10000, 15000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10},
        {1, 1, -1, -1,    0,    15000, 20000, 25000,     0, 40000, 45000, 100000,   1000, 1000, 10, 7},
        {1, 1, -1, -1,    0,    20000, 25000, 25000,     0, 45000, 50000, 100000,   1000, 1000, 10, 8},
    };

    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        // Test 1, the entire chain is resolved.
        {0, 1, "[]"}, {1, 2, "['edge1']"}, {2, 3, "['edge1']"}, {3, 4, "[]"},
        // Test 2, there is no middle anchor, so there are two chains.
        {5, 6, "[]"}, {7, 8, "[]"},
    };

    /*
           ┌-=====-┐
           │       │
    o=====-o-=====-o-=====o
           │       │
           └-=====-┘
    */
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}

TEST(GraphMapper, ConstructAnchorGraph5) {
    std::string test_name("ConstructAnchorGraph5");
    /*
    A very small cycle.
    */
    int32_t max_allowed_dist = 100;
    int32_t predecessor_lookup_window = 100;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        //  This is a very simple small cycle.
        {"E", "edge1", "ref1+", "ref1+", "35002", "35002", "35000", "35000", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // Test 1, a simple very small cycle.
        {0, 0, -1, -1,    0,    0, 5000, 10000,     0, 30000, 34999, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0,    5050, 10000, 10000,     0, 35005, 40000, 100000,   1000, 1000, 10, 1},
    };

    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        {0, 1, "['edge1']"},
        // // This one is the implicit edge.
        // {0, 1},
        // // This one is the segment edge, going through the back edge in the graph.
        // // The developed local graph is acyclic because we can't link any predecessor anchors
        // // within the cycle (there are none). And even if we could, it wouldn't
        // // be cyclic because of query coordinates.
        // {0, 1},
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}


TEST(GraphMapper, ConstructAnchorGraph6) {
    std::string test_name("ConstructAnchorGraph6");
    /*
    Slightly more complicated graph example with cycles.
    */
    int32_t max_allowed_dist = 1000;
    int32_t predecessor_lookup_window = 1000;
    int32_t max_path_edges = -1;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        // Backwards edge.
        {"E", "E1", "ref1+", "ref1+", "35020", "35020", "35005", "35005", "*"},    // Circular edge
        // Forward edge, with a strobe.
        {"E", "E2", "ref1+", "ref1+", "35002", "35002", "35018", "35018", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // Two anchors before the cyclic construct.
        {0, 0, -1, -1,    0,    0, 5000, 10000,     0, 30000, 34999, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0,    5050, 10000, 10000,     0, 35025, 40000, 100000,   1000, 1000, 10, 1},
    };
    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        // For explanation, check notebook on 26.01.2019.
        // Result - unrolled graph, which is acyclic.
        {0, 1, "['E1']"}
        // {0, 1}, // Implicit.
        // {0, 1}, // A -> E1 -> B
        // {0, 1}, // A -> E2 -> B
        // {0, 1}, // A -> E2 -> E1 -> B
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}

TEST(GraphMapper, ConstructAnchorGraph7) {
    std::string test_name("ConstructAnchorGraph7");
    /*
    Slightly more complicated graph example with cycles.
    */
    int32_t max_allowed_dist = 1000;
    int32_t predecessor_lookup_window = 1000;
    int32_t max_path_edges = 10;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        // Backwards edge.
        {"E", "E1", "ref1+", "ref1+", "35020", "35020", "35005", "35005", "*"},    // Circular edge
        // Forward edge, with a strobe.
        {"E", "E2", "ref1+", "ref1+", "35008", "35008", "35018", "35018", "*"},    // Circular edge
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // Two anchors before the cyclic construct.
        {0, 0, -1, -1,    0,    0, 5000, 10000,     0, 30000, 34999, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0,    5050, 10000, 10000,     0, 35025, 40000, 100000,   1000, 1000, 10, 1},
    };
    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        // For explanation, check notebook on 26.01.2019.
        // Result - unrolled graph, which is acyclic.
        {0, 1, "['E1','E2','E1','E2','E1']"}
        // {0, 1}, // Implicit.
        // {0, 1}, // S1 -> E1 -> S2
        // {0, 1}, // S1 -> E2 -> S2
        // {0, 1}, // S1 -> E2 -> E1 -> S2
        // {0, 1}, // S1 -> E1 -> E2 -> S2.    This is currently buggy, we should add this one, but it get's filtered out because E1 is already visited.
        // // {0, 1}, // S1 -> E1 -> E2 -> E1 > S2.    This is currently under consideration. The ConstructAnchorGraph intentionally doesn't allow cycles without nodes in them.
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}

TEST(GraphMapper, ConstructAnchorGraph8) {
    std::string test_name("ConstructAnchorGraph8");
    /*
    Special edge cases for debugging.
    */
    int32_t max_allowed_dist = 1000;
    int32_t predecessor_lookup_window = 1000;
    int32_t max_path_edges = 10;
    int32_t chain_max_skip = 25;
    bool verbose_debug_qid = true;

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "695", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        // Deletion edge.
        {"E", "E1", "ref1+", "ref1+", "463", "463", "469", "469", "*"},
    };
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        // There is a deletion edge jumping over the middle anchor, but the graph shouldn't traverse it.
        // Instead, the query doesn't have the deletion, and the main path should be linear with the reference.
        {0, 0, -1, -1,    0,    302, 451, 690,     0, 302, 451, 695,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0,    472, 473, 690,     0, 467, 468, 695,   1000, 1000, 10, 1},
        {0, 0, -1, -1,    0,    485, 586, 690,     0, 480, 581, 695,   1000, 1000, 10, 2},
    };
    std::vector<std::tuple<int32_t, int32_t, std::string>> expected_edges = {
        {0, 2, "[]"},
        {0, 1, "[]"},
        {1, 2, "[]"}
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                verbose_debug_qid);

    raptor::AnchorGraphPtr local_graph = WrapConstructAnchorGraph(
            test_name,
            seg_graph, anchors_m4_plus,
            max_allowed_dist, predecessor_lookup_window,
            max_path_edges, chain_max_skip, verbose_debug_qid);

    CompareEdges(local_graph, expected_edges);
}



TEST(GraphMapper, CreateAnchorGraphConnectedComponents1) {
    /*
     * Empty input -> empty output.
    */
    std::string test_name("CreateAnchorGraphConnectedComponents1");

    double graph_allowed_score_diff_frac = 0.95;
    bool verbose_debug_qid = false;

    //////////////////////////////////////////////////
    /// Inputs.                                    ///
    //////////////////////////////////////////////////
    // Define anchors for the graph.
    std::vector<raptor::AnchorPtr> anchors = {
        // {{qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds, anchor_id}, {dp_score, dp_predecessor, dp_path_id}}
    };

    // Create the input graph.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = raptor::createAnchorGraph();

    // The graph is identical, but due to backtracking, the order of nodes and edges is different.
    // The graphs are isomorphic though, but we are lacking a method to compare isomorphism.
    std::vector<std::string> expected;

    //////////////////////////////////////////////////
    /// Run the UUT.                               ///
    //////////////////////////////////////////////////
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs =
                    raptor::GraphMapper::CreateAnchorGraphConnectedComponents(
                                            mapped_graph,
                                            graph_allowed_score_diff_frac,
                                            verbose_debug_qid);

    //////////////////////////////////////////////////
    /// Eval.                                      ///
    //////////////////////////////////////////////////

    std::vector<std::string> result;
    for (const auto& pgraph_tuple: pgraphs) {
        const auto& pgraph = std::get<0>(pgraph_tuple);
        result.emplace_back(pgraph->ToJSON());
    }

    EXPECT_EQ(result, expected);
}

TEST(GraphMapper, CreateAnchorGraphConnectedComponents2) {
    /*
     * A simple case of a linear mapping graph of 3 anchors, with only
     * implicit edges.
     * The output should be the same as input, isomorphically.
     * The only difference is that due to backtracking, the order of nodes and
     * edges in the output is different.
    */
    std::string test_name("CreateAnchorGraphConnectedComponents2");

    double graph_allowed_score_diff_frac = 0.95;
    bool verbose_debug_qid = false;

    //////////////////////////////////////////////////
    /// Inputs.                                    ///
    //////////////////////////////////////////////////
    // Define anchors for the graph.
    std::vector<raptor::AnchorPtr> anchors = {
        // {{qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds, anchor_id}, {dp_score, dp_predecessor, dp_path_id}}
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    100, 200, 1000,     0, 100, 200, 10000,   100, 100, 10, -1, 0),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    200, 300, 1000,     0, 200, 300, 10000,   100, 100, 10, -1, 1),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    300, 400, 1000,     0, 300, 400, 10000,   100, 100, 10, -1, 2),
    };

    // Create the input graph.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = raptor::createAnchorGraph();

    // Create input nodes.
    mapped_graph->AddNode(0, raptor::createAnchorGraphNode(anchors[0], 100.0, -1, 0));
    mapped_graph->AddNode(1, raptor::createAnchorGraphNode(anchors[1], 200.0, 0, 0));
    mapped_graph->AddNode(2, raptor::createAnchorGraphNode(anchors[2], 300.0, 1, 0));

    // Create input edges.
    mapped_graph->AddEdge(0, 1, raptor::createImplicitAnchorGraphEdge(false, 200.0, 0));
    mapped_graph->AddEdge(1, 2, raptor::createImplicitAnchorGraphEdge(false, 300.0, 0));

    //////////////////////////////////////////////////
    /// Expected output.                           ///
    //////////////////////////////////////////////////
    // The graph is identical, but due to backtracking, the order of nodes and edges is different.
    // The graphs are isomorphic though, but we are lacking a method to compare isomorphism.
    std::vector<std::string> expected;
    expected.emplace_back(
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 0, \"data\": {\"score\":300,\"predecessor\":1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,300,400,1000,0,300,400,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"score\":200,\"predecessor\":0,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,200,300,1000,0,200,300,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 2, \"data\": {\"score\":100,\"predecessor\":-1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,100,200,1000,0,100,200,10000,100,100,10\"}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 1, \"w\": 2, \"v_int_id\": 1, \"w_int_id\": 0, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":300,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 0, \"w\": 1, \"v_int_id\": 2, \"w_int_id\": 1, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":200,\"path_id\":0,\"segment_edges\":[]}}]\n"
        "],\n"
        "\"adj\": {\n"
        "2: [],\n"
        "1: [2],\n"
        "0: [1]\n"
        "}\n"
        "}\n"
        "]\n"
        "pid=0\n"
        "pid_score=300\n"
    );

    //////////////////////////////////////////////////
    /// Run the UUT.                               ///
    //////////////////////////////////////////////////
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs =
                    raptor::GraphMapper::CreateAnchorGraphConnectedComponents(
                                            mapped_graph,
                                            graph_allowed_score_diff_frac,
                                            verbose_debug_qid);

    //////////////////////////////////////////////////
    /// Eval.                                      ///
    //////////////////////////////////////////////////

    std::vector<std::string> result;
    for (const auto& pgraph_tuple: pgraphs) {
        const auto& pgraph = std::get<0>(pgraph_tuple);
        std::ostringstream ss;
        ss << pgraph->ToJSON()
            << "pid=" << std::get<1>(pgraph_tuple) << "\n"
            << "pid_score=" << std::get<2>(pgraph_tuple) << "\n";
        result.emplace_back(ss.str());
        // std::cerr << pgraph->ToJSON() << "\n";
    }

    EXPECT_EQ(result, expected);
}

TEST(GraphMapper, CreateAnchorGraphConnectedComponents3) {
    /*
     * Same as the previous test (linear path), but 2 separate pid components.
    */
    std::string test_name("CreateAnchorGraphConnectedComponents3");

    double graph_allowed_score_diff_frac = 0.95;
    bool verbose_debug_qid = false;

    //////////////////////////////////////////////////
    /// Inputs.                                    ///
    //////////////////////////////////////////////////
    // Define anchors for the graph.
    std::vector<raptor::AnchorPtr> anchors = {
        // {{qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds, anchor_id}, {dp_score, dp_predecessor, dp_path_id}}
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    100, 200, 1000,     0, 100, 200, 10000,   100, 100, 10, -1, 0),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    200, 300, 1000,     0, 200, 300, 10000,   100, 100, 10, -1, 1),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    300, 400, 1000,     0, 300, 400, 10000,   100, 100, 10, -1, 2),

        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    100, 200, 1000,     0, 100, 200, 10000,   100, 100, 10, -1, 3),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    200, 300, 1000,     0, 200, 300, 10000,   100, 100, 10, -1, 4),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    300, 400, 1000,     0, 300, 400, 10000,   100, 100, 10, -1, 5),
    };

    // Create the input graph.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = raptor::createAnchorGraph();

    // Create input nodes.
    // (node_name, (node_data, score, predecessor, path_id))
    mapped_graph->AddNode(0, raptor::createAnchorGraphNode(anchors[0], 100.0, -1, 0));
    mapped_graph->AddNode(1, raptor::createAnchorGraphNode(anchors[1], 200.0, 0, 0));
    mapped_graph->AddNode(2, raptor::createAnchorGraphNode(anchors[2], 300.0, 1, 0));
    mapped_graph->AddNode(3, raptor::createAnchorGraphNode(anchors[3], 100.0, -1, 1));
    mapped_graph->AddNode(4, raptor::createAnchorGraphNode(anchors[4], 200.0, 3, 1));
    mapped_graph->AddNode(5, raptor::createAnchorGraphNode(anchors[5], 300.0, 4, 1));

    // Create input edges.
    // (source_name, sink_name, (is_reverse, score, path_id))
    mapped_graph->AddEdge(0, 1, raptor::createImplicitAnchorGraphEdge(false, 200.0, 0));
    mapped_graph->AddEdge(1, 2, raptor::createImplicitAnchorGraphEdge(false, 300.0, 0));
    mapped_graph->AddEdge(3, 4, raptor::createImplicitAnchorGraphEdge(false, 200.0, 1));
    mapped_graph->AddEdge(4, 5, raptor::createImplicitAnchorGraphEdge(false, 300.0, 1));

    //////////////////////////////////////////////////
    /// Expected output.                           ///
    //////////////////////////////////////////////////
    std::vector<std::string> expected;
    expected.emplace_back(
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 0, \"data\": {\"score\":300,\"predecessor\":1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,300,400,1000,0,300,400,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"score\":200,\"predecessor\":0,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,200,300,1000,0,200,300,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 2, \"data\": {\"score\":100,\"predecessor\":-1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,100,200,1000,0,100,200,10000,100,100,10\"}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 1, \"w\": 2, \"v_int_id\": 1, \"w_int_id\": 0, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":300,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 0, \"w\": 1, \"v_int_id\": 2, \"w_int_id\": 1, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":200,\"path_id\":0,\"segment_edges\":[]}}]\n"
        "],\n"
        "\"adj\": {\n"
        "2: [],\n"
        "1: [2],\n"
        "0: [1]\n"
        "}\n"
        "}\n"
        "]\n"
        "pid=0\n"
        "pid_score=300\n"
    );
    expected.emplace_back(
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 0, \"data\": {\"score\":300,\"predecessor\":4,\"path_id\":1,\"anchor\":\"0,0,-1,10,0,300,400,1000,0,300,400,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 1, \"data\": {\"score\":200,\"predecessor\":3,\"path_id\":1,\"anchor\":\"0,0,-1,10,0,200,300,1000,0,200,300,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 2, \"data\": {\"score\":100,\"predecessor\":-1,\"path_id\":1,\"anchor\":\"0,0,-1,10,0,100,200,1000,0,100,200,10000,100,100,10\"}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 4, \"w\": 5, \"v_int_id\": 1, \"w_int_id\": 0, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":300,\"path_id\":1,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 3, \"w\": 4, \"v_int_id\": 2, \"w_int_id\": 1, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":200,\"path_id\":1,\"segment_edges\":[]}}]\n"
        "],\n"
        "\"adj\": {\n"
        "5: [],\n"
        "4: [5],\n"
        "3: [4]\n"
        "}\n"
        "}\n"
        "]\n"
        "pid=1\n"
        "pid_score=300\n"
    );

    //////////////////////////////////////////////////
    /// Run the UUT.                               ///
    //////////////////////////////////////////////////
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs =
                    raptor::GraphMapper::CreateAnchorGraphConnectedComponents(
                                            mapped_graph,
                                            graph_allowed_score_diff_frac,
                                            verbose_debug_qid);

    //////////////////////////////////////////////////
    /// Eval.                                      ///
    //////////////////////////////////////////////////

    std::vector<std::string> result;
    for (const auto& pgraph_tuple: pgraphs) {
        const auto& pgraph = std::get<0>(pgraph_tuple);
        std::ostringstream ss;
        ss << pgraph->ToJSON()
            << "pid=" << std::get<1>(pgraph_tuple) << "\n"
            << "pid_score=" << std::get<2>(pgraph_tuple) << "\n";
        result.emplace_back(ss.str());
        // std::cerr << pgraph->ToJSON() << "\n";
    }

    EXPECT_EQ(result, expected);
}

TEST(GraphMapper, CreateAnchorGraphConnectedComponents4) {
    /*
     * Simple linear graph with only implicit edges, but some nodes have multiple
     * input implicit edges with similar score. Only one should be chosen, the best one.
    */
    std::string test_name("CreateAnchorGraphConnectedComponents4");

    double graph_allowed_score_diff_frac = 0.95;
    bool verbose_debug_qid = false;

    //////////////////////////////////////////////////
    /// Inputs.                                    ///
    //////////////////////////////////////////////////
    // Define anchors for the graph.
    std::vector<raptor::AnchorPtr> anchors = {
        // {{qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds, anchor_id}, {dp_score, dp_predecessor, dp_path_id}}
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    100, 200, 1000,     0, 100, 200, 10000,   100, 100, 10, -1, 0),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    200, 300, 1000,     0, 200, 300, 10000,   100, 100, 10, -1, 1),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    300, 400, 1000,     0, 300, 400, 10000,   100, 100, 10, -1, 2),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    400, 500, 1000,     0, 400, 500, 10000,   100, 100, 10, -1, 3),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    500, 600, 1000,     0, 500, 600, 10000,   100, 100, 10, -1, 4),
    };

    // Create the input graph.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = raptor::createAnchorGraph();

    // Create input nodes.
    mapped_graph->AddNode(0, raptor::createAnchorGraphNode(anchors[0], 100.0, -1, 0));
    mapped_graph->AddNode(1, raptor::createAnchorGraphNode(anchors[1], 200.0, 0, 0));
    mapped_graph->AddNode(2, raptor::createAnchorGraphNode(anchors[2], 300.0, 1, 0));
    mapped_graph->AddNode(3, raptor::createAnchorGraphNode(anchors[3], 400.0, 2, 0));
    mapped_graph->AddNode(4, raptor::createAnchorGraphNode(anchors[4], 500.0, 3, 0));

    // Create input edges.
    mapped_graph->AddEdge(0, 1, raptor::createImplicitAnchorGraphEdge(false, 200.0, 0));
    mapped_graph->AddEdge(1, 2, raptor::createImplicitAnchorGraphEdge(false, 300.0, 0));
    mapped_graph->AddEdge(2, 3, raptor::createImplicitAnchorGraphEdge(false, 400.0, 0));
    mapped_graph->AddEdge(3, 4, raptor::createImplicitAnchorGraphEdge(false, 500.0, 0));

    mapped_graph->AddEdge(0, 3, raptor::createImplicitAnchorGraphEdge(false, 399.7, 0));
    mapped_graph->AddEdge(1, 3, raptor::createImplicitAnchorGraphEdge(false, 399.8, 0));
    mapped_graph->AddEdge(2, 3, raptor::createImplicitAnchorGraphEdge(false, 399.9, 0));

    //////////////////////////////////////////////////
    /// Expected output.                           ///
    //////////////////////////////////////////////////
    std::vector<std::string> expected;
    expected.emplace_back(
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 0, \"data\": {\"score\":500,\"predecessor\":3,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,500,600,1000,0,500,600,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 1, \"data\": {\"score\":400,\"predecessor\":2,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,400,500,1000,0,400,500,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"data\": {\"score\":300,\"predecessor\":1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,300,400,1000,0,300,400,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 3, \"data\": {\"score\":200,\"predecessor\":0,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,200,300,1000,0,200,300,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 4, \"data\": {\"score\":100,\"predecessor\":-1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,100,200,1000,0,100,200,10000,100,100,10\"}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 3, \"w\": 4, \"v_int_id\": 1, \"w_int_id\": 0, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":500,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 2, \"w\": 3, \"v_int_id\": 2, \"w_int_id\": 1, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":400,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"v\": 1, \"w\": 2, \"v_int_id\": 3, \"w_int_id\": 2, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":300,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"v\": 0, \"w\": 1, \"v_int_id\": 4, \"w_int_id\": 3, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":200,\"path_id\":0,\"segment_edges\":[]}}]\n"
        "],\n"
        "\"adj\": {\n"
        "4: [],\n"
        "3: [4],\n"
        "2: [3],\n"
        "1: [2],\n"
        "0: [1]\n"
        "}\n"
        "}\n"
        "]\n"
        "pid=0\n"
        "pid_score=500\n"
    );

    //////////////////////////////////////////////////
    /// Run the UUT.                               ///
    //////////////////////////////////////////////////
	// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs =
                    raptor::GraphMapper::CreateAnchorGraphConnectedComponents(
                                            mapped_graph,
                                            graph_allowed_score_diff_frac,
                                            verbose_debug_qid);
	// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    //////////////////////////////////////////////////
    /// Eval.                                      ///
    //////////////////////////////////////////////////

    std::vector<std::string> result;
    for (const auto& pgraph_tuple: pgraphs) {
        const auto& pgraph = std::get<0>(pgraph_tuple);
        std::ostringstream ss;
        ss << pgraph->ToJSON()
            << "pid=" << std::get<1>(pgraph_tuple) << "\n"
            << "pid_score=" << std::get<2>(pgraph_tuple) << "\n";
        result.emplace_back(ss.str());
        // std::cerr << pgraph->ToJSON() << "\n";
    }

    EXPECT_EQ(result, expected);

    // std::cerr << "graph:\n" << mapped_graph->ToJSON() << "\n";
}

TEST(GraphMapper, CreateAnchorGraphConnectedComponents5) {
    /*
     * Test a more complicated example containing a combination of
     * implicit and explicit edges, bubbles, spurs and multiple leaf nodes.
    */
    std::string test_name("CreateAnchorGraphConnectedComponents5");

    double graph_allowed_score_diff_frac = 0.95;
    bool verbose_debug_qid = false;

    //////////////////////////////////////////////////
    /// Inputs.                                    ///
    //////////////////////////////////////////////////
    // Define anchors for the graph.
    // The exact coordinates aren't really important, because the graph is constructed here
    // manually, and not using the methods from GraphMapper.
    std::vector<raptor::AnchorPtr> anchors = {
        // {{qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds, anchor_id}, {dp_score, dp_predecessor, dp_path_id}}
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    100, 200, 1000,     0, 100, 200, 10000,   100, 100, 10, -1, 0),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    300, 500, 1000,     0, 300, 500, 10000,   100, 100, 10, -1, 1),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    700, 800, 1000,     0, 700, 800, 10000,   100, 100, 10, -1, 2),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    1000, 1200, 1000,     0, 1000, 1200, 10000,   100, 100, 10, -1, 3),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    1500, 1600, 1000,     0, 1500, 1600, 10000,   100, 100, 10, -1, 4),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    1800, 1900, 1000,     0, 1800, 1900, 10000,   100, 100, 10, -1, 5),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    2200, 2300, 1000,     0, 2300, 2300, 10000,   100, 100, 10, -1, 6),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    2500, 2600, 1000,     0, 2600, 2600, 10000,   100, 100, 10, -1, 7),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    2800, 2700, 1000,     0, 2800, 2700, 10000,   100, 100, 10, -1, 8),
        ConvertM4PlusToRegionMapped(0, 0, -1, -1,    0,    3300, 3400, 1000,     0, 3300, 3400, 10000,   100, 100, 10, -1, 9),
    };

    // Create the input graph.
    std::shared_ptr<raptor::AnchorGraph> mapped_graph = raptor::createAnchorGraph();

    // Create input nodes.
    mapped_graph->AddNode(0, raptor::createAnchorGraphNode(anchors[0], 100.0, -1, 0));
    mapped_graph->AddNode(1, raptor::createAnchorGraphNode(anchors[1], 200.0, 0, 0));
    mapped_graph->AddNode(2, raptor::createAnchorGraphNode(anchors[2], 300.0, 1, 0));
    mapped_graph->AddNode(3, raptor::createAnchorGraphNode(anchors[3], 400.0, 2, 0));
    mapped_graph->AddNode(4, raptor::createAnchorGraphNode(anchors[4], 500.0, 3, 0));
    mapped_graph->AddNode(5, raptor::createAnchorGraphNode(anchors[5], 600.0, 4, 0));
    mapped_graph->AddNode(6, raptor::createAnchorGraphNode(anchors[6], 234.0, 5, 0));
    mapped_graph->AddNode(7, raptor::createAnchorGraphNode(anchors[7], 800.0, 6, 0));
    mapped_graph->AddNode(8, raptor::createAnchorGraphNode(anchors[8], 100.0, 7, 0));
    mapped_graph->AddNode(9, raptor::createAnchorGraphNode(anchors[9], 1000.0, 8, 0));

    // Create input edges.
    mapped_graph->AddEdge(0, 1, raptor::createImplicitAnchorGraphEdge(false, 200.0, 0));    // Implicit edge, should remain.
    mapped_graph->AddEdge(1, 2, raptor::createImplicitAnchorGraphEdge(false, 300.0, 0));    // Implicit edge, should remain.
    mapped_graph->AddEdge(0, 2, raptor::createImplicitAnchorGraphEdge(false, 299.9, 0));    // Extra implicit edge which should be removed.

    mapped_graph->AddEdge(2, 3, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 400.0, 0));    // Bubble 1, branch 1. Should remain, similar score.
    mapped_graph->AddEdge(3, 4, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 500.0, 0));    // Bubble 1, branch 1. Should remain, similar score.
    mapped_graph->AddEdge(2, 4, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 499.9, 0));    // Bubble 1, branch 2. Should remain, similar score.
    mapped_graph->AddEdge(1, 4, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 123.0, 0));    // Bubble 1, branch 2. Low score to input node 4 (score 500.0). Should be filtered out.

    mapped_graph->AddEdge(4, 5, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 600.0, 0));    // Clean edge, should remain.
    mapped_graph->AddEdge(5, 6, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 234.0, 0));    // Out spur, should be filtered out. Node 6 is a leaf node with lower score than the best leaf will have.
    mapped_graph->AddEdge(5, 7, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 800, 0));    // Clean edge, should remain.
    mapped_graph->AddEdge(8, 7, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 125.0, 0));    // In-spur. Should be filtered out, it's not the highest scoring node.
    mapped_graph->AddEdge(7, 9, raptor::createAnchorGraphEdge({raptor::createSegmentEdge()}, 1000.0, 0));    // Clean edge, should remain.

    //////////////////////////////////////////////////
    /// Expected output.                           ///
    //////////////////////////////////////////////////
    std::vector<std::string> expected;
    expected.emplace_back(
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 9, \"removed\": 0, \"int_id\": 0, \"data\": {\"score\":1000,\"predecessor\":8,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,3300,3400,1000,0,3300,3400,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 7, \"removed\": 0, \"int_id\": 1, \"data\": {\"score\":800,\"predecessor\":6,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,2500,2600,1000,0,2600,2600,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 2, \"data\": {\"score\":600,\"predecessor\":4,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,1800,1900,1000,0,1800,1900,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 3, \"data\": {\"score\":500,\"predecessor\":3,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,1500,1600,1000,0,1500,1600,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 4, \"data\": {\"score\":300,\"predecessor\":1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,700,800,1000,0,700,800,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 5, \"data\": {\"score\":400,\"predecessor\":2,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,1000,1200,1000,0,1000,1200,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 6, \"data\": {\"score\":200,\"predecessor\":0,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,300,500,1000,0,300,500,10000,100,100,10\"}}],\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 7, \"data\": {\"score\":100,\"predecessor\":-1,\"path_id\":0,\"anchor\":\"0,0,-1,10,0,100,200,1000,0,100,200,10000,100,100,10\"}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 7, \"w\": 9, \"v_int_id\": 1, \"w_int_id\": 0, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":1000,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 5, \"w\": 7, \"v_int_id\": 2, \"w_int_id\": 1, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":800,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"v\": 4, \"w\": 5, \"v_int_id\": 3, \"w_int_id\": 2, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":600,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"v\": 2, \"w\": 4, \"v_int_id\": 4, \"w_int_id\": 3, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":499.9,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"v\": 3, \"w\": 4, \"v_int_id\": 5, \"w_int_id\": 3, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":500,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"v\": 2, \"w\": 3, \"v_int_id\": 4, \"w_int_id\": 5, \"data\": {\"implicit\":0,\"implicit_rev\":0,\"score\":400,\"path_id\":0,\"segment_edges\":[{\"symname\":\"\",\"id\":-1,\"v\":,\"v_int_id\":0,\"vr\":0,\"vs\":0,\"ve\":0,\"vl\":0,\"w\":,\"w_int_id\":0,\"wr\":0,\"ws\":0,\"we\":0,\"wl\":0}]}}],\n"
        "[\"eitem\", {\"name\": 6, \"removed\": 0, \"int_id\": 6, \"v\": 1, \"w\": 2, \"v_int_id\": 6, \"w_int_id\": 4, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":300,\"path_id\":0,\"segment_edges\":[]}}],\n"
        "[\"eitem\", {\"name\": 7, \"removed\": 0, \"int_id\": 7, \"v\": 0, \"w\": 1, \"v_int_id\": 7, \"w_int_id\": 6, \"data\": {\"implicit\":1,\"implicit_rev\":0,\"score\":200,\"path_id\":0,\"segment_edges\":[]}}]\n"
        "],\n"
        "\"adj\": {\n"
        "9: [],\n"
        "7: [9],\n"
        "5: [7],\n"
        "4: [5],\n"
        "2: [4,3],\n"
        "3: [4],\n"
        "1: [2],\n"
        "0: [1]\n"
        "}\n"
        "}\n"
        "]\n"
        "pid=0\n"
        "pid_score=1000\n"
    );

    //////////////////////////////////////////////////
    /// Run the UUT.                               ///
    //////////////////////////////////////////////////
	// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);
    std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> pgraphs =
                    raptor::GraphMapper::CreateAnchorGraphConnectedComponents(
                                            mapped_graph,
                                            graph_allowed_score_diff_frac,
                                            verbose_debug_qid);
	// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    //////////////////////////////////////////////////
    /// Eval.                                      ///
    //////////////////////////////////////////////////

    std::vector<std::string> result;
    for (const auto& pgraph_tuple: pgraphs) {
        const auto& pgraph = std::get<0>(pgraph_tuple);
        std::ostringstream ss;
        ss << pgraph->ToJSON()
            << "pid=" << std::get<1>(pgraph_tuple) << "\n"
            << "pid_score=" << std::get<2>(pgraph_tuple) << "\n";
        result.emplace_back(ss.str());
        // std::cerr << pgraph->ToJSON() << "\n";
    }

    EXPECT_EQ(result, expected);

    // std::cerr << "graph:\n" << mapped_graph->ToJSON() << "\n";
}

}
}
