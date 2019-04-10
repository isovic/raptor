#ifndef TESTS_UNIT_GRAPH_MAPPER_H_
#define TESTS_UNIT_GRAPH_MAPPER_H_
#include <raptor/graph_mapper.h>

#include <cstdint>
#include <vector>

namespace raptor {
namespace unit {
/*
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "100000", "0"},
    };
*/
void AddNodesToGraph(raptor::GraphPtr& graph, const std::vector<std::vector<std::string>>& nodes_sts);

/*
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "edge1", "ref1+", "ref1+", "100000$", "100000$", "0", "0", "*"},    // Circular edge
    };
*/
void AddEdgesToGraph(raptor::GraphPtr& graph, const std::vector<std::vector<std::string>>& edges_gfa2);

/*
    std::vector<std::vector<int32_t>> anchors_m4_plus = {
        // {qid, tid, score, idt,   qrev, qstart, qend, qlen,   trev, tstart, tend, tlen   cov_bases_q, cov_bases_t, num_seeds}
        {0, 0, -1, -1,    0,    0, 1000, 2000,     0, 99000, 100000, 100000,   1000, 1000, 10, 0},
        {0, 0, -1, -1,    0, 1000, 2000, 2000,     0,     0,   1000, 100000,   1000, 1000, 10, 1},
    };
*/
raptor::TargetAnchorPtrVector ConstructTargetAnchors(std::vector<std::vector<int32_t>> anchors_m4_plus); // Copies for sorting.

raptor::SegmentGraphPtr WrapConstructSegmentGraph(const std::string& test_name,
                                                const std::vector<std::vector<std::string>>& nodes_str,
                                                const std::vector<std::vector<std::string>>& edges_gfa2,
                                                bool verbose_debug_qid);

raptor::AnchorGraphPtr WrapConstructAnchorGraph(const std::string& test_name,
                                                const raptor::SegmentGraphPtr& graph,
                                                const std::vector<std::vector<int32_t>>& anchors_m4_plus,
                                                int32_t max_allowed_dist,
                                                int32_t predecessor_lookup_window,
                                                int32_t max_path_edges,
                                                int32_t chain_max_skip,
                                                bool verbose_debug_qid);

}
}

#endif
