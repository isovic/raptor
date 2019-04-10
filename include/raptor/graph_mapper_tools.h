/*
 * graph_mapper_tools.h
 *
 *  Created on: Jan 26, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_GRAPH_MAPPER_TOOLS_H_
#define SRC_RAPTOR_GRAPH_MAPPER_TOOLS_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <graph/anchor_graph.h>
#include <graph/anchor_graph_edge.h>
#include <types/typedefs.h>
#include <containers/mapping_result/linear_mapping_result.h>

namespace raptor {
namespace graphmapper {

struct GraphEdgePenalties {
    uint32_t status_flag = 0;

    // Sink data.
    int64_t w_qstart = 0;
    int64_t w_qend = 0;
    int64_t w_tstart = 0;
    int64_t w_tend = 0;
    int64_t w_l = 0;

    int64_t v_qstart = 0;
    int64_t v_qend = 0;
    int64_t v_tstart = 0;
    int64_t v_tend = 0;
    int64_t v_l = 0;
    int64_t v_virtual_tstart = 0;
    int64_t v_virtual_tend = 0;
    int64_t v_virtual_l = 0;
    int64_t v_pivot_tpos = 0;
    int64_t v_score_compensation = 0;

    double score_vw = 0.0;
    double edge_score = 0.0;

    int64_t w_qspan() {
        return w_qend - w_qstart;
    }
    int64_t w_tspan() {
        return w_tend - w_tstart;
    }
    int64_t v_qspan() {
        return v_qend - v_qstart;
    }
    int64_t v_tspan() {
        return v_tend - v_tstart;
    }

    int64_t dist_q() {
        return w_qstart - v_qend;
    }
    // int64_t dist_t() {
    //     return w_tstart - v_tend;
    // }
    int64_t dist_virtual_t() {
        return w_tstart - v_virtual_tend;
    }

    int64_t gap_dist() {
        return std::abs(dist_q() - dist_virtual_t());
        // return (dist_q() < dist_t()) ? (dist_t() - dist_q()) : (dist_q() - dist_t());
    }

};

/*
 * Calculates the penalty of traversing an edge from node_i to node_j.
*/
GraphEdgePenalties CalcAnchorGraphEdgePenalty(const raptor::AnchorPtr& node_i, const raptor::AnchorPtr& node_j,
                    const std::shared_ptr<raptor::AnchorGraphEdge>& local_graph_edge,
                    int32_t node_v_score, int32_t node_w_score, bool verbose_debug);

/*
 * BreakAnchors finds any anchor that has an edge junction in the graph within itself, and breaks
 * it into pieces. The breaking is not performed heuristically, instead the seed hits are grouped
 * into smaller anchors (before and after the edge join), and a new vector of anchors is returned.
 *
 * BreakAnchors looks up graph fork positions within the anchor coordinate span,
 * and breaks the anchor by selecting seed hits before fork as one anchor, and
 * hits after the fork as the other anchor.
 * This is required to make the GraphDP more accurate.
*/
std::shared_ptr<raptor::LinearMappingResult> BreakAnchors(
                    const raptor::GraphPtr graph,
                    const std::shared_ptr<raptor::LinearMappingResult>& input_linear_result,
                    // const std::vector<raptor::TargetAnchorPtr>& target_anchors,
                    // const std::vector<raptor::ChainPtr>& target_hits,
                    int32_t k);

std::vector<std::vector<int64_t>> FindBestLeafNodes(
                    const std::shared_ptr<raptor::AnchorGraph>& graph,
                    double allowed_score_diff_frac,
                    bool verbose_debug_qid);

}
}

#endif