/*
 * local_graph_edge.cc
 *
 *  Created on: Dec 22, 2017
 *      Author: Ivan Sovic
 */

#include <graph/anchor_graph_edge.h>

namespace raptor {

std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(SegmentEdgePtr _segment_edge) {
    return std::shared_ptr<AnchorGraphEdge>(new AnchorGraphEdge(_segment_edge));
}

std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges) {
    return std::shared_ptr<AnchorGraphEdge>(new AnchorGraphEdge(_segment_edges));
}

std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges, double _score, int32_t _path_id) {
    return std::shared_ptr<AnchorGraphEdge>(new AnchorGraphEdge(_segment_edges, false, false, _score, _path_id));
}

std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse) {
    return std::shared_ptr<AnchorGraphEdge>(new AnchorGraphEdge({}, true, is_reverse, 0.0, -1));
}

std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse, double _score, int32_t _path_id) {
    return std::shared_ptr<AnchorGraphEdge>(new AnchorGraphEdge({}, true, is_reverse, _score, _path_id));
}

std::string AnchorGraphEdge::ToDOT() const {
    return std::string();
}

std::string AnchorGraphEdge::ToJSON() const {
    std::ostringstream ss;
    ss << "{"
        << "\"implicit\":" << is_implicit()
        << ",\"implicit_rev\":" << is_implicit_reverse()
        << ",\"score\":" << score()
        << ",\"path_id\":" << path_id()
        << ",\"segment_edges\":[";
    for (const auto& seg_edge: segment_edges_) {
        ss << seg_edge->ToJSON();
    }
    ss << "]"
        << "}";
    return ss.str();
// std::vector<SegmentEdgePtr> segment_edges_;
// bool is_implicit_;
// bool is_implicit_reverse_;
// double score_;
// int32_t path_id_;
}

}
