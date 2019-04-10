/*
 * segment_graph.cc
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#include <graph/segment_graph.h>

namespace raptor {

SegmentGraphPtr createSegmentGraph() {
    return SegmentGraphPtr(new raptor::SegmentGraph());
}

raptor::SegmentGraph::SegmentGraph() : GraphBase() {

}

bool SegmentGraph::AddNode(const int64_t& node_name, const std::shared_ptr<raptor::SegmentNode>& node_data) {
    int64_t num_nodes = nodes().size();

    if (GetNodeByHeader(node_data->header()) != nullptr) {
        LOG_ALL("Node with identical header already exists in the graph. Node is added, but the header will link only to the new node. Header: '%s', seq_id = %ld.\n", node_data->header().c_str(), node_data->seq_id());
    }
    if (GetNodeBySeqId(node_data->seq_id()) != nullptr) {
        LOG_ALL("Node with identical seq_id already exists in the graph. Node is added, but the seq_id will link only to the new node. Header: '%s', seq_id = %ld.\n", node_data->header().c_str(), node_data->seq_id());
    }


    bool ret = GraphBase::AddNode(node_name, node_data);

    header_to_internal_node_id_[node_data->header()] = num_nodes;
    seq_id_to_internal_node_id_[node_data->seq_id()] = num_nodes;

    return ret;
}

SegmentNodePtr SegmentGraph::GetNodeByHeader(const std::string& header) const {
    auto it_node = header_to_internal_node_id_.find(header);
    if (it_node == header_to_internal_node_id_.end()) {
        return nullptr;
    }
    return nodes_[it_node->second]->data();
}

SegmentNodePtr SegmentGraph::GetNodeBySeqId(int64_t seq_id) const {
    auto it_node = seq_id_to_internal_node_id_.find(seq_id);
    if (it_node == seq_id_to_internal_node_id_.end()) {
        return nullptr;
    }
    return nodes_[it_node->second]->data();
}

void raptor::SegmentGraph::BuildSegmentTrees() {
    out_interval_trees_.clear();
    in_interval_trees_.clear();

    out_interval_trees_.reserve(nodes_.size());
    in_interval_trees_.reserve(nodes_.size());

    for (size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
        std::vector<raptor::IntervalInt64> out_intervals;
        for (auto& edge_id: out_edge_list_[node_id]) {
            auto& edge = edges_[edge_id]->data();
            out_intervals.emplace_back(raptor::IntervalInt64(edge->source_start(), edge->source_end(), edge_id));
        }
        out_interval_trees_.emplace_back(raptor::IntervalTreeInt64(out_intervals));

        std::vector<raptor::IntervalInt64> in_intervals;
        for (auto& edge_id: in_edge_list_[node_id]) {
            auto& edge = edges_[edge_id]->data();
            in_intervals.emplace_back(raptor::IntervalInt64(edge->sink_start(), edge->sink_end(), edge_id));
        }
        in_interval_trees_.emplace_back(raptor::IntervalTreeInt64(in_intervals));
    }
}

std::vector<SegmentEdgePtr> raptor::SegmentGraph::FindSegmentOutputEdges(int64_t segment_name, int64_t start, int64_t end) {
    std::vector<SegmentEdgePtr> edges;

    // Find the internal node ID.
    uint64_t node_id = FindNodeId_(segment_name);
    if (node_id == GRAPH_UNDEFINED_NODE_ID) {
        return edges;
    }

    // Sanity check.
    if (out_interval_trees_.size() <= node_id) {
        return edges;
    }

    // Find the out edges within the window.
    std::vector<raptor::IntervalInt64> intervals;
    out_interval_trees_[node_id].findOverlapping(start, end, intervals);

    for (auto& interval: intervals) {
        uint64_t edge_id = interval.value;
        auto& edge = edges_[edge_id]->data();
        edges.emplace_back(edge);
    }

    return edges;
}

std::vector<SegmentEdgePtr> raptor::SegmentGraph::FindSegmentInputEdges(int64_t segment_name, int64_t start, int64_t end) {
    std::vector<SegmentEdgePtr> edges;

    // Find the internal node ID.
    uint64_t node_id = FindNodeId_(segment_name);
    if (node_id == GRAPH_UNDEFINED_NODE_ID) {
        return edges;
    }

    // Sanity check.
    if (in_interval_trees_.size() <= node_id) {
        return edges;
    }

    // Find the out edges within the window.
    std::vector<raptor::IntervalInt64> intervals;
    in_interval_trees_[node_id].findOverlapping(start, end, intervals);

    for (auto& interval: intervals) {
        uint64_t edge_id = interval.value;
        auto& edge = edges_[edge_id]->data();
        edges.emplace_back(edge);
    }

    return edges;
}

std::string raptor::SegmentGraph::Verbose() {
    std::ostringstream oss;
    oss << GraphBase::Verbose();
    oss << std::endl;
    oss << "  All nodes:" << std::endl;
    for (size_t i = 0; i < nodes_.size(); i++) {
        auto& node_item = nodes_[i];
        auto node_id = node_item->name();
        auto& node = node_item->data();
        oss << "    [i = " << i << "] node_name = " << node_id << ", " << node->Verbose() << std::endl;
    }
    oss << std::endl;
    oss << "  All edges:" << std::endl;
    for (size_t i = 0; i < edges_.size(); i++) {
        auto& edge_item = edges_[i];
        auto source_node_id = edge_item->source_name();
        auto sink_node_id = edge_item->sink_name();
        auto& edge = edge_item->data();
        oss << "    [i = " << i << "] source_node_id = " << source_node_id << ", sink_node_id = " << sink_node_id << ", " << edge->Verbose() << std::endl;
    }

    return oss.str();
}

}
