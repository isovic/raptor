/*
 * split_segment_graph.cc
 *
 *  Created on: Feb 1, 2019
 *      Author: Ivan Sovic
 */

#include <graph/split_segment_graph.h>

namespace raptor {

raptor::SplitSegmentGraphPtr createSplitSegmentGraph() {
    raptor::SplitSegmentGraphPtr ret(new raptor::SplitSegmentGraph());
    return ret;
}

raptor::SplitSegmentGraphPtr createSplitSegmentGraph(const raptor::SegmentGraphPtr& seg_graph) {
    raptor::SplitSegmentGraphPtr ret(new raptor::SplitSegmentGraph());
    ret->CreateFromSegmentGraph_(seg_graph);
    return ret;
}

std::string SplitSegmentGraph::Verbose() {
    return SplitSegmentGraphBase::Verbose();
}

void SplitSegmentGraph::Clear() {
    SplitSegmentGraphBase::Clear();
    node_fwd_interval_trees_.clear();
    node_rev_interval_trees_.clear();
    node_pos_map_.clear();
}

void SplitSegmentGraph::FormatNodeKeyStrings(
                                    const std::string& source_name, bool source_is_rev, int64_t source_coord,
                                    const std::string& sink_name, bool sink_is_rev, int64_t sink_coord,
                                    std::string& ret_key_source,
                                    std::string& ret_key_sink) {
    ret_key_source = source_name + ((source_is_rev) ? "R" : "F") + std::to_string(source_coord) + "E";
    ret_key_sink = sink_name + ((sink_is_rev) ? "R" : "F") + std::to_string(sink_coord) + "S";
}

bool SplitSegmentGraph::AddNode(const SplitSegmentGraphNameType& node_name, const std::shared_ptr<SplitSegmentNode>& node_data) {
    bool rv = SplitSegmentGraphBase::AddNode(node_name, node_data);

    if (rv == false) {
        return false;
    }

    std::string seg_node_name_as_str = std::to_string(node_data->seg_name());

    std::string key_start = seg_node_name_as_str + (node_data->is_rev() ? "R" : "F") + std::to_string(node_data->start()) + "S";
    node_pos_map_[key_start] = node_name;

    std::string key_end = seg_node_name_as_str + (node_data->is_rev() ? "R" : "F") + std::to_string(node_data->end()) + "E";
    node_pos_map_[key_end] = node_name;

    return true;
}

void SplitSegmentGraph::CreateFromSegmentGraph_(const SegmentGraphPtr& seg_graph) {

    using SplitEdgeItem = raptor::EdgeItem<int64_t, SegmentEdge>;
    using SplitEdgeItemPtr = std::shared_ptr<SplitEdgeItem>;

    Clear();

    // This string will encode the following: <seq_id>[FR]<position>[SE],
    // where [FR] is forward or reverse, and SE is either the segment start or end.
    // I'm doing this as a shortcut to actually constructing a proper data structure
    // with so many keys.
    node_pos_map_.clear();

    // Break the segments in to pieces and map them.
    for (const auto& seg_item: seg_graph->nodes()) {
        if (seg_item->is_removed()) {
            continue;
        }

        std::string seg_node_name_as_str = std::to_string(seg_item->name());

        std::vector<SplitEdgeItemPtr> in_edges = seg_graph->GetInEdges(seg_item->name());
        std::vector<SplitEdgeItemPtr> out_edges = seg_graph->GetOutEdges(seg_item->name());

        std::vector<int32_t> fwd_positions, rev_positions;
        for (const auto& edge_item: in_edges) {
            // Input edges enter the sequence at the specified position, so we
            // shouldn't offset the position.
            if (edge_item->data()->sink_is_rev() == false) {
                // The "+ 1" is so that the end is not inclusive.
                fwd_positions.emplace_back(edge_item->data()->sink_end());
            } else {
                // No need to reverse the coordinates manually because
                // the edge should have coords in the strand of the node.
                rev_positions.emplace_back(edge_item->data()->sink_end());
            }
        }
        for (const auto& edge_item: out_edges) {
            // Output edges exit the sequence AFTER the specified position.
            // To make that position inclusive, we need to add a +1.
            if (edge_item->data()->source_is_rev() == false) {
                fwd_positions.emplace_back(edge_item->data()->source_start());
            } else {
                // No need to reverse the coordinates manually because
                // the edge should have coords in the strand of the node.
                rev_positions.emplace_back(edge_item->data()->source_start());
            }
        }

        fwd_positions.emplace_back(0);
        fwd_positions.emplace_back(seg_item->data()->len());
        std::sort(fwd_positions.begin(), fwd_positions.end());

        rev_positions.emplace_back(0);
        rev_positions.emplace_back(seg_item->data()->len());
        std::sort(rev_positions.begin(), rev_positions.end());

        // Process the forward coordinates.
        IntervalVectorInt64 split_fwd_node_intervals;
        int64_t prev_fwd_node_name = -1;
        for (size_t pos_id = 1; pos_id < fwd_positions.size(); ++pos_id) {
            if (fwd_positions[pos_id-1] == fwd_positions[pos_id]) {
                continue;
            }

            // Add a new node.
            int32_t start = fwd_positions[pos_id-1];
            int32_t end = fwd_positions[pos_id];
            std::shared_ptr<raptor::SplitSegmentNode> new_node_data = raptor::createSplitSegmentNode(seg_item->name(), seg_item->data()->seq_id(), start, end, false);
            int64_t new_node_name = static_cast<int64_t>(nodes_.size());

            AddNode(new_node_name, new_node_data);

            split_fwd_node_intervals.emplace_back(IntervalInt64(start, end - 1, new_node_name));    // End in IntervalTree is inclusive.

            // Connect neighboring nodes.
            if (prev_fwd_node_name >= 0) {
                std::shared_ptr<raptor::SplitSegmentEdge> new_edge_data = raptor::createSplitSegmentEdge();
                AddEdge(prev_fwd_node_name, new_node_name, new_edge_data);
            }
            prev_fwd_node_name = new_node_name;
        }
        node_fwd_interval_trees_[seg_item->data()->seq_id()] = IntervalTreeInt64(std::move(split_fwd_node_intervals));

        // Process the reverse coordinates.
        IntervalVectorInt64 split_rev_node_intervals;
        int64_t prev_rev_node_name = -1;
        for (size_t pos_id = 1; pos_id < rev_positions.size(); ++pos_id) {
            if (rev_positions[pos_id-1] == rev_positions[pos_id]) {
                continue;
            }

            // Add a new node.
            int32_t start = rev_positions[pos_id-1];
            int32_t end = rev_positions[pos_id];
            std::shared_ptr<raptor::SplitSegmentNode> new_node_data = raptor::createSplitSegmentNode(seg_item->name(), seg_item->data()->seq_id(), start, end, true);
            int64_t new_node_name = static_cast<int64_t>(nodes_.size());

            AddNode(new_node_name, new_node_data);

            split_rev_node_intervals.emplace_back(IntervalInt64(start, end - 1, new_node_name));     // End in IntervalTree is inclusive.

            // Connect neighboring nodes.
            if (prev_rev_node_name >= 0) {
                std::shared_ptr<raptor::SplitSegmentEdge> new_edge_data = raptor::createSplitSegmentEdge();
                AddEdge(prev_rev_node_name, new_node_name, new_edge_data);
            }
            prev_rev_node_name = new_node_name;
        }
        node_rev_interval_trees_[seg_item->data()->seq_id()] = IntervalTreeInt64(std::move(split_rev_node_intervals));
    }

    for (size_t e_id = 0; e_id < seg_graph->edges().size(); ++e_id) {
        const auto& edge_item = seg_graph->edges()[e_id];
        if (edge_item->is_removed()) {
            continue;
        }

        const auto& edge_data = edge_item->data();

        std::string key_source, key_sink;
        SplitSegmentGraph::FormatNodeKeyStrings(std::to_string(edge_item->source_name()), edge_data->source_is_rev(), edge_data->source_start(),
                                    std::to_string(edge_item->sink_name()), edge_data->sink_is_rev(), edge_data->sink_start(),
                                    key_source, key_sink);

        // std::string key_source = std::to_string(edge_item->source_name()) + ((edge_data->source_is_rev()) ? "R" : "F") + std::to_string(edge_data->source_start()) + "E";
        int64_t new_source_name = -1;
        bool rv_find_source = FindSplitNodeByKeyString(key_source, new_source_name);
        if (rv_find_source == false) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge source node '%s'! This is strange, because all should have been indexed. Skipping edge: '%s'.\n", key_source.c_str(), edge_data->Verbose().c_str());
            continue;
        }

        // std::string key_sink = std::to_string(edge_item->sink_name()) + ((edge_data->sink_is_rev()) ? "R" : "F") + std::to_string(edge_data->sink_start()) + "S";
        int64_t new_sink_name = -1;
        bool rv_find_sink = FindSplitNodeByKeyString(key_sink, new_sink_name);
        if (rv_find_sink == false) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge sink node '%s'! This is strange, because all should have been indexed. Skipping edge: '%s'.\n", key_sink.c_str(), edge_data->Verbose().c_str());
            continue;
        }

        // auto it_map_source = node_pos_map_.find(key_source);
        // if (it_map_source == node_pos_map_.end()) {
        //     WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge source node '%s'! This is strange, because all should have been indexed. Skipping edge: '%s'.\n", key_source.c_str(), edge_data->Verbose().c_str());
        //     continue;
        // }

        // auto it_map_sink = node_pos_map_.find(key_sink);
        // if (it_map_sink == node_pos_map_.end()) {
        //     WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge sink node '%s'! This is strange, because all should have been indexed. Skipping edge: '%s'.\n", key_sink.c_str(), edge_data->Verbose().c_str());
        //     continue;
        // }

        // int64_t new_source_name = it_map_source->second;
        // int64_t new_sink_name = it_map_sink->second;

        std::shared_ptr<raptor::SplitSegmentEdge> new_edge_data = raptor::createSplitSegmentEdge();
        AddEdge(new_source_name, new_sink_name, new_edge_data);

    }
}

bool SplitSegmentGraph::FindSplitNodeByKeyString(const std::string& node_key_string, raptor::SplitSegmentGraphNameType& ret_node_name) const {
    /*
     * Example of how a key string is formatted:
     *  std::string key_source = std::to_string(edge_item->source_name()) + ((edge_data->source_is_rev()) ? "R" : "F") + std::to_string(edge_data->source_start()) + "E";
    */
    auto it_map_node = node_pos_map_.find(node_key_string);
    if (it_map_node == node_pos_map_.end()) {
        return false;
    }
    ret_node_name = it_map_node->second;
    return true;
}

std::vector<raptor::SplitSegmentGraphNameType> SplitSegmentGraph::FindSplitNodesByTargetId(int64_t t_id, bool t_rev, int64_t start, int64_t end) const {
    IntervalVectorInt64 intervals;

    if (t_rev == false) {
        auto it_trees = node_fwd_interval_trees_.find(t_id);
        if (it_trees == node_fwd_interval_trees_.end()) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find t_id = %ld when looking up fwd split node intervals. Returning empty.", t_id);
            return {};
        }
        // The IntervalTree API uses INCLUSIVE end coordinates for some reason.
        // We dislike that and our end coordinates
        // are always non-inclusive.
        intervals = it_trees->second.findOverlapping(start, end - 1);
    } else {
        auto it_trees = node_rev_interval_trees_.find(t_id);
        if (it_trees == node_rev_interval_trees_.end()) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find t_id = %ld when looking up rev_ split node intervals. Returning empty.", t_id);
            return {};
        }
        // The IntervalTree API uses INCLUSIVE end coordinates for some reason.
        // We dislike that and our end coordinates
        // are always non-inclusive.
        intervals = it_trees->second.findOverlapping(start, end - 1);
    }

    std::vector<raptor::SplitSegmentGraphNameType> ret;
    ret.reserve(intervals.size());
    for (const auto& interval: intervals) {
        ret.emplace_back(interval.value);
    }

    return ret;
}

bool SplitSegmentGraph::FindSegmentEdgeSourceAndSink(const raptor::SegmentEdgePtr& seg_edge, int64_t& ssg_source, int64_t& ssg_sink) const {
    ssg_source = -1;
    ssg_sink = -1;

    std::string key_source, key_sink;
    SplitSegmentGraph::FormatNodeKeyStrings(std::to_string(seg_edge->source_id()), seg_edge->source_is_rev(), seg_edge->source_start(),
                                std::to_string(seg_edge->sink_id()), seg_edge->sink_is_rev(), seg_edge->sink_start(),
                                key_source, key_sink);

    // std::string key_source = std::to_string(seg_edge->source_id()) + ((seg_edge->source_is_rev()) ? "R" : "F") + std::to_string(seg_edge->source_start()) + "E";
    bool rv_find_source = this->FindSplitNodeByKeyString(key_source, ssg_source);
    if (rv_find_source == false) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge source node '%s'!\n", seg_edge->Verbose().c_str());
        return true;
    }

    // std::string key_sink = std::to_string(seg_edge->sink_id()) + ((seg_edge->sink_is_rev()) ? "R" : "F") + std::to_string(seg_edge->sink_start()) + "S";
    bool rv_find_sink = this->FindSplitNodeByKeyString(key_sink, ssg_sink);
    if (rv_find_sink == false) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Cannot find edge sink node '%s'!\n", seg_edge->Verbose().c_str());
        return false;
    }

    return true;
}

raptor::SplitSegmentGraphPtr SplitSegmentGraph::ExtractSubgraph(
                const std::unordered_map<raptor::SplitSegmentGraphNameType, int8_t>& selected_nodes,
                const std::unordered_map<raptor::SplitSegmentGraphNameType, int8_t>& selected_edges) {

    // Create the subgraph composed only of the selected nodes and edges.
    raptor::SplitSegmentGraphPtr subgraph = raptor::createSplitSegmentGraph();

    for (const auto& node_item: this->nodes()) {
        if (node_item->is_removed()) {
            continue;
        }
        auto it = selected_nodes.find(node_item->name());
        if (it == selected_nodes.end()) {
            continue;
        }
        auto color = it->second;
        if (color == 0) {
            continue;
        }
        subgraph->AddNode(node_item->name(), node_item->data());
    }
    for (const auto& edge_item: this->edges()) {
        if (edge_item->is_removed()) {
            continue;
        }
        auto it = selected_edges.find(edge_item->name());
        if (it == selected_edges.end()) {
            continue;
        }
        auto color = it->second;
        if (color == 0) {
            continue;
        }

        // No need for anything fancy, just add it. Node names are identical as in the original graph.
        subgraph->AddEdge(edge_item->source_name(), edge_item->sink_name(), edge_item->data());
    }

    return subgraph;
}

}
