/*
 * split_segment_graph.h
 *
 *  Created on: Feb 01, 2019
 *      Author: Ivan Sovic
 *
 * Split segment graph is a broken-down representation of the SegmentGraph.
 * It doesn't allow internal edges. Instead, each out- or -in edge is
 * a location were segments are broken into pecies.
 */

#ifndef SRC_SPLIT_SEGMENT_GRAPH_H_
#define SRC_SPLIT_SEGMENT_GRAPH_H_

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>

#include <algorithm/graph.hpp>
#include <graph/segment_graph.h>
#include <types/typedefs_interval_tree.h>
#include <algorithm/node_data_base.hpp>
#include <algorithm/edge_data_base.hpp>

namespace raptor {

class SplitSegmentGraph;
class SplitSegmentNode;
class SplitSegmentEdge;

using SplitSegmentGraphNameType = int64_t;
using SplitSegmentGraphPtr = std::shared_ptr<raptor::SplitSegmentGraph>;
using SplitSegmentGraphBase = raptor::Graph<SplitSegmentGraphNameType, raptor::SplitSegmentNode, raptor::SplitSegmentEdge>; // SplitSegmentGraphBase;

SplitSegmentGraphPtr createSplitSegmentGraph();
SplitSegmentGraphPtr createSplitSegmentGraph(const SegmentGraphPtr& seg_graph);

// The node should be minimalistic. There can be many
// junctions, and we can really fragment the input segments a lot.
// std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode();
std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode(int64_t _seg_name, int32_t _t_id, int32_t _start, int32_t _end, bool _is_rev);
class SplitSegmentNode : public NodeDataBase {
  public:
    // friend std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode();
    friend std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode(int64_t _seg_name, int32_t _t_id, int32_t _start, int32_t _end, bool _is_rev);

    SplitSegmentNode()
        :   seg_name_(0), t_id_(0), start_(0), end_(0), is_rev_(false)
    { }
    SplitSegmentNode(int64_t _seg_name, int32_t _t_id, int32_t _start, int32_t _end, bool _is_rev)
        :   seg_name_(_seg_name), t_id_(_t_id), start_(_start), end_(_end), is_rev_(_is_rev)
    { }

    std::string ToDOT() const {
        return std::string("[]");
    }

    std::string ToJSON() const {
        std::ostringstream ss;
        ss << "{"
            << "\"seg_name\":" << seg_name()
            << ",\"tid\":" << t_id()
            << ",\"tr\":" << is_rev()
            << ",\"ts\":" << start()
            << ",\"te\":" << end()
            << "}";
        return ss.str();
    }

    std::string Verbose() const {
        std::ostringstream oss;
        oss << "(SplitSegmentNode) [t_id:" << t_id_ << ", s:" << start_ << ", e:" << end_ << ", is_rev:" << is_rev_ << "]";
        return oss.str();
    }

    void seg_name(int64_t val) { seg_name_ = val; }
    void t_id(int32_t val) { t_id_ = val; }
    void start(int32_t val) { start_ = val; }
    void end(int32_t val) { end_ = val; }
    void is_rev(bool val) { is_rev_ = val; }

    int64_t seg_name() const { return seg_name_; }
    int32_t t_id() const { return t_id_; }
    int32_t start() const { return start_; }
    int32_t end() const { return end_; }
    bool is_rev() const { return is_rev_; }

  private:
    int64_t seg_name_;      // Segment name from the original segment graph.
    int32_t t_id_;          // Target ID, as indexed from the reference sequences.
    int32_t start_;         // Start in the target. Zero based. In the strand of the sequence.
    int32_t end_;           // End in the target (Non-inclusive).
    bool is_rev_;           // Strand where this node originates from.
};
// inline std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode() {
//     return std::shared_ptr<raptor::SplitSegmentNode>(new SplitSegmentNode());
// }
inline std::shared_ptr<raptor::SplitSegmentNode> createSplitSegmentNode(int64_t _seg_name, int32_t _t_id, int32_t _start, int32_t _end, bool _is_rev) {
    return std::shared_ptr<raptor::SplitSegmentNode>(new raptor::SplitSegmentNode(_seg_name, _t_id, _start, _end, _is_rev));
}

std::shared_ptr<raptor::SplitSegmentEdge> createSplitSegmentEdge();
class SplitSegmentEdge : public EdgeDataBase {
  public:
    friend std::shared_ptr<raptor::SplitSegmentEdge> createSplitSegmentEdge();

    std::string ToDOT() const {
        return std::string();
    }

    std::string ToJSON() const {
        return std::string("[]");
    }

    std::string Verbose() const {
        return std::string();
    }

    SplitSegmentEdge() { }
};
inline std::shared_ptr<raptor::SplitSegmentEdge> createSplitSegmentEdge() {
    return std::shared_ptr<raptor::SplitSegmentEdge>(new raptor::SplitSegmentEdge());
}



class SplitSegmentGraph : public SplitSegmentGraphBase {
public:
    friend SplitSegmentGraphPtr createSplitSegmentGraph();
    friend SplitSegmentGraphPtr createSplitSegmentGraph(const SegmentGraphPtr& seg_graph);
    ~SplitSegmentGraph() override = default;

    bool AddNode(const SplitSegmentGraphNameType& node_name, const std::shared_ptr<SplitSegmentNode>& node_data) override;
    std::string Verbose() override;
    void Clear() override;

    static void FormatNodeKeyStrings(
                                const std::string& source_name, bool source_is_rev, int64_t source_coord,
                                const std::string& sink_name, bool sink_is_rev, int64_t sink_coord,
                                std::string& ret_key_source,
                                std::string& ret_key_sink);

    /*
     * Finds the names of all nodes within the given coordinates on a target SegmentGraphNode (before splitting).
     * t_id - Target ID.
     * t_rev - Is target reversed?
     * start - 0-based start position for lookup.
     * end - 0-based, non-inclusive end position for lookup.
    */
    std::vector<raptor::SplitSegmentGraphNameType> FindSplitNodesByTargetId(int64_t t_id, bool t_rev, int64_t start, int64_t end) const;
    // std::vector<SplitSegmentGraphNameType> FindSplitNodesBySegmentName(int64_t segment_name, int64_t start, int64_t end);

    raptor::SplitSegmentGraphPtr ExtractSubgraph(
                    const std::unordered_map<raptor::SplitSegmentGraphNameType, int8_t>& selected_nodes,
                    const std::unordered_map<raptor::SplitSegmentGraphNameType, int8_t>& selected_edges);
    /*
     * This function is used to locate the new SplitSegmentGraph node for an original SegmentGraphEdge object.
     * Example of how a key string is formatted:
     *  std::string key_source = std::to_string(edge_item->source_name()) + ((edge_data->source_is_rev()) ? "R" : "F") + std::to_string(edge_data->source_start()) + "E";
     * Returns true if the node is found, otherwise false. The node name is returned via parameter. If the return value
     * is false, then ret_node_name is undefined.
    */
    bool FindSplitNodeByKeyString(const std::string& node_key_string, raptor::SplitSegmentGraphNameType& ret_node_name) const;

    bool FindSegmentEdgeSourceAndSink(const raptor::SegmentEdgePtr& seg_edge, int64_t& ssg_source, int64_t& ssg_sink) const;

private:
    SplitSegmentGraph() { }
    SplitSegmentGraph(const SplitSegmentGraph&) = delete;
    SplitSegmentGraph& operator=(const SplitSegmentGraph&) = delete;
    void CreateFromSegmentGraph_(const SegmentGraphPtr& seg_graph);

    // Key is the "node name" (which in our case is int64_t, and is the name of the new split node.
    // This value can then be used in combination with either GetNode or GetEdge.
    std::unordered_map<SplitSegmentGraphNameType, raptor::IntervalTreeInt64> node_fwd_interval_trees_;
    std::unordered_map<SplitSegmentGraphNameType, raptor::IntervalTreeInt64> node_rev_interval_trees_;
    std::unordered_map<std::string, int32_t> node_pos_map_;
};





// using SplitSegmentGraphPtr = std::shared_ptr<raptor::SplitSegmentGraph>; // SplitSegmentGraphPtr;

// SplitSegmentGraphPtr createSplitSegmentGraph();

// class SplitSegmentGraph : public SplitSegmentGraphBase {
// public:
//     friend SplitSegmentGraphPtr createSplitSegmentGraph();
//     ~SplitSegmentGraph() = default;

//     bool AddNode(const int64_t& node_name, const std::shared_ptr<raptor::SplitSegmentNode>& node_data);

//     // SegmentNodePtr GetNodeByHeader(const std::string& header) const;
//     // SegmentNodePtr GetNodeBySeqId(int64_t seq_id) const;

//     // void BuildSegmentTrees();

//     std::string Verbose();

// private:
//     SplitSegmentGraph();
//     SplitSegmentGraph(const SplitSegmentGraph&) = delete;
//     SplitSegmentGraph& operator=(const SplitSegmentGraph&) = delete;
// };

}

#endif
