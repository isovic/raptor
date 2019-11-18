/*
 * segment_graph.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEGMENT_GRAPH_H_
#define SRC_SEGMENT_GRAPH_H_

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>

#include <algorithm/graph.hpp>
#include <graph/segment_node.h>
#include <graph/segment_edge.h>
#include <types/typedefs_interval_tree.h>

namespace raptor {

class SegmentGraph;

typedef raptor::Graph<int64_t, raptor::SegmentNode, raptor::SegmentEdge> GraphBase;
typedef std::shared_ptr<raptor::SegmentGraph> SegmentGraphPtr;

SegmentGraphPtr createSegmentGraph();

class SegmentGraph : public GraphBase {
public:
    friend SegmentGraphPtr createSegmentGraph();
    ~SegmentGraph() = default;

    bool AddNode(const int64_t& node_name, const std::shared_ptr<raptor::SegmentNode>& node_data);

    SegmentNodePtr GetNodeByHeader(const std::string& header) const;
    SegmentNodePtr GetNodeBySeqId(int64_t seq_id) const;

    void BuildSegmentTrees();

    std::vector<SegmentEdgePtr> FindSegmentOutputEdges(int64_t segment_name, int64_t start, int64_t end);
    std::vector<SegmentEdgePtr> FindSegmentInputEdges(int64_t segment_name, int64_t start, int64_t end);

    std::string Verbose();

private:
    SegmentGraph();
    SegmentGraph(const SegmentGraph&) = delete;
    SegmentGraph& operator=(const SegmentGraph&) = delete;

    std::vector<IntervalTreeInt64> out_interval_trees_;
    std::vector<IntervalTreeInt64> in_interval_trees_;

    std::unordered_map<std::string, int64_t> header_to_internal_node_id_;
    std::unordered_map<int64_t, int64_t> seq_id_to_internal_node_id_;
};

}

#endif
