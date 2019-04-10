/*
 * local_graph.h
 *
 *  Created on: Dec 29, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_LOCAL_GRAPH_H_
#define SRC_RAPTOR_LOCAL_GRAPH_H_

#include <graph/anchor_graph_edge.h>
#include <graph/anchor_graph_node.hpp>
#include <algorithm/graph.hpp>
#include <algorithm/path.hpp>
#include <algorithm/node_item.hpp>
#include <algorithm/edge_item.hpp>
#include <containers/region/region_mapped.h>
#include <cstdint>

namespace raptor {

using AnchorGraph = raptor::Graph<int64_t, raptor::AnchorGraphNode, raptor::AnchorGraphEdge>;
using NodeItemAnchorGraph = raptor::NodeItem<int64_t, raptor::AnchorGraphNode>;
using EdgeItemAnchorGraph = raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>;

using AnchorGraphPtr = std::shared_ptr<raptor::AnchorGraph>;
using NodeItemAnchorGraphPtr = std::shared_ptr<raptor::NodeItemAnchorGraph>;
using EdgeItemAnchorGraphPtr = std::shared_ptr<raptor::EdgeItemAnchorGraph>;

inline std::shared_ptr<AnchorGraph> createAnchorGraph() {
    return AnchorGraph::createGraph();
}

}

#endif
