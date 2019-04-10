/*
 * backtrack_list_node.h
 *
 *  Created on: Jan 08, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_BACKTRACK_LIST_NODE_H_
#define SRC_GRAPH_BACKTRACK_LIST_NODE_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <algorithm/edge_item.hpp>
#include <graph/anchor_graph.h>
#include <graph/local_path.h>

namespace raptor {

/*
 * The BacktrackListNode is a single-linked-list node used
 * to make the backtracking path tracking more efficient.
 */
class BacktrackListNode;

std::shared_ptr<BacktrackListNode> createBacktrackListNode(
    int64_t _dp_id, int64_t _curr_node_id, int64_t _pred_node_id,
    std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> _edge);

class BacktrackListNode {
   public:
    friend std::shared_ptr<BacktrackListNode> createBacktrackListNode(
        int64_t _dp_id, int64_t _curr_node_id, int64_t _pred_node_id,
        std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> _edge);

    ~BacktrackListNode() = default;

    int64_t dp_id() const { return dp_id_; }
    int64_t curr_node_id() const { return curr_node_id_; }
    int64_t pred_node_id() const { return pred_node_id_; }
    const std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> edge() const {
        return edge_;
    }

    void dp_id(int64_t _dp_id) { dp_id_ = _dp_id; }
    void curr_node_id(int64_t _curr_node_id) { curr_node_id_ = _curr_node_id; }
    void pred_node_id(int64_t _pred_node_id) { pred_node_id_ = _pred_node_id; }

   private:
    BacktrackListNode(int64_t _dp_id, int64_t _curr_node_id, int64_t _pred_node_id,
                      std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> _edge)
        : dp_id_(_dp_id),
          curr_node_id_(_curr_node_id),
          pred_node_id_(_pred_node_id),
          edge_(_edge) {}

    BacktrackListNode(const BacktrackListNode&) = delete;
    BacktrackListNode& operator=(const BacktrackListNode&) = delete;

    int64_t dp_id_;
    int64_t curr_node_id_;
    int64_t pred_node_id_;
    std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> edge_;
};

}  // namespace raptor

#endif
