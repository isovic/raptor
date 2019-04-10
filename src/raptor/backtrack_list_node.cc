/*
 * backtrack_list_node.cc
 *
 *  Created on: Jan 08, 2018
 *      Author: Ivan Sovic
 */

#include <raptor/backtrack_list_node.h>

namespace raptor {

std::shared_ptr<BacktrackListNode> createBacktrackListNode(
    int64_t _dp_id, int64_t _curr_node_id, int64_t _pred_node_id,
    std::shared_ptr<raptor::EdgeItem<int64_t, raptor::AnchorGraphEdge>> _edge) {
    return std::shared_ptr<BacktrackListNode>(
        new BacktrackListNode(_dp_id, _curr_node_id, _pred_node_id, _edge));
}

}  // namespace raptor
