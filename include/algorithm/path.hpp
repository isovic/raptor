/*
 * path.hpp
 *
 *  Created on: Dec 26, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_PATH_H_
#define SRC_PATH_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <algorithm/graph.hpp>
#include <algorithm/node_item.hpp>
#include <algorithm/edge_item.hpp>

namespace raptor {

template<class NameType, class NodeDataType, class EdgeDataType>
class Path;

template<class NameType, class NodeDataType, class EdgeDataType>
class Path {
public:
    typedef std::shared_ptr<NodeItem<NameType, NodeDataType>> NodeItemPtr;
    typedef std::shared_ptr<EdgeItem<NameType, EdgeDataType>> EdgeItemPtr;

    /*
     * Factory function.
    */
    static std::shared_ptr<raptor::Path<NameType, NodeDataType, EdgeDataType>> createPath(const NameType& node_name, const std::shared_ptr<NodeDataType> node_data) {
        return std::shared_ptr<raptor::Path<NameType, NodeDataType, EdgeDataType>>(
                                    new raptor::Path<NameType, NodeDataType, EdgeDataType>(node_name, node_data));
    }

    virtual ~Path() = default;

    /*
     * Adds a node to the path, linked with it's predecessor with edge_data.
    */
    bool Add(const NameType& node_name, const std::shared_ptr<NodeDataType> node_data, const std::shared_ptr<EdgeDataType> edge_data) {
        int64_t new_node_id = (int64_t) nodes_.size();
        name_to_node_id_[node_name] = nodes_.size();

        auto new_node = NodeItem<NameType, NodeDataType>::createNodeItem(new_node_id, node_name, node_data);
        nodes_.emplace_back(new_node);

        auto new_edge = EdgeItem<NameType, EdgeDataType>::createEdgeItem(edges_.size(), static_cast<int64_t>(edges_.size()), last_added_node_, node_name, new_node_id - 1, new_node_id, edge_data);
        edges_.emplace_back(new_edge);

        last_added_node_ = node_name;

        return true;
    }

    const std::vector<NodeItemPtr>& nodes() const {
        return nodes_;
    }

    const std::vector<EdgeItemPtr>& edges() const {
        return edges_;
    }

    int64_t score() const {
        return score_;
    }

    void nodes(const std::vector<NodeItemPtr>& _nodes) {
        nodes_ = _nodes;
    }

    void edges(const std::vector<EdgeItemPtr>& _edges) {
        edges_ = _edges;
    }

    void score(int64_t _score) {
        score_ = _score;
    }

    std::string Verbose() const {
        std::ostringstream oss;

        oss << "Path info:" << std::endl;
        oss << "  Score: " << score_ << std::endl;
        oss << "  Nodes: " << nodes_.size() << std::endl;
        oss << "  Edges: " << edges_.size() << std::endl;
        oss << std::endl;
        if (nodes_.size() == 0) {
            oss << "  Path is empty." << std::endl;
            return oss.str();
        }

        oss << "  Exact path:" << std::endl;
        oss << "    [0] Node: name = " << nodes_[0]->name() << std::endl;
        oss << "            Node data: " << nodes_[0]->data()->Verbose() << std::endl;

        for (size_t i = 1; i < nodes_.size(); i++) {
            auto& edge_item = edges_[i - 1];
            auto edge_source = edge_item->source_name();
            auto edge_sink = edge_item->sink_name();
            auto& edge = edge_item->data();
            oss << "    '" << nodes_[i-1]->name() << "' -> '" << nodes_[i]->name() << "'" << std::endl;
            oss <<"            Edge data: " << edge->Verbose() << std::endl;

            auto& node_item = nodes_[i];
            auto node_id = node_item->name();
            auto& node = node_item->data();
            oss << "    [" << i << "] Node: name = " << node_id << std::endl;
            oss << "            Node data: " << node->Verbose() << std::endl;
        }

        oss << std::endl;

        return oss.str();
    }

private:
    Path(const Path&) = delete;
    Path& operator=(const Path&) = delete;

    Path(const NameType& node_name, const std::shared_ptr<NodeDataType> node_data)
                                                    : last_added_node_(node_name),
                                                        nodes_{},
                                                        edges_{},
                                                        name_to_node_id_{},
                                                        score_(-1) {

        int64_t new_node_id = (int64_t) nodes_.size();
        name_to_node_id_[node_name] = nodes_.size();

        auto new_node = NodeItem<NameType, NodeDataType>::createNodeItem(new_node_id, node_name, node_data);
        nodes_.emplace_back(new_node);
    }

    // Keep track of the last added node, so that we can add an edge.
    NameType last_added_node_;
    // Vector of all nodes.
    std::vector<NodeItemPtr> nodes_;
    // Vector of all edges.
    std::vector<EdgeItemPtr> edges_;
    // Relation map from node name to the internal ID of the node.
    std::unordered_map<NameType, uint64_t> name_to_node_id_;

    int64_t score_;

};

}

#endif
