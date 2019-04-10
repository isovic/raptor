/*
 * graph.hpp
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_H_
#define SRC_GRAPH_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <log/log_tools.h>
#include <algorithm/node_item.hpp>
#include <algorithm/edge_item.hpp>

namespace raptor {

static constexpr __int128 GRAPH_EDGE_KEY_64bit_MASK = (((__int128) 0x0FFFFFFFFFFFFFFFF));
static constexpr uint64_t GRAPH_UNDEFINED_EDGE_ID = ((uint64_t) 0xFFFFFFFFFFFFFFFF);
static constexpr uint64_t GRAPH_UNDEFINED_NODE_ID = ((uint64_t) 0xFFFFFFFFFFFFFFFF);

template<class NameType, class NodeDataType, class EdgeDataType>
class Graph;

template<class NameType, class NodeDataType, class EdgeDataType>
class Graph {
public:
    using EdgeListType = std::vector<int64_t>;
    using NodeItemPtr = std::shared_ptr<NodeItem<NameType, NodeDataType>>;
    using EdgeItemPtr = std::shared_ptr<EdgeItem<NameType, EdgeDataType>>;

    static std::shared_ptr<raptor::Graph<NameType, NodeDataType, EdgeDataType>> createGraph() {
        return std::shared_ptr<raptor::Graph<NameType, NodeDataType, EdgeDataType>>(new raptor::Graph<NameType, NodeDataType, EdgeDataType>());
    }

    virtual ~Graph() = default;

    virtual bool AddNode(const NameType& node_name, const std::shared_ptr<NodeDataType>& node_data) {
        // Do not allow duplicate nodes.
        auto it = name_to_node_id_.find(node_name);
        if (it != name_to_node_id_.end()) {
            LOG_ALL("Warning: Node with the same name already exists, not adding again.\n");
            return false;
        }

        // Add a new node.
        name_to_node_id_[node_name] = nodes_.size();
        auto new_node = NodeItem<NameType, NodeDataType>::createNodeItem(nodes_.size(), node_name, node_data);
        nodes_.emplace_back(new_node);
        in_edge_list_.emplace_back(EdgeListType());
        out_edge_list_.emplace_back(EdgeListType());

        return true;
    }

    virtual bool AddEdge(const NameType& source_node_name, const NameType& sink_node_name, const std::shared_ptr<EdgeDataType>& edge_data) {
        uint64_t id_source = FindNodeId_(source_node_name);
        if (id_source == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Edge not added because the source node is undefined!\n");
            return false;
        }

        uint64_t id_sink = FindNodeId_(sink_node_name);
        if (id_sink == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Edge not added because the sinknode is undefined!\n");
            return false;
        }

        // Add the in and out edge to the nodes.
        out_edge_list_[id_source].emplace_back(edges_.size());
        in_edge_list_[id_sink].emplace_back(edges_.size());

        // // Annotate the edge for quick lookup.
        // auto edge_key = MakeEdgeKey_(id_source, id_sink);
        // auto edges_it = FindEdgesIterator_(id_source, id_sink);
        // if (edges_it == edge_to_edge_ids_.end()) {
        //     edge_to_edge_ids_[edge_key] = std::vector<uint64_t>{};
        // }
        // edge_to_edge_ids_[edge_key].emplace_back(edges_.size());

        // Add the edge data.
        auto new_edge = EdgeItem<NameType, EdgeDataType>::createEdgeItem(edges_.size(), static_cast<int64_t>(edges_.size()), source_node_name, sink_node_name, id_source, id_sink, edge_data);
        edges_.emplace_back(new_edge);

        return true;
    }

    inline std::shared_ptr<NodeDataType> GetNode(const NameType& node_name) const {
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the requested node.\n");
            return nullptr;
        }
        if (nodes_[node_id]->is_removed()) {
            return nullptr;
        }
        return nodes_[node_id]->data();
    }

    inline std::vector<EdgeItemPtr> GetInEdges(const NameType& node_name) const {
        std::vector<EdgeItemPtr> edges;
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested in edges.\n");
            return edges;
        }
        for (auto& edge_id: in_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false) {
                edges.emplace_back(edges_[edge_id]);
            }
        }
        return edges;
    }

    inline int64_t CountInEdges(const NameType& node_name) const {
        int64_t num_edges = 0;
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested in edges.\n");
            return num_edges;
        }
        for (auto& edge_id: in_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == true) {
                continue;
            }
            ++num_edges;
        }
        return num_edges;
    }

    inline std::vector<EdgeItemPtr> GetOutEdges(const NameType& node_name) const {
        std::vector<EdgeItemPtr> edges;
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested out edges.\n");
            return edges;
        }
        for (auto& edge_id: out_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false) {
                edges.emplace_back(edges_[edge_id]);
            }
        }
        return edges;
    }

    inline int64_t CountOutEdges(const NameType& node_name) const {
        int64_t num_edges;
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested out edges.\n");
            return num_edges;
        }
        for (auto& edge_id: out_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == true) {
                continue;
            }
            ++num_edges;
        }
        return num_edges;
    }

    inline std::vector<EdgeItemPtr> GetEdges(const NameType& node_name) const {
        /*
         * Gets all edges for a particular node, both in and out.
        */
        std::vector<EdgeItemPtr> edges;
        auto node_id = FindNodeId_(node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested out edges.\n");
            return edges;
        }
        for (auto& edge_id: out_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false) {
                edges.emplace_back(edges_[edge_id]);
            }
        }
        for (auto& edge_id: in_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false) {
                edges.emplace_back(edges_[edge_id]);
            }
        }
        return edges;
    }

    inline std::vector<EdgeItemPtr> GetEdges(const NameType& source_node_name, const NameType& sink_node_name) const {
        std::vector<EdgeItemPtr> edges;
        auto node_id = FindNodeId_(sink_node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested in edges.\n");
            return edges;
        }
        for (auto& edge_id: in_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false &&
                    edges_[edge_id]->source_name() == source_node_name) {
                edges.emplace_back(edges_[edge_id]);
            }
        }
        return edges;
    }

    inline bool HasEdge(const NameType& source_node_name, const NameType& sink_node_name) const {
        auto node_id = FindNodeId_(sink_node_name);
        if (node_id == GRAPH_UNDEFINED_NODE_ID) {
            LOG_ALL("Warning: Can not find the node for requested in edges.\n");
            return false;
        }
        for (auto& edge_id: in_edge_list_[node_id]) {
            if (edges_[edge_id]->is_removed() == false && edges_[edge_id]->source_name() == source_node_name) {
                return true;
            }
        }
        return false;
    }

    // bool RemoveEdge(const NameType& edge_name) {

    // }
    bool RemoveEdges(const NameType& source_node_name, const NameType& sink_node_name) {
        // Filter the in edges.
        {
            auto node_id = FindNodeId_(sink_node_name);
            if (node_id == GRAPH_UNDEFINED_NODE_ID) {
                LOG_ALL("Warning: Can not find the node for requested in edges.\n");
                return false;
            }

            EdgeListType new_in_edge_list;
            for (auto& edge_id: in_edge_list_[node_id]) {
                if (edges_[edge_id]->is_removed() == true) {
                    continue;
                }
                if (edges_[edge_id]->source_name() == source_node_name) {
                    edges_[edge_id]->is_removed(true);
                    continue;
                }
                new_in_edge_list.emplace_back(edge_id);
            }
            in_edge_list_[node_id] = new_in_edge_list;
        }

        // Filter the out edges.
        {
            auto node_id = FindNodeId_(source_node_name);
            if (node_id == GRAPH_UNDEFINED_NODE_ID) {
                LOG_ALL("Warning: Can not find the node for requested in edges.\n");
                return false;
            }

            EdgeListType new_out_edge_list;
            for (auto& edge_id: out_edge_list_[node_id]) {
                if (edges_[edge_id]->is_removed() == true) {
                    continue;
                }
                if (edges_[edge_id]->source_name() == sink_node_name) {
                    edges_[edge_id]->is_removed(true);
                    continue;
                }
                new_out_edge_list.emplace_back(edge_id);
            }
            out_edge_list_[node_id] = new_out_edge_list;
        }

        return false;

    }

    // std::shared_ptr<Graph<NameType, NodeDataType, EdgeDataType>> ExtractSubgraph(
    //                 const std::unordered_map<NameType, int8_t>& selected_nodes,
    //                 const std::unordered_map<NameType, int8_t>& selected_edges) {

    //     // Create the subgraph composed only of the selected nodes and edges.
    //     std::shared_ptr<Graph<NameType, NodeDataType, EdgeDataType>> subgraph = Graph<NameType, NodeDataType, EdgeDataType>::createGraph();

    //     for (const auto& node_item: this->nodes()) {
    //         if (node_item->is_removed()) {
    //             continue;
    //         }
    //         auto it = selected_nodes.find(node_item->name());
    //         if (it == selected_nodes.end()) {
    //             continue;
    //         }
    //         auto color = it->second;
    //         if (color == 0) {
    //             continue;
    //         }
    //         subgraph->AddNode(node_item->name(), node_item->data());
    //     }
    //     for (const auto& edge_item: this->edges()) {
    //         if (edge_item->is_removed()) {
    //             continue;
    //         }
    //         auto it = selected_edges.find(edge_item->name());
    //         if (it == selected_edges.end()) {
    //             continue;
    //         }
    //         auto color = it->second;
    //         if (color == 0) {
    //             continue;
    //         }

    //         // No need for anything fancy, just add it. Node names are identical as in the original graph.
    //         subgraph->AddEdge(edge_item->source_name(), edge_item->sink_name(), edge_item->data());
    //     }

    //     return subgraph;
    // }

    // inline std::vector<EdgeItemPtr> GetEdges(const NameType& source_node_name, const NameType& sink_node_name) const {
    //     std::vector<EdgeItemPtr> edges;
    //     // auto it = FindEdgesIterator_(source_node_name, sink_node_name);
    //     auto edge_key = MakeEdgeKey_(source_node_name, sink_node_name);
    //     auto it = edge_to_edge_ids_.find(edge_key);
    //     if (it == edge_to_edge_ids_.end()) {
    //         LOG_ALL("Warning: Can not find the node for requested edges.\n");
    //         return edges;
    //     }
    //     for (auto& edge_id: it->second) {
    //         edges.emplace_back(edges_[edge_id]);
    //     }
    //     return edges;
    // }

    // bool EdgeExists(const NameType& source_node_name, const NameType& sink_node_name) const {
    //     auto it = FindEdgesIterator_(source_node_name, sink_node_name);
    //     if (it == edge_to_edge_ids_.end()) {
    //         return false;
    //     }
    //     return true;
    // }

    //////////////////////////////////////////////
    ////////////// Graph algorithms //////////////
    //////////////////////////////////////////////

    /*
     * An implementation of Khan's algorithm
     * for topological sorting. Allow inheritance
     * for specialization to different graph types.
     * Returns false if the graph is not a DAG.
     * The sorted node names are reeturned via the only parameter.
    */
    virtual bool TopologicalSort(std::vector<NameType>& ret_sorted_node_names) const {
        std::vector<int64_t> sorted_nodes;
        ret_sorted_node_names.clear();

        std::deque<int64_t> q;                                      // Maintain the queue to process.
        std::unordered_map<int64_t, bool> is_edge_visited;         // Keep track of visited edges and
        std::vector<int64_t> in_edge_count(nodes_.size(), 0);       // processed in edges instead of deleting.
        int64_t num_processed_edges = 0;                            // Needed to detrmine if cyclic.
        int64_t num_nonremoved_nodes = 0;

        // Initialize the queue.
        for (size_t i = 0; i < nodes_.size(); ++i) {
            if (nodes_[i]->is_removed()) {
                continue;
            }
            ++num_nonremoved_nodes;
            // The in_edge_list_ should contain only non-filtered edges.
            int64_t num_in_edges = CountInEdges(nodes_[i]->name());
            if (num_in_edges == 0) {
                q.push_back((int64_t) i);
            }
        }

        // Empty graph, return true because it theoretically
        // already is fully sorted.
        if (num_nonremoved_nodes == 0) {
            return true;
        }

        // So, there are non-removed nodes in the graph, but
        // there are no nodes with zero in-degree. This means non-DAG.
        if (q.size() == 0) {
            return false;
        }

        int64_t num_non_removed_edges = 0;
        for (size_t i = 0; i < edges_.size(); ++i) {
            if (edges_[i]->is_removed()) {
                continue;
            }
            ++num_non_removed_edges;
        }

        while(q.size() > 0) {
            auto v = q.front();
            q.pop_front();
            sorted_nodes.emplace_back(v);

            for (auto& e: out_edge_list_[v]) {
                auto& edge_item = edges_[e];
                if (edge_item->is_removed()) {
                    continue;
                }
                auto w = edge_item->internal_sink_id(); // Get the destination node.
                // "Delete" the edge from the graph.
                if (is_edge_visited.find(e) != is_edge_visited.end()) {
                    continue;
                }
                is_edge_visited[e] = true;
                num_processed_edges += 1;
                in_edge_count[w] += 1;

                // Check if this node has any more in-edges. If not,
                // add it to the queue.
                if (in_edge_count[w] >= in_edge_list_[w].size()) {
                    q.push_back(w);
                }
            }
        }

        ret_sorted_node_names.clear();
        ret_sorted_node_names.reserve(sorted_nodes.size());

        // If true, the graph is actually cyclic!
        // We could throw, but as a library, we choose not to.
        if (num_processed_edges != num_non_removed_edges) {
            // LOG_ALL("Warning: Graph is cyclic, returning empty handed. num_processed_edges = %ld, num_non_removed_edges = %ld.\n", num_processed_edges, num_non_removed_edges);
            return false;
        }

        // Convert the internal node IDs to external names.
        for (auto& node_id: sorted_nodes) {
            auto& node_item = nodes_[node_id];
            ret_sorted_node_names.emplace_back(node_item->name());
        }

        return true;
    }

    bool CheckIsDAG() const {
        std::vector<NameType> sorted_nodes;
        bool is_dag = TopologicalSort(sorted_nodes);
        return is_dag;
    }

    bool CheckIsEmpty() const {
        if (nodes_.size() == 0 && edges_.size() == 0) {
            return true;
        }
        for (const auto& item: nodes_) {
            if (item->is_removed() == false) {
                return false;
            }
        }
        for (const auto& item: edges_) {
            if (item->is_removed() == false) {
                return false;
            }
        }
        return false;
    }

    std::vector<int64_t> FindSourceNodes() {
        return FindSourceOrSinkNodes_(false);
    }

    std::vector<int64_t> FindSinkNodes() {
        return FindSourceOrSinkNodes_(true);
    }



    //////////////////////////////////////////////
    /////////// Outputting and verbose ///////////
    //////////////////////////////////////////////

    virtual std::string ToJSON() const {
        std::ostringstream ss;
        ss << "[\"graph\", {\n"
            << "\"nitems\": [\n";
            // << "\t\"nitems\": [\n";
        for (size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
            // ss << "\t\t" << nodes_[node_id]->ToJSON();
            ss << nodes_[node_id]->ToJSON();
            if ((node_id + 1) < nodes_.size()) {
                ss << ",";
            }
            ss << "\n";
        }
        // ss << "\t],\n"
        //     << "\t\"eitems\": [\n";
        ss << "],\n"
            << "\"eitems\": [\n";
        for (size_t edge_id = 0; edge_id < edges_.size(); ++edge_id) {
            // ss << "\t\t" << edges_[edge_id]->ToJSON();
            ss << edges_[edge_id]->ToJSON();
            if ((edge_id + 1) < edges_.size()) {
                ss << ",";
            }
            ss << "\n";
        }
        // ss << "\t],\n"
        //     << "\t\"adj\": {\n";
        ss << "],\n"
            << "\"adj\": {";
        bool is_first_source = true;
        for (size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
            if (nodes_[node_id]->is_removed()) {
                continue;
            }
            if (is_first_source == false) {
                ss << ",";
            }
            // ss << "\n\t\t";
            ss << "\n";
            is_first_source = false;
            ss << nodes_[node_id]->name() << ": [";

            bool is_first_sink = true;
            for (auto& edge_id: out_edge_list_[node_id]) {
                if (edges_[edge_id]->is_removed() == true) {
                    continue;
                }
                if (is_first_sink == false) {
                    ss << ",";
                }
                is_first_sink = false;
                ss << edges_[edge_id]->sink_name();
            }
            ss << "]";
        }
        ss << "\n";

        ss << "}\n"
            << "}\n"
            << "]\n";
        // ss << "\t}\n"
        //     << "}\n"
        //     << "]\n";
        return ss.str();
    }

    /*
     * Allow writing the graph into the standard
     * DOT format which can be visualized.
    */
    virtual std::string GetDotString(std::string delimiter = " ") const {
        std::ostringstream oss;

        oss << "digraph raptor_graph {" << delimiter;

        for (size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
            if (nodes_[node_id]->is_removed()) {
                continue;
            }
            // const auto& node = nodes_[node_id]->data();
            oss << "    " << nodes_[node_id]->name() << ";" << delimiter; // Attributes can go in brackets, e.g. [label = "Test"].
        }

        for (size_t edge_id = 0; edge_id < edges_.size(); ++edge_id) {
            if (edges_[edge_id]->is_removed()) {
                continue;
            }
            oss << "    " << edges_[edge_id]->source_name() << " -> " << edges_[edge_id]->sink_name() << ";" << delimiter;
        }

        oss << "}";

        return oss.str();
    }

    /*
     * Writes the graph in the HTML format which
     * uses Vis.js to draw interactive nodes.
    */
    virtual void WriteAsHTML(std::ostream& oss) const {
        oss << "<!doctype html>" << "\n"
            << "    <html>" << "\n"
            << "        <head>" << "\n"
            << "            <title>Raptor Graph</title>" << "\n"
            << "            <script type=\"text/javascript\" src=\"http://visjs.org/dist/vis.js\"></script>" << "\n"
            << "        </head>" << "\n"
            << "        <body>" << "\n";

        // Define the CSS style.
        oss << "            <style media=\"screen\" type=\"text/css\">" << "\n"
            << "            " << "\n"
            << "                code {" << "\n"
            << "                    background: #2db34a;" << "\n"
            << "                    border-radius: 6px;" << "\n"
            << "                    color: #fff;" << "\n"
            << "                    display: block;" << "\n"
            << "                    font: 14px/24px \"Source Code Pro\", Inconsolata, \"Lucida Console\", Terminal, \"Courier New\", Courier;" << "\n"
            << "                    padding: 24px 15px;" << "\n"
            << "                    text-align: center;" << "\n"
            << "                }" << "\n"
            << "                header," << "\n"
            << "                section," << "\n"
            << "                aside," << "\n"
            << "                footer {" << "\n"
            << "                    margin: 0 1.5% 24px 1.5%;" << "\n"
            << "                }" << "\n"
            << "                section {" << "\n"
            << "                    float: left;" << "\n"
            << "                    width: 63%;" << "\n"
            << "                }" << "\n"
            << "                aside {" << "\n"
            << "                    float: right;" << "\n"
            << "                    width: 30%;" << "\n"
            << "                }" << "\n"
            << "                footer {" << "\n"
            << "                    clear: both;" << "\n"
            << "                    margin-bottom: 0;" << "\n"
            << "                }" << "\n"
            << "            " << "\n"
            << "            </style>" << "\n";

        oss << "            <aside>" << "\n";
        for (size_t node_id = 0; node_id < nodes_.size(); ++node_id) {
            if (nodes_[node_id]->is_removed()) {
                continue;
            }
            oss << "           Node: '" << nodes_[node_id]->name() << " = '" << nodes_[node_id]->data()->Verbose() << "'<br>\n";
        }
        oss << "            </aside>\n";

        std::string dot = GetDotString();

        oss
            << "            <section>" << "\n"
            << "            <div id=\"raptor_graph_div\"></div>" << "\n"
            << "            <script type=\"text/javascript\">" << "\n"
            << "            dot_string = '" << dot << "';" << "\n"
            << "            parsed_data = vis.network.convertDot(dot_string);" << "\n"
            << "            var container = document.getElementById('raptor_graph_div');" << "\n"
            << "            var data= {" << "\n"
            << "                nodes: parsed_data.nodes," << "\n"
            << "                edges: parsed_data.edges," << "\n"
            << "            };" << "\n"
            << "            var options = {" << "\n"
            << "                width: '100%'," << "\n"
            << "                height: '800px'" << "\n"
            << "            };" << "\n"
            << "            var network = new vis.Network(container, data, options);" << "\n"
            << "            </script>" << "\n"
            << "            </section>" << "\n"
            << "        </body>" << "\n"
            << "    </html>" << "\n";
    }

    virtual std::string Verbose() {
        std::ostringstream oss;
        oss << "Graph info:" << std::endl;
        oss << "  Nodes: " << nodes_.size() << std::endl;
        oss << "  Edges: " << edges_.size() << std::endl;
        oss << std::endl;
        oss << "  All nodes:" << std::endl;
        for (size_t i = 0; i < nodes_.size(); i++) {
            auto& node_item = nodes_[i];
            if (node_item->is_removed()) {
                continue;
            }
            auto node_id = node_item->name();
            auto& node = node_item->data();
            oss << "    [i = " << i << "] Node: name = " << node_id << ", is_removed = " << static_cast<int32_t>(node_item->is_removed()) << std::endl;
            oss << "            Node data: " << node->Verbose() << std::endl;
        }
        oss << std::endl;
        oss << "  All edges:" << std::endl;
        for (size_t i = 0; i < edges_.size(); i++) {
            auto& edge_item = edges_[i];
            if (edge_item->is_removed()) {
                continue;
            }
            auto source_node_id = edge_item->internal_source_id();
            auto sink_node_id = edge_item->internal_sink_id();
            auto& edge = edge_item->data();
            oss << "    [i = " << i << "] Edge: source_name = " << edge_item->source_name() << ", sink_name = " << edge_item->sink_name()
                << "; (source_node_id =  " << source_node_id << ", sink_node_id = " << sink_node_id
                << ", is_removed = " << static_cast<int32_t>(edge_item->is_removed()) << std::endl;
            oss << "            Edge data: " << edge->Verbose() << std::endl;
        }

        return oss.str();
    }

    // virtual std::string Verbose() {
    //     std::ostringstream oss;
    //     oss << "Graph info:" << std::endl;
    //     oss << "  Nodes: " << nodes_.size() << std::endl;
    //     oss << "  Edges: " << edges_.size() << std::endl;
    //     return oss.str();
    // }

    const std::vector<NodeItemPtr>& nodes() {
        return nodes_;
    }

    const std::vector<EdgeItemPtr>& edges() {
        return edges_;
    }

    virtual int64_t CountEdges() const {
        int64_t ret = 0;
        for (const auto& e: edges_) {
            if (e->is_removed()) {
                continue;
            }
            ++ret;
        }
        return ret;
    }

    virtual void Clear() {
        nodes_.clear();
        edges_.clear();
        in_edge_list_.clear();
        out_edge_list_.clear();
        name_to_node_id_.clear();
    }

protected:
    Graph() {

    }

// private:
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;

    // Vector of all nodes.
    std::vector<NodeItemPtr> nodes_;
    // Vector of all edges.
    std::vector<EdgeItemPtr> edges_;
    // List of in-edges for every node, in the same order.
    std::vector<EdgeListType> in_edge_list_;
    // List of in-edges for every node, in the same order.
    std::vector<EdgeListType> out_edge_list_;

    // Relation map from node name to the internal ID of the node.
    std::unordered_map<NameType, uint64_t> name_to_node_id_;

    /*
     * Relation map between an edge encoded by source and sink node ID, and the actual
     * objects in the edges_ vector. Multiedges between two nodes are supported.
    */
    // std::unordered_map<__int128, std::vector<uint64_t>> edge_to_edge_ids_;

    /*
     * Returns the node ID if there is a node with the specified name, otherwise
     * returns GRAPH_UNDEFINED_NODE_ID value.
    */
    uint64_t FindNodeId_(const NameType& node_name) const {
        auto it = name_to_node_id_.find(node_name);
        if (it == name_to_node_id_.end()) {
            return GRAPH_UNDEFINED_NODE_ID;
        }
        return it->second;
    }

    static inline __int128 FormatEdgeKey_(uint64_t source_node_id, uint64_t sink_node_id) {
        __int128 edge_key = (((__int128) source_node_id) << 64) | (((__int128) sink_node_id) & GRAPH_EDGE_KEY_64bit_MASK);
        return edge_key;
    }

    /*
     * Supports multiedges.
    */
    inline __int128 MakeEdgeKey_(const NameType& source_node_name, const NameType& sink_node_name) const {
        // Node ID will be GRAPH_UNDEFINED_NODE_ID if the node name was not found.
        auto source_node_id = FindNodeId_(source_node_name);
        auto sink_node_id = FindNodeId_(sink_node_name);
        return FormatEdgeKey_(source_node_id, sink_node_id);
    }

    std::vector<int64_t> FindSourceOrSinkNodes_(bool find_sinks) {
        std::vector<int64_t> ret;
        for (const auto& node_item: nodes()) {
            if (node_item->is_removed()) {
                continue;
            }
            auto num_edges = (find_sinks == false) ? CountInEdges(node_item->name()) : CountOutEdges(node_item->name());
            if (num_edges == 0) {
                ret.emplace_back(node_item->name());
            }
        }
        return ret;
    }



};

// /*
//  * MakeGraphShallowCopy creates new NodeItem and EdgeItem objects and copies their original values,
//  * but the underlying data() points to the same object as the original graph.
// */
// template<class NameType, class NodeDataType, class EdgeDataType>
// std::shared_ptr<raptor::Graph<NameType, NodeDataType, EdgeDataType>> MakeGraphShallowCopy(const std::shared_ptr<raptor::Graph<NameType, NodeDataType, EdgeDataType>>& in_graph) {
//     std::shared_ptr<raptor::Graph<NameType, NodeDataType, EdgeDataType>> out_graph;

//     std::vector<int64_t> node_id_after_removal(in_graph->nodes().size(), -1);
//     size_t num_removed = 0;
//     for (size_t i = 0; i < in_graph->nodes().size(); ++i) {
//         const auto& v_item = in_graph->nodes()[i];
//         if (v_item->is_removed()) {
//             ++num_removed;
//             continue;
//         }
//         node_id_after_removal[i] = static_cast<int64_t>(i - num_removed);


//         out_graph->AddNode(v_item->name(), new_v_item);
//     }

// }

}

#endif
