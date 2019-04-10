/*
 * segment_graph_parser.h
 *
 *  Created on: Dec 17, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEGMENT_GRAPH_PARSER_H_
#define SRC_SEGMENT_GRAPH_PARSER_H_

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>

#include <graph/segment_graph.h>
#include <index/minimizer_index.h>
#include <types/typedefs.h>
#include <graph/segment_graph_parser_enums.h>

namespace raptor {

typedef std::shared_ptr<raptor::SegmentGraph> SegmentGraphPtr;

class GraphLoader {
public:
    static std::shared_ptr<raptor::SegmentGraph> FromGFA(const std::string& graph_path, mindex::IndexPtr index, bool add_symmetric_arcs);
    static void AssertGraphValidity(const mindex::IndexPtr& index, const std::shared_ptr<raptor::SegmentGraph>& graph);
    /*
     * For every sequence in the index, this method adds a single node to
     * the segment graph. It also creates a query_name->query_id lookup.
    */
    static void AddGraphNodesFromIndex(
                            mindex::IndexPtr index, std::shared_ptr<raptor::SegmentGraph> graph);

    static raptor::SegmentEdgePtr ParseGFA1Edge(const std::shared_ptr<raptor::SegmentGraph>& graph,
                                                        const std::vector<std::string>& params);

    static raptor::SegmentEdgePtr ParseGFA2Edge(const std::shared_ptr<raptor::SegmentGraph>& graph,
                                                        const std::vector<std::string>& params);

private:

    /*
     * For each edge in the graph, this method adds it's reverse complement.
    */
    static void AddSymmetricArcs_(std::shared_ptr<raptor::SegmentGraph> graph);

    static float FindGFAVersion_(std::ifstream& ifs);

    static void LoadEdgesFromGFA1_(
                            std::ifstream& ifs,
                            mindex::IndexPtr index,
                            std::shared_ptr<raptor::SegmentGraph> graph);
    static void LoadEdgesFromGFA2_(
                            std::ifstream& ifs,
                            mindex::IndexPtr index,
                            std::shared_ptr<raptor::SegmentGraph> graph);
};

}

#endif
