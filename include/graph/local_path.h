/*
 * local_path.h
 *
 *  Created on: Feb 1, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_LOCAL_PATH_H_
#define SRC_RAPTOR_LOCAL_PATH_H_

#include <graph/anchor_graph_edge.h>
#include <algorithm/graph.hpp>
#include <algorithm/path.hpp>
#include <containers/region/region_mapped.h>
#include <cstdint>
#include <string>

namespace raptor {

typedef raptor::Path<int64_t, raptor::RegionMapped, raptor::AnchorGraphEdge> LocalPath;

typedef std::shared_ptr<raptor::NodeItem<int64_t, raptor::RegionMapped>> LocalPathNodeItemPtr;

class LocalPathTools {
    public:
    /*
     * This method generates a string representation of the full path, where all
     * nodes and edges are listed. Every edge is described as a list of segment edge IDs,
     * e.g. "(1=2=83)"
     * If an edge is implicit, it will be marked as "(i)"
     * The output will be in the form of:
     *  1(i)2(i)3(1=2=3)4(5)6(i)7(i)8(i)9(i)10
     * Explanation:
     *   Nodes and edges are listed interleaved. Each edge is surrounded with matching
     *   brackets. There cannot be consecutive edges listed in the string. Composite
     *   edges composed of several successive segment edges are listed in the same
     *   brackets, but separated with a '=' character.
    */
    static std::string SummarizeFullPathAsString(const std::shared_ptr<raptor::LocalPath> path);

    /*
     * This method generates a string representation of the path, but unlike the
     * SummarizeFullPathAsString method, this method does not summarize nodes
     * between the LocalEdges.
     * Reason: Let N be a node with one input LocalEdge E1 and one output local edge E2.
     * A path which enters node N through E1 and exits it through E2 will spell out the
     * full sequence between E1 and E2, regardless of the nodes in between E1 and E2.
     * The first node in the path and the last node in the path will be listed, and
     * other nodes ignored.
     * The output will be in the form of:
     *  1(1)(2)(3)(5)10
     * In this case - the string can have several consecutive edges listed.
    */
    static std::string SummarizeJumpPathAsString(const std::shared_ptr<raptor::LocalPath> path);

    /*
     * Creates a new LocalPath with neighboring nodes, which are connected exclusively by implicit edges,
     * are merged together into larger anchors (mapped regions).
     * This functionality is needed for e.g. writing to PAF.
    */
    static std::shared_ptr<raptor::LocalPath> MergeImplicitEdges(const std::shared_ptr<raptor::LocalPath> path);

};

}

#endif
