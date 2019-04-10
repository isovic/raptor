/*
 * anchor_graph_edge.h
 *
 *  Created on: Dec 22, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_LOCAL_GRAPH_EDGE_H_
#define SRC_RAPTOR_LOCAL_GRAPH_EDGE_H_

#include <memory>
#include <sstream>
#include <vector>
#include <graph/segment_edge.h>
#include <log/log_tools.h>
#include <algorithm/edge_data_base.hpp>

namespace raptor {

class AnchorGraphEdge;

using AnchorGraphEdgePtr = std::shared_ptr<AnchorGraphEdge>;
std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(SegmentEdgePtr _segment_edge);
std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges);
std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges, double _score, int32_t _path_id);
std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse);
std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse, double _score, int32_t _path_id);

/*
 * Each edge in the AnchorGraph can actually contain multiple successive
 * edges from the SegmentGraph.
 * Interpretation: transcriptome mapping, but some small exons do not have
 * any hits. We need to be able to align through them.
*/
class AnchorGraphEdge : public EdgeDataBase {
public:
    friend std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(SegmentEdgePtr _segment_edge);
    friend std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges);
    friend std::shared_ptr<AnchorGraphEdge> createAnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges, double _score, int32_t _path_id);
    friend std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse);
    friend std::shared_ptr<AnchorGraphEdge> createImplicitAnchorGraphEdge(bool is_reverse, double _score, int32_t _path_id);

    ~AnchorGraphEdge() = default;

    std::string ToDOT() const;
    std::string ToJSON() const;

    /*
     * In case there are multiple segment_edges_, their junction
     * locations can be separated by a certain amount of bases.
     * This method returns the sum of target bases between all
     * segment_edges_.
     * It returns 0 if there is <= segment edge.
    */
    inline int64_t CalcJumpLength() const {
        int64_t jump_length = 0;
        for (size_t i = 1; i < segment_edges_.size(); i++) {
            // Neighboring edges need to link the same segment.
            if (segment_edges_[i-1]->sink_id() != segment_edges_[i]->source_id()) {
                LOG_ALL("Warning: incorrect definition of a AnchorGraphEdge! Neighboring segment edges do not link same nodes.\n");
                return 0;
            }
            jump_length += (segment_edges_[i]->source_start() - segment_edges_[i-1]->sink_end());
        }
        return jump_length;
    }

    std::string SummarizeSegmentEdgesAsString() const {
        /*
         * Returns a JSON formatted list of segment edges.
        */
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < segment_edges_.size(); i++) {
            if (i > 0) {
                oss << ",";
            }
            oss << "'" << segment_edges_[i]->symbolic_edge_name() << "'";
        }
        oss << "]";
        return oss.str();
    }

    /*
     * Implicit edges would be those which link two chains/anchors
     * on the same segment. E.g. anchors which would form a chain on the same
     * segment, or chains which would form a path on the same segment.
     * If any of the SegmentEdgePtr objects is nullptr, or the segment_edges_
     * vector is empty, the edge is considered implicit.
    */
    inline bool IsImplicit() const {
        return is_implicit_ || is_implicit_reverse_;
        // if (segment_edges_.size() == 0) {
        //     return true;
        // }
        // for (auto& edge: segment_edges_) {
        //     if (edge == nullptr) {
        //         return true;
        //     }
        // }
        // return false;
    }

    /* Start position of the source junction point of the first edge.
     * Returns < 0 if no segment_edges_ exist.
    */
    inline int64_t source_start() const {
        if (IsImplicit()) { return -1; }
        return segment_edges_.front()->source_start();
    }

    /* End position of the source junction point of the first edge.
     * Returns < 0 if no segment_edges_ exist.
    */
    inline int64_t source_end() const {
        if (IsImplicit()) { return -1; }
        return segment_edges_.front()->source_end();
    }

    /* Start position of the sink junction point of the last edge.
     * Returns < 0 if no segment_edges_ exist.
    */
    inline int64_t sink_start() const {
        if (IsImplicit()) { return -1; }
        return segment_edges_.back()->sink_start();
    }

    /* End position of the sink junction point of the last edge.
     * Returns < 0 if no segment_edges_ exist.
    */
    inline int64_t sink_end() const {
        if (IsImplicit()) { return -1; }
        return segment_edges_.back()->sink_end();
    }

    double score() const {
        return score_;
    }
    int32_t path_id() const {
        return path_id_;
    }

    inline const std::vector<SegmentEdgePtr> segment_edges() const {
        return segment_edges_;
    }
    bool is_implicit() const {
        return is_implicit_;
    }
    bool is_implicit_reverse() const {
        return is_implicit_reverse_;
    }
    void is_implicit(bool _is_implicit) {
        is_implicit_ = _is_implicit;
    }
    void is_implicit_reverse(bool _is_implicit_reverse) {
        is_implicit_reverse_ = _is_implicit_reverse;
    }
    void score(double val) {
        score_ = val;
    }
    void path_id(int32_t val) {
        path_id_ = val;
    }

    std::string Verbose() const {
        std::ostringstream oss;
        int64_t jump_length = CalcJumpLength();
        oss << "AnchorGraphEdge: is_implicit = " << is_implicit() << ", is_implicit_reverse = " << is_implicit_reverse() <<
                ", score = " << score_ << ", path_id = " << path_id_ <<
                ", num_segment_edges = " << segment_edges_.size() << ", jump_len = " << jump_length;
        if (segment_edges_.size() > 0) {
            oss << ":" << std::endl;
            for (size_t i = 0; i < segment_edges_.size(); i++) {
                auto& edge = segment_edges_[i];
                if (edge != nullptr) {
                    oss << "  " << edge->Verbose() << std::endl;
                } else {
                    oss << "  (implicit edge)" << std::endl;
                }
            }
        }

        oss << std::endl;

        return oss.str();
    }

private:
    /*
     * Initialize the LocalEdge with a single SegmentEdge.
     * This is therefore not an implicit edge.
    */
    AnchorGraphEdge(SegmentEdgePtr _segment_edge)
        : segment_edges_{_segment_edge}, is_implicit_(false), is_implicit_reverse_(false), score_(0.0), path_id_(-1) {

    }

    /*
     * Initialize the LocalEdge with a vector of segment edges.
     * This is therefore not an implicit edge.
    */
    AnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges)
        : segment_edges_{_segment_edges}, is_implicit_(false), is_implicit_reverse_(false), score_(0.0), path_id_(-1) {

    }

    AnchorGraphEdge(const std::vector<SegmentEdgePtr>& _segment_edges, bool _is_implicit, bool _is_implicit_reverse, double _score, int32_t _path_id)
        : segment_edges_{_segment_edges}, is_implicit_(_is_implicit), is_implicit_reverse_(_is_implicit_reverse), score_(_score), path_id_(_path_id) {

    }

    std::vector<SegmentEdgePtr> segment_edges_;
    bool is_implicit_;
    bool is_implicit_reverse_;
    double score_;
    int32_t path_id_;
};



}

#endif
