/*
 * local_graph_node.h
 *
 *  Created on: Jan 27, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_LOCAL_GRAPH_NODE_H_
#define SRC_RAPTOR_LOCAL_GRAPH_NODE_H_

#include <memory>
#include <sstream>
#include <vector>
#include <types/typedefs.h>
#include <algorithm/node_data_base.hpp>

namespace raptor {

class AnchorGraphNode;

std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor);
std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor, double _score, int64_t _predecessor, int32_t _path_id);

/*
 * Each edge in the AnchorGraph can actually contain multiple successive
 * edges from the SegmentGraph.
 * Interpretation: transcriptome mapping, but some small exons do not have
 * any hits. We need to be able to align through them.
*/
class AnchorGraphNode : public NodeDataBase {
public:
    friend std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor);
    friend std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor, double _score, int64_t _predecessor, int32_t _path_id);

    ~AnchorGraphNode() = default;

    std::string ToDOT() const {
        return std::string();
    }

    std::string ToJSON() const {
        std::ostringstream ss;
        ss << "{\"score\":" << score()
            << ",\"predecessor\":" << predecessor()
            << ",\"path_id\":" << path_id()
            << ",\"anchor\":\"" << anchor()->WriteAsCSV(',') << "\""
            << "}";
        return ss.str();
    }

    double score() const {
        return score_;
    }
    raptor::AnchorPtr anchor() const {
        return anchor_;
    }
    int64_t predecessor() const {
        return predecessor_;
    }
    int32_t path_id() const {
        return path_id_;
    }

    void score(double val) {
        score_ = val;
    }
    void predecessor(int64_t val) {
        predecessor_ = val;
    }
    void path_id(int32_t val) {
        path_id_ = val;
    }
    void anchor(const raptor::AnchorPtr& val) {
        anchor_ = val;
    }

    std::string Verbose() const {
        std::ostringstream oss;
        oss << "AnchorGraphNode: score = " << score_ << ", predecessor = " << predecessor_ << ", path_id = " << path_id_ << ":=> ";
        oss << anchor_->Verbose();
        return oss.str();
    }

private:
    AnchorGraphNode(const raptor::AnchorPtr& _anchor)
        : anchor_(_anchor), score_(0.0),
            predecessor_(-1), path_id_(-1) { }
    AnchorGraphNode(const raptor::AnchorPtr& _anchor, double _score, int64_t _predecessor, int32_t _path_id)
        : anchor_(_anchor), score_(_score),
            predecessor_(_predecessor), path_id_(_path_id) { }

    raptor::AnchorPtr anchor_;
    double score_;
    int64_t predecessor_;
    int32_t path_id_;
};

inline std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor) {
    return std::shared_ptr<AnchorGraphNode>(new AnchorGraphNode(_anchor));
}

inline std::shared_ptr<AnchorGraphNode> createAnchorGraphNode(const raptor::AnchorPtr& _anchor, double _score, int64_t _predecessor, int32_t _path_id) {
    return std::shared_ptr<AnchorGraphNode>(new AnchorGraphNode(_anchor, _score, _predecessor, _path_id));
}

}

#endif
