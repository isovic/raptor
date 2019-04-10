/*
 * segment_node.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEGMENT_NODE_H_
#define SRC_SEGMENT_NODE_H_

#include <cstdint>
#include <memory>
#include <string>
#include <algorithm/node_data_base.hpp>

namespace raptor {

class SegmentNode;

typedef std::shared_ptr<raptor::SegmentNode> SegmentNodePtr;

// SegmentNodePtr createSegmentNode();
// SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len_, bool _is_rev);
SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len, bool _is_rev, int64_t _seq_id);

class SegmentNode : public NodeDataBase {
public:
    // friend SegmentNodePtr createSegmentNode();
    // friend SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len, bool _is_rev);
    friend SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len, bool _is_rev, int64_t _seq_id);
    ~SegmentNode() = default;

    std::string Verbose() const;

    std::string ToDOT() const;
    std::string ToJSON() const;

    const std::string& header() const { return header_; }
    int64_t len() const { return len_; }
    bool is_rev() const { return is_rev_; }
    int64_t seq_id() const { return seq_id_; }
    int64_t internal_id() const { return internal_id_; }
    float score() const { return score_;}

    void header(const std::string& val) { header_ = val; }
    void len(int64_t val) { len_ = val; }
    void is_rev(bool val) { is_rev_ = val; }
    void seq_id(int64_t val) { seq_id_ = val; }
    void internal_id(int64_t val) { internal_id_ = val; }
    void score(float val) { score_ = val; }

private:
    SegmentNode()
        :   header_(), len_(0), is_rev_(false), seq_id_(-1), internal_id_(-1),
            score_(0.0f)
    { }

    SegmentNode(const std::string& _header, int64_t _len, bool _is_rev)
        :   header_(_header), len_(_len), is_rev_(_is_rev), seq_id_(-1), internal_id_(-1),
            score_(0.0f)
    { }

    SegmentNode(const std::string& _header, int64_t _len, bool _is_rev, int64_t _seq_id)
        :   header_(_header), len_(_len), is_rev_(_is_rev), seq_id_(_seq_id), internal_id_(-1),
            score_(0.0f)
    { }

    SegmentNode(const SegmentNode&) = delete;
    SegmentNode& operator=(const SegmentNode&) = delete;

    std::string header_;
    int64_t len_;
    bool is_rev_;
    int64_t seq_id_;
    int64_t internal_id_;
    float score_;
};

}

#endif
