/*
 * segment_node.cc
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#include <graph/segment_node.h>
#include <sstream>

namespace raptor {

// SegmentNodePtr createSegmentNode() {
//     return SegmentNodePtr(new raptor::SegmentNode());
// }

// SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len, bool _is_rev) {
//     return SegmentNodePtr(new raptor::SegmentNode(_name, _len, _is_rev));
// }

SegmentNodePtr createSegmentNode(const std::string& _name, int64_t _len, bool _is_rev, int64_t _seq_id) {
    return SegmentNodePtr(new raptor::SegmentNode(_name, _len, _is_rev, _seq_id));
}

std::string raptor::SegmentNode::Verbose() const {
    return std::string();
}

std::string raptor::SegmentNode::ToDOT() const {
    return std::string();
}

std::string raptor::SegmentNode::ToJSON() const {
    std::ostringstream ss;
        ss << "{\"seq_id\":" << seq_id()
            << ",\"len\":" << len()
            << ",\"rev\":" << is_rev()
            << ",\"score\":" << score()
            << ",\"header\":\"" << header() << "\""
            << ",\"int_id\":" << internal_id()
            << "}";
    return ss.str();
}

}
