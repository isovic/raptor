/*
 * segment_edge.cc
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#include <graph/segment_edge.h>
#include <sstream>

namespace raptor {

SegmentEdgePtr createSegmentEdge() {
    return SegmentEdgePtr(new raptor::SegmentEdge());
}

raptor::SegmentEdge::SegmentEdge()
                :   symbolic_edge_name_(),
                    id_(-1),
                    source_name_(),
                    source_is_rev_(false),
                    source_id_(0),
                    sink_name_(),
                    sink_is_rev_(false),
                    sink_id_(0),
                    source_start_(0),
                    source_end_(0),
                    source_len_(0),
                    sink_start_(0),
                    sink_end_(0),
                    sink_len_(0),
                    alignment_(),
                    tags_() {

}

std::string raptor::SegmentEdge::Verbose() const {
    std::ostringstream oss;

    // oss << "symbolic_edge_name = " << symbolic_edge_name_ << ", "
    //     << "id = " << id_ << ", "
    //     << "source_name = " << source_name_ << ", "
    //     << "source_id = " << source_id_ << ", "
    //     << "source_is_rev = " << source_is_rev_ << ", "
    //     << "sink_name = " << sink_name_ << ", "
    //     << "sink_id = " << sink_id_ << ", "
    //     << "sink_is_rev = " << sink_is_rev_ << ", "
    //     << "source_start = " << source_start_ << ", "
    //     << "source_end = " << source_end_ << ", "
    //     << "source_len = " << source_len_ << ","
    //     << "sink_start = " << sink_start_ << ", "
    //     << "sink_end = " << sink_end_ << ", "
    //     << "sink_len = " << sink_len_;

    oss << "name_str = '" << symbolic_edge_name_ << "', "
        << "id = " << id_ << ", "
        << "source = {name:'" << source_name_ << "', s:" << source_start_ << ", e:" << source_end_ << ", id:" << source_id_ << ", len:" << source_len_ << ", rev:" << source_is_rev_ << "}"
        << ", "
        << "sink = {name:'" << sink_name_ << "', s:" << sink_start_ << ", e:" << sink_end_ << ", id:" << sink_id_ << ", len:" << sink_len_ << ", rev:" << sink_is_rev_ << "}";

    return oss.str();
}

std::string raptor::SegmentEdge::ToDOT() const {
    return std::string();
}

std::string raptor::SegmentEdge::ToJSON() const {
    std::ostringstream ss;
    ss << "{"
        << "\"symname\":\"" << symbolic_edge_name() << "\""
        << ",\"id\":" << id()
        << ",\"v\":" << source_name()
        << ",\"v_int_id\":" << source_id()
        << ",\"vr\":" << source_is_rev()
        << ",\"vs\":" << source_start()
        << ",\"ve\":" << source_end()
        << ",\"vl\":" << source_len()

        << ",\"w\":" << sink_name()
        << ",\"w_int_id\":" << sink_id()
        << ",\"wr\":" << sink_is_rev()
        << ",\"ws\":" << sink_start()
        << ",\"we\":" << sink_end()
        << ",\"wl\":" << sink_len()
        << "}";
    return ss.str();
}

}
