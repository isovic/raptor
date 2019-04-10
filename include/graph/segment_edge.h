/*
 * segment_edge.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEGMENT_EDGE_H_
#define SRC_SEGMENT_EDGE_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <algorithm/edge_data_base.hpp>

namespace raptor {

class SegmentEdge;

typedef std::shared_ptr<raptor::SegmentEdge> SegmentEdgePtr;

SegmentEdgePtr createSegmentEdge();

class SegmentEdge : public EdgeDataBase {
public:
    friend SegmentEdgePtr createSegmentEdge();
    ~SegmentEdge() = default;

    std::string Verbose() const;

    std::string ToDOT() const;
    std::string ToJSON() const;

    // Getters
    const std::string& symbolic_edge_name() const {
        return symbolic_edge_name_;
    }
    int64_t id() const {
        return id_;
    }
    const std::string& source_name() const {
        return source_name_;
    }
    int64_t source_id() const {
        return source_id_;
    }
    const std::string& sink_name() const {
        return sink_name_;
    }
    int64_t sink_id() const {
        return sink_id_;
    }
    bool source_is_rev() const {
        return source_is_rev_;
    }
    bool sink_is_rev() const {
        return sink_is_rev_;
    }
    int64_t source_start() const {
        return source_start_;
    }
    int64_t source_end() const {
        return source_end_;
    }
    int64_t source_len() const {
        return source_len_;
    }
    int64_t sink_start() const {
        return sink_start_;
    }
    int64_t sink_end() const {
        return sink_end_;
    }
    int64_t sink_len() const {
        return sink_len_;
    }
    const std::string& alignment() const {
        return alignment_;
    }
    const std::vector<std::string>& tags() const {
        return tags_;
    }

    // Setters.
    void symbolic_edge_name(const std::string& _symbolic_edge_name) {
        symbolic_edge_name_ = _symbolic_edge_name;
    }
    void id(int64_t _id) {
        id_ = _id;
    }
    void source_name(const std::string& _source_name) {
        source_name_ = _source_name;
    }
    void source_id(int64_t _source_id) {
        source_id_ = _source_id;
    }
    void sink_name(const std::string& _sink_name) {
        sink_name_ = _sink_name;
    }
    void sink_id(int64_t _sink_id) {
        sink_id_ = _sink_id;
    }
    void source_is_rev(bool _source_is_rev) {
        source_is_rev_ = _source_is_rev;
    }
    void sink_is_rev(bool _sink_is_rev) {
        sink_is_rev_ = _sink_is_rev;
    }
    void source_start(int64_t _source_start) {
        source_start_ = _source_start;
    }
    void source_end(int64_t _source_end) {
        source_end_ = _source_end;
    }
    void source_len(int64_t _source_len) {
        source_len_ = _source_len;
    }
    void sink_start(int64_t _sink_start) {
        sink_start_ = _sink_start;
    }
    void sink_end(int64_t _sink_end) {
        sink_end_ = _sink_end;
    }
    void sink_len(int64_t _sink_len) {
        sink_len_ = _sink_len;
    }
    void alignment(const std::string& _alignment) {
        alignment_ = _alignment;
    }
    void tags(const std::vector<std::string>& _tags) {
        tags_ = _tags;
    }

private:
    SegmentEdge();
    SegmentEdge(const SegmentEdge&) = delete;
    SegmentEdge& operator=(const SegmentEdge&) = delete;

    std::string symbolic_edge_name_;
    int64_t id_;    // ID of this SegmentEdge object. Needs to be unique.

    std::string source_name_;
    bool source_is_rev_;
    int64_t source_id_;
    std::string sink_name_;
    bool sink_is_rev_;
    int64_t sink_id_;

    int64_t source_start_;
    int64_t source_end_;
    int64_t source_len_;

    int64_t sink_start_;
    int64_t sink_end_;
    int64_t sink_len_;

    std::string alignment_;
    std::vector<std::string> tags_;

};

}

#endif
