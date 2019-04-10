/*
 * edge_item.hpp
 *
 *  Created on: Jan 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_EDGE_ITEM_H_
#define SRC_GRAPH_EDGE_ITEM_H_

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace raptor {

template<class NameType, class NodeDataType>
class NodeItem;

template<class NameType, class EdgeDataType>
class EdgeItem {
public:
    // friend std::shared_ptr<EdgeItem<NameType, EdgeDataType>> make_shared();

    EdgeItem(int64_t _internal_id, const NameType& _name, const NameType& _source_name, const NameType& _sink_name,
             int64_t _internal_source_id, int64_t _internal_sink_id,
             const std::shared_ptr<EdgeDataType> _data)
                            :   internal_id_(_internal_id),
                                name_(_name),
                                source_name_(_source_name),
                                sink_name_(_sink_name),
                                internal_source_id_(_internal_source_id),
                                internal_sink_id_(_internal_sink_id),
                                data_(_data),
                                is_removed_(false) { }

    static std::shared_ptr<EdgeItem<NameType, EdgeDataType>> createEdgeItem(
                                                                int64_t _internal_id,
                                                                const NameType& _name,
                                                                const NameType& _source_name,
                                                                const NameType& _sink_name,
                                                                int64_t _internal_source_id,
                                                                int64_t _internal_sink_id,
                                                                const std::shared_ptr<EdgeDataType> _data) {
        return std::make_shared<EdgeItem<NameType, EdgeDataType>>(_internal_id, _name, _source_name, _sink_name, _internal_source_id, _internal_sink_id, _data);
    }

    std::string ToDOT() const {
        return "";
    }
    std::string ToJSON() const {
        std::ostringstream ss;
        ss << "[\"eitem\", {\"name\": " << name()
            << ", \"removed\": " << is_removed()
            << ", \"int_id\": " << internal_id()
            << ", \"v\": " << source_name()
            << ", \"w\": " << sink_name()
            << ", \"v_int_id\": " << internal_source_id()
            << ", \"w_int_id\": " << internal_sink_id()
            << ", \"data\": " << data()->ToJSON()
            << "}]";
        return ss.str();
    }

    int64_t internal_id() const {
        return internal_id_;
    }
    const NameType& name() const {
        return name_;
    }
    const NameType& source_name() const {
        return source_name_;
    }
    const NameType& sink_name() const {
        return sink_name_;
    }
    int64_t internal_source_id() const {
        return internal_source_id_;
    }
    int64_t internal_sink_id() const {
        return internal_sink_id_;
    }
    const std::shared_ptr<EdgeDataType> data() const {
        return data_;
    }

    bool is_removed() const {
        return is_removed_;
    }
    void is_removed(bool val) {
        is_removed_ = val;
    }

private:
    EdgeItem(const EdgeItem&) = delete;
    EdgeItem& operator=(const EdgeItem&) = delete;

    int64_t internal_id_;
    NameType name_;
    NameType source_name_;
    NameType sink_name_;
    int64_t internal_source_id_;
    int64_t internal_sink_id_;
    const std::shared_ptr<EdgeDataType> data_;
    bool is_removed_;
};

}

#endif
