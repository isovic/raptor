/*
 * node_item.hpp
 *
 *  Created on: Jan 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_NODE_ITEM_H_
#define SRC_GRAPH_NODE_ITEM_H_

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace raptor {

template<class NameType, class NodeDataType>
class NodeItem {
public:
    // friend std::shared_ptr<NodeItem<NameType, NodeDataType>> std::make_shared(
    //                                                             int64_t _internal_id,
    //                                                             const NameType& _name,
    //                                                             const std::shared_ptr<NodeDataType> _data);
    NodeItem(int64_t _internal_id, const NameType& _name, const std::shared_ptr<NodeDataType> _data)
                            :   internal_id_(_internal_id), name_(_name), data_(_data), is_removed_(false) { }

    static std::shared_ptr<NodeItem<NameType, NodeDataType>> createNodeItem(
                                                                int64_t _internal_id,
                                                                const NameType& _name,
                                                                const std::shared_ptr<NodeDataType> _data) {
        return std::make_shared<NodeItem<NameType, NodeDataType>>(_internal_id, _name, _data);
    }

    ~NodeItem() = default;

    std::string ToDOT() const {
        return "";
    }
    std::string ToJSON() const {
        std::ostringstream ss;
        ss << "[\"nitem\", {\"name\": " << name()
            << ", \"removed\": " << is_removed()
            << ", \"int_id\": " << internal_id()
            << ", \"data\": " << data()->ToJSON()
            << "}]";
        return ss.str();
    }

    const std::shared_ptr<NodeDataType> data() const {
        return data_;
    }

    int64_t internal_id() const {
        return internal_id_;
    }

    const NameType& name() const {
        return name_;
    }

    bool is_removed() const {
        return is_removed_;
    }
    void is_removed(bool val) {
        is_removed_ = val;
    }

private:
    NodeItem(const NodeItem&) = delete;
    NodeItem& operator=(const NodeItem&) = delete;

    int64_t internal_id_;
    NameType name_;
    const std::shared_ptr<NodeDataType> data_;
    bool is_removed_;
};

}

#endif
