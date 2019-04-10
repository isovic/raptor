/*
 * node_data_base.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_NODE_BASE_H_
#define SRC_GRAPH_NODE_BASE_H_

#include <string>

namespace raptor {

class NodeDataBase {
public:
    virtual ~NodeDataBase() = 0;
    virtual std::string ToDOT() const = 0;
    virtual std::string ToJSON() const = 0;
};

inline NodeDataBase::~NodeDataBase() {
}

}

#endif
