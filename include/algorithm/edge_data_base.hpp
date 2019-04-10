/*
 * edge_data_base.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_EDGE_BASE_H_
#define SRC_GRAPH_EDGE_BASE_H_

#include <string>

namespace raptor {

class EdgeDataBase {
public:
    virtual ~EdgeDataBase() = 0;
    virtual std::string ToDOT() const = 0;
    virtual std::string ToJSON() const = 0;
};

inline EdgeDataBase::~EdgeDataBase() {
}

}

#endif
