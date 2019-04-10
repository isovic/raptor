/*
 * typedefs_interval_tree.h
 *
 *  Created on: Dec 18, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_TYPES_TYPEDEFS_INTERVAL_TREE_H_
#define SRC_TYPES_TYPEDEFS_INTERVAL_TREE_H_

#include <intervaltree/IntervalTree.h>

namespace raptor {

typedef Interval<int64_t> IntervalInt64;
typedef IntervalTree<int64_t> IntervalTreeInt64;

}

#endif
