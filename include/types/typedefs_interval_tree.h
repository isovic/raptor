/*
 * typedefs_interval_tree.h
 *
 *  Created on: Dec 18, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_TYPES_TYPEDEFS_INTERVAL_TREE_H_
#define SRC_TYPES_TYPEDEFS_INTERVAL_TREE_H_

#include <intervaltree/IntervalTree.h>

using IntervalTreeInt64 = IntervalTree<int64_t, size_t>;
using IntervalVectorInt64 = IntervalTreeInt64::interval_vector;
using IntervalInt64 = IntervalTreeInt64::interval;

#endif
