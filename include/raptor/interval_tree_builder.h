/*
 * interval_tree_builder.h
 *
 *  Created on: Jan 09, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_INTERVAL_TREE_BUILDER_H_
#define SRC_RAPTOR_INTERVAL_TREE_BUILDER_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <unordered_map>
#include <containers/region/region_mapped.h>
#include <types/typedefs.h>
#include <types/typedefs_interval_tree.h>

namespace raptor {

std::unordered_map<int32_t, IntervalTreeInt64> BuildAnchorIntervalTrees(
    const std::vector<raptor::AnchorPtr>& anchors);

std::unordered_map<int32_t, IntervalTreeInt64> BuildTargetAnchorIntervalTrees(
    const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors);

}  // namespace raptor

#endif
