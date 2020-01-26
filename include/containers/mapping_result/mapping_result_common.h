/*
 * mapping_result_common.h
 *
 *  Created on: Aug 04, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_MAPPING_RESULT_COMMON_H_
#define SRC_CONTAINERS_MAPPING_RESULT_COMMON_H_

#include <memory>
#include <vector>
#include <unordered_map>
#include <containers/mapping_result/mapping_result_base.h>
#include <containers/region/region_base.h>
#include <containers/region/region_aligned.h>
#include <containers/region/region_mapped.h>
#include <types/typedefs_interval_tree.h>

namespace raptor {

void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& ret, double min_sec_to_prim_ratio);
void FindSupplementary(
                IntervalTreeInt64& query_trees, IntervalVectorInt64& query_intervals,
                std::unordered_map<int64_t, IntervalTreeInt64>& target_trees,
                std::unordered_map<int64_t, IntervalVectorInt64>& target_intervals,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio);
void FindSecondary(
                IntervalTreeInt64& query_trees, IntervalVectorInt64& query_intervals,
                std::unordered_map<int64_t, IntervalTreeInt64>& target_trees,
                std::unordered_map<int64_t, IntervalVectorInt64>& target_intervals,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio);

}

#endif
