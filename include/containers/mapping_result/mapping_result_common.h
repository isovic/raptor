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

void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& ret, double min_sec_to_prim_ratio, int32_t grace_dist_for_mapq_scaling);

/*
 * Helper function.
 *
 * For a given set of regions and a priority comparison function, create
 * interval trees in both query and target spaces.
 * The CompPriority function takes the RegionPriority() value and returns true
 * if the region should be accepted and added to the interval tree.
 * This function enables the client to build the trees from only the primary,
 * only the supplementary, or any other combination of regions.
 *
 * Returns results via function arguments:
 * - qi - Query intervals.
 * - qt - Query interval trees.
 * - ti - Target intervals. It's an unordered map, where the key is the TargetID() of a region.
 * - tt - Target interval trees. The key of the map is the TargetID() of the tree.
*/
void CreateRegionIntervalTrees(
        const std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
        std::function<bool(int32_t a)> CompPriority,
        IntervalVectorInt64& qi,
        IntervalTreeInt64& qt,
        std::unordered_map<int64_t, IntervalVectorInt64>& ti,
        std::unordered_map<int64_t, IntervalTreeInt64>& tt);

void FindSupplementary(
                IntervalTreeInt64& query_trees, IntervalVectorInt64& query_intervals,
                std::unordered_map<int64_t, IntervalTreeInt64>& target_trees,
                std::unordered_map<int64_t, IntervalVectorInt64>& target_intervals,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio);
void FindSecondary(
                const IntervalTreeInt64& qt_prim,
                const std::unordered_map<int64_t, IntervalTreeInt64>& tt_prim,
                const IntervalTreeInt64& qt_sec,
                const std::unordered_map<int64_t, IntervalTreeInt64>& tt_sec,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio);
void AssignMappingQuality(std::vector<std::shared_ptr<raptor::RegionBase>>& alns, int32_t grace_dist);

}

#endif
