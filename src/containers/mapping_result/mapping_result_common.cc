/*
 * mapping_result_common.cc
 *
 *  Created on: Aug 04, 2018
 *      Author: Ivan Sovic
 */

#include <containers/mapping_result/mapping_result_common.h>
#include <types/typedefs_interval_tree.h>
#include <unordered_map>
#include <log/log_tools.h>

namespace raptor {

void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& ret) {

    // Build the query and target interval trees for the primary alignments.
    std::vector<raptor::IntervalInt64> query_intervals;
    std::unordered_map<int64_t, std::vector<raptor::IntervalInt64>> encoded_target_intervals;    // The key is (target_id * 2 + (is_rev ? 1 : 0))
    for (int64_t i = 0; i < static_cast<int64_t>(ret.size()); ++i) {
        const auto& aln = ret[i];
        // Only expand the interval list if the alignment is a primary one.
        if (aln->GetRegionPriority() == 0) {
            query_intervals.emplace_back(raptor::IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
            // The target interval trees have to be a map, and take care of the reverse complement too.
            int64_t target_key = aln->TargetID() * 2 + (aln->TargetRev() ? 1 : 0);
            if (encoded_target_intervals.find(target_key) != encoded_target_intervals.end()) {
                encoded_target_intervals[target_key].emplace_back(raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
            } else {
                // encoded_target_intervals[target_key] = std::vector<raptor::IntervalInt64>(raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
                encoded_target_intervals[target_key] = {raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i)};
            }
        }
    }

    // Construct the interval trees from the intervals.
    raptor::IntervalTreeInt64 query_trees(query_intervals);
    std::unordered_map<int64_t, raptor::IntervalTreeInt64> encoded_target_trees;
    for (auto it = encoded_target_intervals.begin(); it != encoded_target_intervals.end(); ++it) {
        encoded_target_trees[it->first] = raptor::IntervalTreeInt64(encoded_target_intervals[it->first]);
    }

    // Linear pass through the non-primary alignments to relabel them.
    for (int64_t i = 0; i < static_cast<int64_t>(ret.size()); ++i) {
        // Not const, because the primary/supplementary tags will be updated below.
        auto& aln = ret[i];

        // #ifdef RAPTOR_TESTING_MODE
        DEBUG_RUN(true,
            LOG_ALL("[Supp. i = %ld] Checking: aln = '%s'\naln->GetRegionPriority() = %ld\naln->GetRegionIsSupplementary() = %s\n", i, aln->WriteAsCSV(',').c_str(), aln->GetRegionPriority(), (aln->GetRegionIsSupplementary() ? "true" : "false"))
        );
        // #endif

        // Skip any primary alignment.
        if (aln->GetRegionPriority() == 0) {
            DEBUG_RUN(true,
                LOG_ALL("    - Skipping, aln->GetRegionPriority() == 0.")
            );

            continue;
        }

        // Check if there are any overlaps with other regions in the query coordinates.
        std::vector<raptor::IntervalInt64> query_intervals;
        query_trees.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1, query_intervals);
        // If there are any overlapping regions, this is not a valid candidate.
        if (query_intervals.size() > 0) {
            DEBUG_RUN(true,
                LOG_ALL("    - Skipping, query_intervals.size() > 0 (query_intervals.size() = %lu)\n", query_intervals.size())
            );
            continue;
        }

        // Check if there are any overlaps with other regions in the target coordinates.
        int64_t target_key = aln->TargetID() * 2 + (aln->TargetRev() ? 1 : 0);
        // Find the target among the trees, if it exists.
        auto it_trees = encoded_target_trees.find(target_key);
        if (it_trees != encoded_target_trees.end()) {
            // Only lookup the intervals if the tree exists (of course).
            std::vector<raptor::IntervalInt64> target_intervals;
            it_trees->second.findOverlapping(aln->TargetStart(), aln->TargetEnd() - 1, target_intervals);
            // If there are any overlapping regions, this is not a valid candidate.
            if (target_intervals.size() > 0) {
                DEBUG_RUN(true,
                    LOG_ALL("    - Skipping, there is an overlapping interval in target coordinates. target_intervals.size() = %lu\n", target_intervals.size())
                );
                continue;
            }
        }

        DEBUG_RUN(true,
            LOG_ALL("    - Updating the supplementary info.")
        );

        // All tests passed. This is a valid supplementary alignment. Add it to the tree, and update it's priority.
        aln->SetRegionPriority(0);
        aln->SetRegionIsSupplementary(true);

        // Update the intervals and trees so that the next candidate region can be evaluated.
        {
            query_intervals.emplace_back(raptor::IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
            // The target interval trees have to be a map, and take care of the reverse complement too.
            int64_t target_key = aln->TargetID() * 2 + (aln->TargetRev() ? 1 : 0);
            if (encoded_target_intervals.find(target_key) != encoded_target_intervals.end()) {
                encoded_target_intervals[target_key].emplace_back(raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
            } else {
                // encoded_target_intervals[target_key] = std::vector<raptor::IntervalInt64>(raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
                encoded_target_intervals[target_key] = {raptor::IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i)};
            }
            query_trees = raptor::IntervalTreeInt64(query_intervals);
            encoded_target_trees[target_key] = raptor::IntervalTreeInt64(encoded_target_intervals[target_key]);
        }

    }
}

}  // namespace raptor
