/*
 * mapping_result_common.cc
 *
 *  Created on: Aug 04, 2019
 *      Author: Ivan Sovic
 */

#include <containers/mapping_result/mapping_result_common.h>
#include <types/typedefs_interval_tree.h>
#include <unordered_map>
#include <log/log_tools.h>

namespace raptor {

/*
 * Paths should ideally be sorted by score, so that the newly assigned region priority reflects
 * the score.
*/
void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& sorted_paths) {
    /*
        General algorithm:
        1. Sort the paths by score.
        2. Find the best scoring path. This one will be marked as the primary.
        3. Collect all alignments.
        4. Mark the primary ones.
        5. For any unmarked, check if they do not overlap the primary path. These will be the supplementary ones.
        6. Mark the rest as secondary.
    */

    /*
     * First, build the query and target interval trees for primary alignments (priority == 0).
    */
    IntervalVectorInt64 query_intervals;
    std::unordered_map<int64_t, IntervalVectorInt64> encoded_target_intervals;    // The key is (target_id * 2 + (is_rev ? 1 : 0))
    for (int64_t i = 0; i < static_cast<int64_t>(sorted_paths.size()); ++i) {
        const auto& aln = sorted_paths[i];
        // Only expand the interval list if the alignment is a primary one.
        if (aln->GetRegionPriority() == 0) {
            query_intervals.emplace_back(IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
            // Encode the target key to contain the strand info as the LSB bit.
            // The target interval trees have to be a map, and take care of the reverse complement too.
            int64_t target_key = aln->TargetID();
            int64_t tstart_fwd = aln->TargetFwdStart();
            int64_t tend_fwd = aln->TargetFwdEnd();
            if (encoded_target_intervals.find(target_key) != encoded_target_intervals.end()) {
                encoded_target_intervals[target_key].emplace_back(IntervalInt64(tstart_fwd, tend_fwd - 1, i));
            } else {
                // encoded_target_intervals[target_key] = std::vector<IntervalInt64>(IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
                encoded_target_intervals[target_key] = {IntervalInt64(tstart_fwd, tend_fwd - 1, i)};
            }
        }
    }
    // Debug print before constructing the tree, because the IntervalTree moves the data.
    LOG_ALL("Query intervals:\n");
    for (size_t i = 0; i < query_intervals.size(); ++i) {
        LOG_NOHEADER("[query interval i = %ld] (%ld, %ld, %ld)\n", i, query_intervals[i].start, query_intervals[i].stop, query_intervals[i].value);
    }
    // Construct the query interval trees from the intervals.
    IntervalTreeInt64 query_trees;
    {
        auto query_intervals_copy = query_intervals;
        query_trees = IntervalTreeInt64(std::move(query_intervals_copy));
    }
    // Construct the target interval trees from the intervals.
    std::unordered_map<int64_t, IntervalTreeInt64> encoded_target_trees;
    for (auto it = encoded_target_intervals.begin(); it != encoded_target_intervals.end(); ++it) {
        auto target_intervals_copy = encoded_target_intervals[it->first];
        encoded_target_trees[it->first] = IntervalTreeInt64(std::move(target_intervals_copy));
    }

    /*
     * Find all overlaps with the primary mappings from the IntervalTrees.
     * If there are no overlaps for a mapping, relabel it as a supplementary.
     * Otherwise, mark it as a secondary alignment and re-count the region's priority.
    */
    // Mark any non-primary alignment with a separate priority.
    int64_t next_region_priority = 1;
    // Linear pass through the non-primary alignments to relabel them.
    for (int64_t i = 0; i < static_cast<int64_t>(sorted_paths.size()); ++i) {
        // Not const, because the primary/supplementary tags will be updated below.
        auto& aln = sorted_paths[i];
        DEBUG_RUN(true,
            LOG_ALL("[Supp. i = %ld] Checking: aln = '%s'\n(before) aln->GetRegionPriority() = %ld\n(before) aln->GetRegionIsSupplementary() = %s\n", i, aln->WriteAsCSV(',').c_str(), aln->GetRegionPriority(), (aln->GetRegionIsSupplementary() ? "true" : "false"))
        );
        // Skip any primary alignment.
        if (aln->GetRegionPriority() == 0) {
            DEBUG_RUN(true,
                LOG_ALL("    - Skipping, aln->GetRegionPriority() == 0.\n")
            );

            continue;
        }

        // SECONDARY ALIGNMENTS based on the query coordinate overlaps.
        DEBUG_RUN(true,
            LOG_ALL("    - Checking query interval: (%ld, %ld)\n", aln->QueryStart(), aln->QueryEnd() - 1)
        );
        auto found_query_intervals = query_trees.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        // If there are any overlapping regions, this is not a valid candidate.
        if (found_query_intervals.size() > 0) {
            DEBUG_RUN(true,
                LOG_ALL("    - Skipping, found_query_intervals.size() > 0 (found_query_intervals.size() = %lu)\n", found_query_intervals.size())
            );
            aln->SetRegionPriority(next_region_priority);
            aln->SetRegionIsSupplementary(false);
            ++next_region_priority;
            continue;
        }

        // SECONDARY ALIGNMENTS based on the target coordinate overlaps.
        // Check if there are any overlaps with other regions in the target coordinates.
        // We need to check the fwd strand, because query could have a missing adapter and align
        // over the same spot, but on a different strand.
        int64_t target_key = aln->TargetID();
        int64_t tstart_fwd = aln->TargetFwdStart();
        int64_t tend_fwd = aln->TargetFwdEnd();
        // Find the target among the trees, if it exists.
        auto it_trees = encoded_target_trees.find(target_key);
        if (it_trees != encoded_target_trees.end()) {
            DEBUG_RUN(true,
                LOG_ALL("    - Found the target in target trees. Looking for overlaps.\n")
            );
            // Only lookup the intervals if the tree exists (of course).
            // std::vector<IntervalInt64> target_intervals;
            auto found_target_intervals = it_trees->second.findOverlapping(tstart_fwd, tend_fwd - 1);
            DEBUG_RUN(true,
                LOG_ALL("    - There are %lu overlapping target intervals.\n", found_target_intervals.size())
            );
            // If there are any overlapping regions, this is not a valid candidate.
            if (found_target_intervals.size() > 0) {
                DEBUG_RUN(true,
                    LOG_ALL("    - Skipping, there is an overlapping interval in target coordinates. found_target_intervals.size() = %lu\n", found_target_intervals.size())
                );
                aln->SetRegionPriority(next_region_priority);
                aln->SetRegionIsSupplementary(false);
                ++next_region_priority;
                continue;
            }
        }

        // SUPPLEMENTARY ALIGNMENTS. If the region was not marked as a secondary, it has no overlaps
        // with the primary regino. Mark it as a supplementary.
        // Add it to the tree, and update it's priority.
        aln->SetRegionPriority(0);
        aln->SetRegionIsSupplementary(true);

        DEBUG_RUN(true,
            LOG_ALL("    - Updating the supplementary info.\n")
            LOG_ALL("\n(after) aln->GetRegionPriority() = %ld\n(after) aln->GetRegionIsSupplementary() = %s\n", aln->GetRegionPriority(), (aln->GetRegionIsSupplementary() ? "true" : "false"))
        );

        // Update the intervals and trees so that the next candidate region can be evaluated.
        {
            query_intervals.emplace_back(IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
            // The target interval trees have to be a map, and take care of the reverse complement too.
            int64_t target_key = aln->TargetID() * 2 + (aln->TargetRev() ? 1 : 0);
            if (encoded_target_intervals.find(target_key) != encoded_target_intervals.end()) {
                encoded_target_intervals[target_key].emplace_back(IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
            } else {
                // encoded_target_intervals[target_key] = std::vector<IntervalInt64>(IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i));
                encoded_target_intervals[target_key] = {IntervalInt64(aln->TargetStart(), aln->TargetEnd() - 1, i)};
            }

            // Update the query interval tree.
            {
                auto query_intervals_copy = query_intervals;
                query_trees = IntervalTreeInt64(std::move(query_intervals_copy));

                auto target_intervals_copy = encoded_target_intervals[target_key];
                encoded_target_trees[target_key] = IntervalTreeInt64(std::move(target_intervals_copy));
            }
        }

    }
}

}  // namespace raptor
