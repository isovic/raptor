/*
 * mapping_result_common.cc
 *
 *  Created on: Aug 04, 2019
 *      Author: Ivan Sovic
 */

#include <containers/mapping_result/mapping_result_common.h>
#include <unordered_map>
#include <log/log_tools.h>
#include <utility/mapq.hpp>
#include <set>

namespace raptor {

/*
 * Paths should ideally be sorted by score, so that the newly assigned region priority reflects
 * the score.
*/
void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& alns, double min_sec_to_prim_ratio) {
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
    std::unordered_map<int64_t, IntervalVectorInt64> target_intervals;    // The key is (target_id * 2 + (is_rev ? 1 : 0))
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        const auto& aln = alns[i];
        // Only expand the interval list if the alignment is a primary one.
        if (aln->GetRegionPriority() == 0) {
            query_intervals.emplace_back(IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
            target_intervals[aln->TargetID()].emplace_back(IntervalInt64(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1, i));
        }
    }
    // // Debug print before constructing the tree, because the IntervalTree moves the data.
    // LOG_ALL("Query intervals:\n");
    // for (size_t i = 0; i < query_intervals.size(); ++i) {
    //     LOG_NOHEADER("[query interval i = %ld] (%ld, %ld, %ld)\n", i, query_intervals[i].start, query_intervals[i].stop, query_intervals[i].value);
    // }
    // Construct the query interval trees from the intervals.
    IntervalTreeInt64 query_trees;
    {
        auto query_intervals_copy = query_intervals;
        query_trees = IntervalTreeInt64(std::move(query_intervals_copy));
    }
    // Construct the target interval trees from the intervals.
    std::unordered_map<int64_t, IntervalTreeInt64> target_trees;
    for (auto it = target_intervals.begin(); it != target_intervals.end(); ++it) {
        auto target_intervals_copy = target_intervals[it->first];
        target_trees[it->first] = IntervalTreeInt64(std::move(target_intervals_copy));
    }

    FindSupplementary(query_trees, query_intervals, target_trees, target_intervals, alns, min_sec_to_prim_ratio);
    FindSecondary(query_trees, query_intervals, target_trees, target_intervals, alns, min_sec_to_prim_ratio);
}

void FindSupplementary(
                IntervalTreeInt64& query_trees, IntervalVectorInt64& query_intervals,
                std::unordered_map<int64_t, IntervalTreeInt64>& target_trees,
                std::unordered_map<int64_t, IntervalVectorInt64>& target_intervals,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio) {

    // Linear pass through the non-primary alignments to relabel them.
    // Here, we find only the SUPPLEMENTARY alignments. Secondary are ignored until
    // all suppl. regions are found (otherwise we might miss overlaps).
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        // Not const, because the primary/supplementary tags will be updated below.
        auto& aln = alns[i];
        DEBUG_RUN(true,
            LOG_ALL("[FindSupplementary i = %ld] Checking: aln = '%s'\n(before) aln->GetRegionPriority() = %ld\n(before) aln->GetRegionIsSupplementary() = %s\n", i, aln->WriteAsCSV(',').c_str(), aln->GetRegionPriority(), (aln->GetRegionIsSupplementary() ? "true" : "false"))
        );
        // Skip any primary alignment.
        if (aln->GetRegionPriority() == 0) {
            DEBUG_RUN(true, LOG_ALL("    - Skipping, aln->GetRegionPriority() == 0.\n"));
            continue;
        }

        // SECONDARY ALIGNMENTS based on the query coordinate overlaps.
        DEBUG_RUN(true, LOG_ALL("    - Checking query interval: (%ld, %ld)\n", aln->QueryStart(), aln->QueryEnd() - 1));
        auto found_query_intervals = query_trees.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        if (found_query_intervals.size()) {
            DEBUG_RUN(true, LOG_ALL("    - Skipping, overlaps in query.\n"));
            continue;
        }
        // SECONDARY ALIGNMENTS based on the target coordinate overlaps.
        auto it_trees = target_trees.find(aln->TargetID());
        if (it_trees != target_trees.end()) {
            auto found_target_intervals = it_trees->second.findOverlapping(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1);
            if (found_target_intervals.size()) {
                DEBUG_RUN(true, LOG_ALL("    - Skipping, overlaps in target.\n"));
                continue;
            }
        }

        DEBUG_RUN(true, LOG_ALL("    - Supplementary region found.\n"));
        aln->SetRegionPriority(0);
        aln->SetRegionIsSupplementary(true);

        // Update the query interval trees.
        query_intervals.emplace_back(IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
        auto query_intervals_copy = query_intervals;
        query_trees = IntervalTreeInt64(std::move(query_intervals_copy));
        // Update the target interval trees.
        target_intervals[aln->TargetID()].emplace_back(IntervalInt64(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1, i));
        auto target_intervals_copy = target_intervals[aln->TargetID()];
        target_trees[aln->TargetID()] = IntervalTreeInt64(std::move(target_intervals_copy));
    }
}

void FindSecondary(
                IntervalTreeInt64& query_trees, IntervalVectorInt64& query_intervals,
                std::unordered_map<int64_t, IntervalTreeInt64>& target_trees,
                std::unordered_map<int64_t, IntervalVectorInt64>& target_intervals,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio) {

    /*
     * Find all overlaps with the primary mappings from the IntervalTrees.
     * If there are no overlaps for a mapping, relabel it as a supplementary.
     * Otherwise, mark it as a secondary alignment and re-count the region's priority.
    */
    // Mark any non-primary alignment with a separate priority.
    int64_t next_region_priority = 1;

    // Reset all region counts for all alignments first.
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        aln->SetAltRegionCount(1);
    }

    // Linear pass through the non-primary alignments to relabel them.
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        // Not const, because the primary/supplementary tags will be updated below.
        auto& aln = alns[i];
        DEBUG_RUN(true,
            LOG_ALL("[Supp. i = %ld] Checking: aln = '%s'\n(before) aln->GetRegionPriority() = %ld\n(before) aln->GetRegionIsSupplementary() = %s\n", i, aln->WriteAsCSV(',').c_str(), aln->GetRegionPriority(), (aln->GetRegionIsSupplementary() ? "true" : "false"))
        );
        // Skip any primary alignment.
        if (aln->GetRegionPriority() == 0) {
            DEBUG_RUN(true, LOG_ALL("    - Skipping, aln->GetRegionPriority() == 0.\n"));
            continue;
        }

        std::set<int64_t> ovl_ids;
        raptor::RegionType type(raptor::RegionType::Undefined);

        // SECONDARY ALIGNMENTS based on the query coordinate overlaps.
        DEBUG_RUN(true, LOG_ALL("    - Checking query interval: (%ld, %ld)\n", aln->QueryStart(), aln->QueryEnd() - 1));
        auto found_query_intervals = query_trees.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        for (auto& interval: found_query_intervals) {
            auto& primary_aln = alns[interval.value];
            auto ratio = static_cast<double>(aln->Score()) / static_cast<double>(primary_aln->Score());
            ratio = std::max(0.0, ratio);
            // std::cerr << "(query interval) Current aln: " << aln->WriteAsCSV(',') << "\n";
            // std::cerr << "(query interval) Primary aln: " << primary_aln->WriteAsCSV(',') << "\n";
            // std::cerr << "(query interval)ratio = " << ratio << "\n";
            if (ratio < min_sec_to_prim_ratio) {
                continue;
            }
            ovl_ids.emplace(interval.value);
        }

        // SECONDARY ALIGNMENTS based on the target coordinate overlaps.
        auto it_trees = target_trees.find(aln->TargetID());
        if (it_trees != target_trees.end()) {
            auto found_target_intervals = it_trees->second.findOverlapping(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1);
            for (auto& interval: found_target_intervals) {
                auto& primary_aln = alns[interval.value];
                auto ratio = static_cast<double>(aln->Score()) / static_cast<double>(primary_aln->Score());
                ratio = std::max(0.0, ratio);
                // std::cerr << "(target interval) Current aln: " << aln->WriteAsCSV(',') << "\n";
                // std::cerr << "(target interval) Primary aln: " << primary_aln->WriteAsCSV(',') << "\n";
                // std::cerr << "(target interval)ratio = " << ratio << "\n";
                if (ratio < min_sec_to_prim_ratio) {
                    continue;
                }
                ovl_ids.emplace(interval.value);
            }
        }

        // Resolve the type of this region.
        type = (ovl_ids.size()) ? raptor::RegionType::Secondary : raptor::RegionType::Undefined;

        if (type == raptor::RegionType::Secondary) {
            DEBUG_RUN(true, LOG_ALL("    - Secondary. Found %ld overlaps with other regions.\n", ovl_ids.size()));
            aln->SetRegionPriority(next_region_priority);
            aln->SetRegionIsSupplementary(false);
            ++next_region_priority;

            // Count the alternative hits covering the same region, and increase the count of the Primary.
            for (const auto ovl_id: ovl_ids) {
                auto& primary_aln = alns[ovl_id];
                int32_t count = primary_aln->GetAltRegionCount();
                primary_aln->SetAltRegionCount(count + 1);
            }
            // The current secondary alignment's alt count is the same as the number ovl_ids.
            aln->SetAltRegionCount(static_cast<int32_t>(ovl_ids.size() + 1));

        } else {
            // These should be filtered out eventually. This means that the score overlap is below a specified threshold.
            DEBUG_RUN(true, LOG_ALL("    - Region does not satisfy conditions.\n"));
            aln->SetRegionPriority(-1);
            aln->SetRegionIsSupplementary(false);
        }
    }

    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        int32_t mapq = CalcMapq(aln->GetAltRegionCount());
        aln->SetMappingQuality(mapq);
    }

    // std::cerr << "---------------\n";
}

}  // namespace raptor
