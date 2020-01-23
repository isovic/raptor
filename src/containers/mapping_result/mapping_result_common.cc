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
#include <functional>

namespace raptor {

/*
 * Paths should ideally be sorted by score, so that the newly assigned region priority reflects
 * the score.
*/
void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& alns, double min_sec_to_prim_ratio, int32_t grace_dist_for_mapq_scaling) {
    /*
        General algorithm:
        1. Collect all regions and build 2 interval trees: one for query coordintes, and one for target coordinates.
        2. Collect only the primary regions and build 2 interval trees: one for query coordintes, and one for target coordinates.
        3. Use the trees of primary regions to find supplementary alignments. Supp. alignments should not overlap
            any of the primary alignments. Once a new suppl. alignment is found, expand the interval tree.
        4. Anything left over is a secondary alignment.
        5. Count the occurrences of all secondary alignments by looking up their hits in the query and
            target interval trees of ALL regions.
        6. Compute the mapping quality for all alignments.
    */

    // Compute the interval trees for the Primary regions.
    IntervalVectorInt64 qi_prim;
    std::unordered_map<int64_t, IntervalVectorInt64> ti_prim;
    IntervalTreeInt64 qt_prim;
    std::unordered_map<int64_t, IntervalTreeInt64> tt_prim;
    // The comparison function compares the region priority.
    // If it's 0, it's a primary region.
    CreateRegionIntervalTrees(alns, [](int32_t a){ return a == 0; },
                                qi_prim, qt_prim, ti_prim, tt_prim);

    // Label the supplementary, and mark them as primary too.
    FindSupplementary(qt_prim, qi_prim, tt_prim, ti_prim, alns, min_sec_to_prim_ratio);



    // Compute the interval trees for the Secondary regions.
    // This requires a second pass, because the primary regions were just
    // updated in FindSupplementary.
    IntervalVectorInt64 qi_sec;
    std::unordered_map<int64_t, IntervalVectorInt64> ti_sec;
    IntervalTreeInt64 qt_sec;
    std::unordered_map<int64_t, IntervalTreeInt64> tt_sec;
    // The comparison function compares the region priority.
    // Values > 0 are the secondary regions.
    CreateRegionIntervalTrees(alns, [](int32_t a){ return a > 0; },
                                qi_sec, qt_sec, ti_sec, tt_sec);

    // Label the secondary and count the occurrences.
    FindSecondary(qt_prim, tt_prim, qt_sec, tt_sec, alns, min_sec_to_prim_ratio);

    // Compute and assign the mapping quality based on the occurrence counts.
    AssignMappingQuality(alns, grace_dist_for_mapq_scaling);
}

void CreateRegionIntervalTrees(
        const std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
        std::function<bool(int32_t a)> CompPriority,
        IntervalVectorInt64& qi,
        IntervalTreeInt64& qt,
        std::unordered_map<int64_t, IntervalVectorInt64>& ti,
        std::unordered_map<int64_t, IntervalTreeInt64>& tt) {

    qi.clear();
    ti.clear();
    tt.clear();

    // IntervalVectorInt64 query_intervals_secondary;
    // std::unordered_map<int64_t, IntervalVectorInt64> target_intervals_secondary;
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        const auto& aln = alns[i];
        // Skip regions that are not of interest.
        if (CompPriority(aln->GetRegionPriority()) == false) {
            continue;
        }
        // Accumulate the intervals.
        qi.emplace_back(IntervalInt64(aln->QueryStart(), aln->QueryEnd() - 1, i));
        ti[aln->TargetID()].emplace_back(IntervalInt64(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1, i));
    }

    // Construct the query interval trees from the intervals.
    auto qi_copy = qi;
    qt = IntervalTreeInt64(std::move(qi_copy));

    // Construct the target interval trees from the intervals.
    for (auto it = ti.begin(); it != ti.end(); ++it) {
        auto ti_copy = ti[it->first];
        tt[it->first] = IntervalTreeInt64(std::move(ti_copy));
    }
}

void FindSupplementary(
                IntervalTreeInt64& query_trees,
                IntervalVectorInt64& query_intervals,
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
                const IntervalTreeInt64& qt_prim,
                const std::unordered_map<int64_t, IntervalTreeInt64>& tt_prim,
                const IntervalTreeInt64& qt_sec,
                const std::unordered_map<int64_t, IntervalTreeInt64>& tt_sec,
                std::vector<std::shared_ptr<raptor::RegionBase>>& alns,
                double min_sec_to_prim_ratio) {

    /*
     * Find all overlaps with the primary mappings from the IntervalTrees.
     * If there are no overlaps for a mapping, relabel it as a supplementary.
     * Otherwise, mark it as a secondary alignment and re-count the region's priority.
    */
    // Reset all region counts for all alignments first.
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        aln->SetAltRegionCount(1);
    }

    // First pass. For every secondary alignment, verify that it satisfies the
    // min_sec_to_prim_ratio criteria. If it does satisfy it, we keep the
    // secondary alignment, otherwise we need to filter it out.
    // Pair: <ID, score>.
    std::vector<std::pair<int32_t, int32_t>> valid_secondary;
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        // Skip any primary alignment, and all the filtered alignments.
        if (aln->GetRegionPriority() <= 0) {
            continue;
        }
        bool is_valid = false;

        auto found_qi_prim = qt_prim.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        for (auto& interval: found_qi_prim) {
            auto& other_aln = alns[interval.value];
            auto ratio = std::max(0.0, static_cast<double>(aln->Score()) / static_cast<double>(other_aln->Score()));
            if (ratio < min_sec_to_prim_ratio) {
                continue;
            }
            is_valid = true;
            break;
        }
        // The condition is not needed here, but it could save some compute cycles if we already
        // determined that this is valid.
        if (is_valid == false) {
            auto it_trees_prim = tt_prim.find(aln->TargetID());
            if (it_trees_prim != tt_prim.end()) {
                auto found_ti = it_trees_prim->second.findOverlapping(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1);
                for (auto& interval: found_ti) {
                    auto& other_aln = alns[interval.value];
                    auto ratio = std::max(0.0, static_cast<double>(aln->Score()) / static_cast<double>(other_aln->Score()));
                    if (ratio < min_sec_to_prim_ratio) {
                        continue;
                    }
                    is_valid = true;
                    break;
                }
            }
        }
        // If no valid primary regions, deactivate this region
        // (mark it's priority to 0).
        if (is_valid == false) {
            aln->SetRegionPriority(-1);
        } else {
            valid_secondary.emplace_back(std::make_pair(i, aln->Score()));
        }
    }

    // Sort by score in reverse, so that secondary regions will be marked in ascending priority by score.
    std::sort(valid_secondary.begin(), valid_secondary.end(), [](const auto& a, const auto& b) { return std::get<1>(a) > std::get<1>(b); });

    // Mark any non-primary alignment with a separate priority.
    int64_t next_region_priority = 1;

    // Second pass: go through all valid secondary regions, and count the overlaps.
    for (const auto& vals: valid_secondary) {
        int32_t aln_id = std::get<0>(vals);
        int32_t score = std::get<1>(vals);
        auto& aln = alns[aln_id];

        std::set<int64_t> prim_hits;
        std::set<int64_t> sec_hits;

        // Check the query intervals.
        auto found_qi_prim = qt_prim.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        for (auto& interval: found_qi_prim) {
            if (interval.value == aln_id || alns[interval.value]->GetRegionPriority() < 0) {
                continue;
            }
            prim_hits.emplace(interval.value);
        }
        auto found_qi_sec = qt_sec.findOverlapping(aln->QueryStart(), aln->QueryEnd() - 1);
        for (auto& interval: found_qi_sec) {
            // Skip self hits, and hits to filtered regions.
            if (interval.value == aln_id || alns[interval.value]->GetRegionPriority() < 0) {
                continue;
            }
            sec_hits.emplace(interval.value);
        }

        // Check the target intervals.
        auto it_trees_prim = tt_prim.find(aln->TargetID());
        if (it_trees_prim != tt_prim.end()) {
            auto found_ti = it_trees_prim->second.findOverlapping(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1);
            for (auto& interval: found_ti) {
                if (interval.value == aln_id || alns[interval.value]->GetRegionPriority() < 0) {
                    continue;
                }
                prim_hits.emplace(interval.value);
            }
        }
        auto it_trees_sec = tt_sec.find(aln->TargetID());
        if (it_trees_sec != tt_sec.end()) {
            auto found_ti = it_trees_sec->second.findOverlapping(aln->TargetFwdStart(), aln->TargetFwdEnd() - 1);
            for (auto& interval: found_ti) {
                // Skip self hits, and hits to filtered regions.
                if (interval.value == aln_id || alns[interval.value]->GetRegionPriority() < 0) {
                    continue;
                }
                sec_hits.emplace(interval.value);
            }
        }

        aln->SetRegionPriority(next_region_priority);
        aln->SetRegionIsSupplementary(false);
        ++next_region_priority;

        // Count the alternative hits covering the same region, and increase the count of the Primary.
        for (const auto ovl_id: prim_hits) {
            int32_t count = alns[ovl_id]->GetAltRegionCount();
            alns[ovl_id]->SetAltRegionCount(count + 1);
        }

        std::set<int64_t> all_hits = prim_hits;
        all_hits.insert(sec_hits.begin(), sec_hits.end());
        aln->SetAltRegionCount(static_cast<int32_t>(all_hits.size() + 1));
    }
}

void AssignMappingQuality(std::vector<std::shared_ptr<raptor::RegionBase>>& alns, int32_t grace_dist) {
    if (alns.empty()) {
        return;
    }

    auto sigmoid = [](double x){ double L = 1.0, k = 12.0, x0 = 0.5; return L / (1.0 + exp(-k * (x - x0))); };

    // Compute a scaling factor for  the largest mapping quality range (from 3 to 60).
    // This will allow more resolution for some of the sketchy alignments.
    // In this particular case, the scaling is based on the span of the overlaps.
    int32_t qlen = alns.front()->QueryLen();
    int32_t qspan = 0;
    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        if (aln->GetRegionPriority() == 0) {
            qspan += aln->QuerySpan();
        }
    }
    double factor = static_cast<double>(std::min(qlen, qspan + grace_dist)) / static_cast<double>(qlen);
    factor = sigmoid(factor);

    for (int64_t i = 0; i < static_cast<int64_t>(alns.size()); ++i) {
        auto& aln = alns[i];
        int32_t mapq_clean = CalcMapq(aln->GetAltRegionCount());
        int32_t next_mapq = CalcMapq(aln->GetAltRegionCount() + 1);

        int32_t mapq = next_mapq + static_cast<int32_t>(std::ceil((mapq_clean - next_mapq) * factor));

        aln->SetMappingQuality(mapq);
    }
}

}  // namespace raptor
