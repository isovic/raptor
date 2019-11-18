/*
 * interval_tree_builder.cc
 *
 *  Created on: Jan 09, 2018
 *      Author: Ivan Sovic
 */

#include <raptor/interval_tree_builder.h>

namespace raptor {

std::unordered_map<int32_t, IntervalTreeInt64> BuildAnchorIntervalTrees(
    const std::vector<raptor::AnchorPtr>& anchors) {
    // This maps the target ID to the tree of all chain intervals.
    std::unordered_map<int32_t, IntervalTreeInt64> trees;

    // We need to poll intervals first, so that the trees can be constructed.
    std::unordered_map<int32_t, IntervalVectorInt64> intervals;

    for (size_t i = 0; i < anchors.size(); i++) {
        auto& anchor = anchors[i];
        int32_t t_id = anchor->TargetID();
        int32_t t_start = anchor->TargetStart();
        int32_t t_end = anchor->TargetEnd();
        bool t_rev = anchor->TargetRev();

        auto it = intervals.find(t_id);
        if (it == intervals.end()) {
            intervals[t_id] = IntervalVectorInt64{};
        }

        intervals[t_id].emplace_back(IntervalInt64(t_start, t_end, i));
    }

    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
        trees[it->first] = IntervalTreeInt64(std::move(it->second));
    }

    return trees;
}

std::unordered_map<int32_t, IntervalTreeInt64> BuildTargetAnchorIntervalTrees(
    const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors) {
    // This maps the target ID to the tree of all chain intervals.
    std::unordered_map<int32_t, IntervalTreeInt64> trees;

    // We need to poll intervals first, so that the trees can be constructed.
    std::unordered_map<int32_t, IntervalVectorInt64> intervals;

    for (size_t i = 0; i < target_anchors.size(); i++) {
        // Uses `i` instead of `t_id` so that we know where exactly the
        // anchors came from.
        // int32_t t_id = target_anchors[i]->env()->t_id;
        bool t_rev = target_anchors[i]->env()->t_rev;

        auto it = intervals.find(i);
        if (it == intervals.end()) {
            intervals[i] = IntervalVectorInt64{};
            it = intervals.find(i);
        }

        for (size_t j = 0; j < target_anchors[i]->hits().size(); j++) {
            auto& hit = target_anchors[i]->hits()[j];
            int32_t t_start = hit->TargetStart();
            int32_t t_end = hit->TargetEnd();
            it->second.emplace_back(IntervalInt64(t_start, t_end, j));
        }
    }

    for (auto it = intervals.begin(); it != intervals.end(); ++it) {
        trees[it->first] = IntervalTreeInt64(std::move(it->second));
    }

    return trees;
}

}  // namespace raptor
