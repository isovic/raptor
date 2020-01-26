/*
 * graph_mapping_result.cc
 *
 *  Created on: Dec 18, 2017
 *      Author: Ivan Sovic
 */

#include <containers/mapping_result/graph_mapping_result.h>
#include <containers/mapping_result/mapping_result_common.h>
#include <algorithm>
#include <tuple>
#include <limits>

namespace raptor {

std::shared_ptr<raptor::GraphMappingResult> createGraphMappingResult(int64_t qseq_id,
                                                       int64_t qseq_len,
                                                       std::string qseq_header,
                                                       mindex::IndexPtr index) {
    return std::shared_ptr<raptor::GraphMappingResult>(new GraphMappingResult(qseq_id, qseq_len, qseq_header, index));
}

GraphMappingResult::GraphMappingResult(int64_t _qseq_id,
                             int64_t _qseq_len,
                             std::string _qseq_header,
                             mindex::IndexPtr _index) :
                                      qseq_id_(_qseq_id),
                                      qseq_len_(_qseq_len),
                                      qseq_header_(_qseq_header),
                                      index_(_index),
                                      rv_(MapperReturnValueBase::NotRunYet),
                                      paths_(),
                                      all_similar_scores_(12, 0),
                                      fraction_query_covered_(0.0) {

}

// Interface implementation.
std::vector<std::shared_ptr<raptor::RegionBase>> GraphMappingResult::CollectRegions(bool one_hit_per_target, bool do_relabel_sec_supp) const {

    // // Sort paths by score.
    // auto sorted = paths_;
    // std::sort(sorted.begin(), sorted.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );

    std::vector<std::shared_ptr<raptor::RegionBase>> ret;
    std::vector<std::pair<std::shared_ptr<raptor::RegionBase>, int32_t>> secondary;

    // Used for filtering multiple hits to the same target.
    std::unordered_map<std::string, int8_t> query_target_pairs;

    // Collect the PRIMARY and SECONDARY alignments separately.
    for (size_t i = 0; i < paths_.size(); ++i) {
        const auto& path = paths_[i];
        if (path == nullptr) {
            continue;
        }
        // Merge neighboring nodes connected with implicit paths into a single node.
        auto merged_path = LocalPathTools::MergeImplicitEdges(path);
        if (merged_path == nullptr) {
            continue;
        }
        // Annotate all chunks as either Secondary or SecondarySupplementary.
        for (size_t aln_id = 0; aln_id < merged_path->nodes().size(); ++aln_id) {
            const auto& aln = merged_path->nodes()[aln_id]->data();
            if (one_hit_per_target) {
                std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
                if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
                    continue;
                }
                query_target_pairs[pair_name] = 1;
            }
            // Priority of 0 are the primary alignments.
            auto priority = aln->GetRegionPriority();
            if (priority == 0) {
                ret.emplace_back(aln);
            } else if (priority > 0) {
                secondary.emplace_back(std::make_pair(aln, path->score()));
            }
        }
    }

    // For now, output absolutely any alignment. Filtering should have happened earlier.
    for (const auto& vals: secondary) {
        ret.emplace_back(std::get<0>(vals));
    }

    return ret;
}

int64_t GraphMappingResult::QueryId() const {
    return qseq_id_;
}

int64_t GraphMappingResult::QueryLen() const {
    return qseq_len_;
}

std::string GraphMappingResult::QueryHeader() const {
    return qseq_header_;
}

const mindex::IndexPtr GraphMappingResult::Index() const {
    return index_;
}

MapperReturnValueBase GraphMappingResult::ReturnValue() const {
    return rv_;
}

const std::unordered_map<std::string, double>& GraphMappingResult::Timings() const {
    return timings_;
}



// Implementation-specific interfaces.
std::vector<int32_t> GraphMappingResult::CountSimilarMappings_(const std::vector<std::shared_ptr<raptor::LocalPath>>& paths, int32_t min_score_diff_margin) {

    std::vector<int32_t> score_counts(12, 0); // For the distribution between 0.00% and 10.0% in 1% steps.

    if (paths.size() == 0) {
        return score_counts;
    }

    auto sorted = paths;

    std::sort(sorted.begin(), sorted.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );

    double best_score = sorted[0]->score();
    double best_score_f = (double) std::abs(best_score);

    if (best_score == 0.0) {
        return score_counts;
    }

    for (size_t i = 0; i < sorted.size(); i++) {
        int32_t score = sorted[i]->score();
        int32_t score_diff = std::abs(best_score - score);
        double score_diff_f = (double) score_diff;
        double fraction_diff = score_diff_f / best_score_f;
        int32_t bucket = ((int32_t) std::floor(fraction_diff * 100)) + 1;

        // Allow a small numerical margin. Should be of the order of magnitude of a seed len.
        if (score_diff > 0 && score_diff <= min_score_diff_margin) {
            bucket = 1;
        }

        // Count all of the above ones in a single bucket.
        if (bucket > 10) {
            bucket = 11;
        }

        // Count only exact hits in bucket 0.
        // All hits between <0.00, 0.01> go into bucket 1.
        if (fraction_diff == 0.00) {
            score_counts[0] += 1;

        } else {
            score_counts[bucket] += 1;
        }
    }

    // Create the cumulative sum.
    for (size_t i = 1; i < score_counts.size(); i++) {
        score_counts[i] += score_counts[i - 1];
    }

    return score_counts;
}

double GraphMappingResult::CalcBestFractionQueryCovered_(const std::vector<std::shared_ptr<raptor::LocalPath>>& paths, int64_t qseq_len) {
    double ret = 0.0f;

    if (paths.size() == 0) {
        return ret;
    }

    int32_t best_i = -1;
    for (int32_t i = 0; i < (int32_t) paths.size(); i++) {
        if (best_i == -1 || paths[i]->score() > paths[best_i]->score()) {
            best_i = i;
        }
    }

    auto& best_path = paths[best_i];

    if (best_path->nodes().size() == 0) {
        return ret;
    }

    ret = ((double) (best_path->nodes().back()->data()->QueryEnd() - best_path->nodes().front()->data()->QueryStart())) / ((double) qseq_len);

    return ret;
}

std::vector<std::shared_ptr<raptor::LocalPath>> GraphMappingResult::GenerateFiltered(
                                                const std::vector<std::shared_ptr<raptor::LocalPath>>& paths,
                                                int32_t bestn,
                                                double max_fraction_diff,
                                                int32_t min_score_diff_margin,
                                                int32_t min_map_len
                                                ) {

    /*
     * This function keeps any path which has at least one region within the bestn
     * priority. Also, any path within the max_fraction_diff from the best scoring
     * path is kept even if it's beyond bestn priority.
     *
    */

    std::vector<std::shared_ptr<raptor::LocalPath>> ret;

    if (paths.empty()) {
        return ret;
    }

    // Find the best score.
    int64_t best_score = std::numeric_limits<int32_t>::lowest();
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        if (path == nullptr) { continue; }
        best_score = std::max(best_score, path->score());
    }
    // Keep any score above this fraction difference even if it's above bestn regions.
    int64_t min_score = best_score * (1.0 - max_fraction_diff);
    bool ret_contains_primary = false;

    // Filter the bestn alternative paths.
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        if (path == nullptr) { continue; }
        if (path->nodes().empty()) { continue; }

        // Keep path if any region has priority less than bestn.
        bool keep_path = false;
        bool contains_primary = false;
        for (size_t aln_id = 0; aln_id < path->nodes().size(); ++aln_id) {
            auto& aln = path->nodes()[aln_id]->data();
            auto priority = aln->GetRegionPriority();
            if (priority == 0) {
                contains_primary = true;
            }
            // If bestn == 1, then no secondary alignments should explicitly be taken.
            // If priority < 0, this region needs to be filtered because it's
            // secondary-to-primary score ratio is above the limit.
            if (bestn <= 0 || (bestn == 1 && priority == 0) || (bestn > 1 && priority >= 0 && priority < bestn)) {
                keep_path = true;
            }
        }

        // Keep any score above this fraction difference, but only if bestn != 1.
        // If bestn == 1, then only primary alignments should be kept.
        if (path->score() >= min_score && (bestn != 1 || (bestn == 1 && contains_primary))) {
            keep_path = true;
        }

        // Keep a path if it's score is within some absolute difference from the best.
        // This is intended to keep those paths which may differ in a k-mer or two.
        if ((best_score - path->score()) <= min_score_diff_margin && (bestn != 1 || (bestn == 1 && contains_primary))) {
            keep_path = true;
        }

        // Filter paths that are too short.
        int32_t qspan = path->nodes().front()->data()->QueryEnd() - path->nodes().front()->data()->QueryStart();
        if (qspan < min_map_len) {
            keep_path = false;
        }

        if (keep_path) {
            ret.emplace_back(path);
            if (contains_primary) {
                ret_contains_primary = true;
            }
        }
    }

    // If the primary alignment was not added to the "ret" vector, then do not
    // report any secondary or supplementary alignments either.
    if (ret_contains_primary == false) {
        ret = {};
    }

    return ret;

}

void GraphMappingResult::Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, bool just_sort) {
    paths_ = GenerateFiltered(paths_, bestn, max_fraction_diff, index_->params()->k * 3, min_map_len);
}

}
