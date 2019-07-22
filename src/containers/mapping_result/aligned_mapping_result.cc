#include <containers/mapping_result/aligned_mapping_result.h>
#include <aligner/aligner_util.hpp>
#include <utility/revcmp.hpp>
#include <math.h>

namespace raptor {

std::shared_ptr<raptor::AlignedMappingResult> createRaptorAlignmentResult(const mindex::SequencePtr& qseq,
                                                                       mindex::IndexPtr index) {
    return std::shared_ptr<raptor::AlignedMappingResult>(new raptor::AlignedMappingResult(qseq, index));
}

AlignedMappingResult::AlignedMappingResult(const mindex::SequencePtr& qseq, mindex::IndexPtr _index)
    : qseq_(qseq),
      index_(_index),
      rv_(raptor::MapperReturnValueBase::NotRunYet),
      path_alignments_(),
      all_similar_scores_(12, 0),
      fraction_query_covered_(0.0) {}

// Interface implementation.
std::vector<std::shared_ptr<raptor::RegionBase>> AlignedMappingResult::CollectRegions(bool one_hit_per_target) const {
    std::vector<std::shared_ptr<raptor::RegionBase>> ret;

    // 1. Sort the paths by score.
    // 2. Find the best scoring path. This one will be marked as the primary.
    // 3. Collect all alignments.
    // 4. Mark the primary ones.
    // 5. For any unmarked, check if they do not overlap the primary path. These will be the supplementary ones.
    // 6. Mark the rest as secondary.
    std::vector<std::pair<int64_t, size_t>> path_scores;
    for (size_t path_id = 0; path_id < path_alignments().size(); ++path_id) {
        const std::shared_ptr<raptor::PathAlignment>& path_aln = path_alignments()[path_id];
        path_scores.emplace_back(std::make_pair(path_aln->path_score(), path_id));
    }
    std::sort(path_scores.begin(), path_scores.end());



    // Used for filtering multiple hits to the same target.
    std::unordered_map<std::string, int32_t> query_target_pairs;

    // Go through each path.
    for (size_t i = 0; i < path_alignments().size(); i++) {
        const std::shared_ptr<raptor::PathAlignment>& path_aln = path_alignments()[i];

        for (size_t j = 0; j < path_aln->alns().size(); j++) {
            auto& aln = path_aln->alns()[j];

            if (one_hit_per_target) {
                std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
                if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
                    continue;
                }
                query_target_pairs[pair_name] = 1;
            }
            ret.emplace_back(aln);
        }
    }

//    std::cerr << "[AlignedMappingResult::CollectRegions]:\n";
//    for (size_t i = 0; i < ret.size(); ++i) {
//        std::cerr << "[" << i << "] " << ret[i]->WriteAsCSV('\t') << "\n";
//    }

    return ret;
}

// // Interface implementation.
// std::vector<std::shared_ptr<raptor::RegionBase>> AlignedMappingResult::CollectRegions(bool one_hit_per_target) const {
//     std::vector<std::shared_ptr<raptor::RegionBase>> ret;

//     // Used for filtering multiple hits to the same target.
//     std::unordered_map<std::string, int32_t> query_target_pairs;

//     // Go through each path.
//     for (size_t i = 0; i < path_alignments().size(); i++) {
//         const std::shared_ptr<raptor::PathAlignment>& path_aln = path_alignments()[i];

//         for (size_t j = 0; j < path_aln->alns().size(); j++) {
//             auto& aln = path_aln->alns()[j];

//             if (one_hit_per_target) {
//                 std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
//                 if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
//                     continue;
//                 }
//                 query_target_pairs[pair_name] = 1;
//             }
//             ret.emplace_back(aln);
//         }
//     }

// //    std::cerr << "[AlignedMappingResult::CollectRegions]:\n";
// //    for (size_t i = 0; i < ret.size(); ++i) {
// //        std::cerr << "[" << i << "] " << ret[i]->WriteAsCSV('\t') << "\n";
// //    }

//     return ret;
// }

int64_t AlignedMappingResult::QueryId() const {
    return qseq_->abs_id();
}

int64_t AlignedMappingResult::QueryLen() const {
    return qseq_->data().size();
}

std::string AlignedMappingResult::QueryHeader() const {
    return qseq_->header();
}

const mindex::IndexPtr AlignedMappingResult::Index() const {
    return index_;
}

MapperReturnValueBase AlignedMappingResult::ReturnValue() const {
  return rv_;
}

const std::unordered_map<std::string, double>& AlignedMappingResult::Timings() const {
    return timings_;
}




std::vector<int32_t> AlignedMappingResult::CountSimilarMappings_(const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths) {

    std::vector<int32_t> score_counts(12, 0); // For the distribution between 0.00% and 10.0% in 1% steps.

    if (paths.size() == 0) {
        return score_counts;
    }

    auto sorted = paths;

    std::sort(sorted.begin(), sorted.end(),
                [](const std::shared_ptr<raptor::PathAlignment>& a,
                    const std::shared_ptr<raptor::PathAlignment>& b) {
                    return a->path_score() > b->path_score();
                });

    double best_score = sorted[0]->path_score();
    double best_score_f = (double) abs(best_score);

    if (best_score == 0.0) {
        return score_counts;
    }

    for (size_t i = 0; i < sorted.size(); i++) {
        int32_t score = sorted[i]->path_score();
        double score_diff_f = abs((double) (best_score - score));
        double fraction_diff = score_diff_f / best_score_f;
        int32_t bucket = ((int32_t) std::floor(fraction_diff * 100)) + 1;

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

double AlignedMappingResult::CalcBestFractionQueryCovered_(const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths, int64_t qseq_len) {
  double ret = 0.0f;

  if (paths.size() == 0) {
    return ret;
  }

  int32_t best_i = -1;
  for (int32_t i = 0; i < (int32_t) paths.size(); i++) {
    if (best_i == -1 || paths[i]->path_score() > paths[best_i]->path_score()) {
      best_i = i;
    }
  }

  auto& best_path = paths[best_i];

  if (best_path->alns().size() == 0) {
    return ret;
  }

  ret = ((double) (best_path->alns().back()->QueryEnd() - best_path->alns().front()->QueryStart())) / ((double) qseq_len);

  return ret;
}

std::vector<std::shared_ptr<raptor::PathAlignment>> AlignedMappingResult::GenerateFiltered(
    const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths, int32_t bestn,
    double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, double min_idt,
    bool just_sort) {

    std::vector<std::shared_ptr<raptor::PathAlignment>> sorted;

    // First, filter nullptr paths.
    for (size_t i = 0; i < paths.size(); ++i) {
        if (paths[i] == nullptr) {
            continue;
        }
        sorted.emplace_back(paths[i]);
    }

    // Then, sort.
    std::sort(sorted.begin(), sorted.end(),
              [](const std::shared_ptr<raptor::PathAlignment>& a,
                 const std::shared_ptr<raptor::PathAlignment>& b) {
                  return a->path_score() > b->path_score();
              });

    if (sorted.size() == 0) {
        return sorted;
    }

    /*
     * If only sorting is needed, just return.
     */
    if (just_sort) {
        return sorted;
    }

    /*
     * Filter the scores according to the parameters.
     */
    std::vector<std::shared_ptr<raptor::PathAlignment>> ret;
    ret.reserve(sorted.size());

    int32_t best_score = sorted[0]->path_score();

    auto similar_score_count = CountSimilarMappings_(sorted);

    for (size_t i = 0; i < sorted.size(); ++i) {
        if (sorted[i]->entire_alignment() != nullptr &&
            sorted[i]->entire_alignment()->op_counts().identity < min_idt) {
            continue;
        }

        // Safeguard against empty alignments.
        if (sorted[i]->alns().size() == 0) {
            continue;
        }

        if (bestn > 0 && i >= bestn) {
            break;
        }

        if (max_fraction_diff >= 0.0 && i > 0) {
            int32_t score = sorted[i]->path_score();

            double fraction_diff = abs((double)(best_score - score)) / abs((double)best_score);

            if (fraction_diff > max_fraction_diff) {
                break;
            }
        }

        if (sorted[i]->alns().size() > 0 && min_map_len > 0) {
            int32_t map_len =
                sorted[i]->alns().back()->QueryEnd() -
                sorted[i]->alns().front()->QueryStart();
                // TODO: Handling query boundaries needs to be handled better.
            if (map_len < min_map_len) {
                continue;
            }
        }

        ret.emplace_back(sorted[i]);
    }

    // Mapping quality threshold.
    int32_t mapq = raptor::CalcMapqFromScoreCounts(similar_score_count);
    if (mapq < min_mapq) {
        ret.clear();
    }

    ret.shrink_to_fit();

    return ret;
}

void AlignedMappingResult::Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len,
                                   int32_t min_mapq, double min_idt, bool just_sort) {
    path_alignments_ = GenerateFiltered(path_alignments_, bestn, max_fraction_diff, min_map_len,
                                        min_mapq, min_idt, just_sort);
}

}  // namespace raptor
