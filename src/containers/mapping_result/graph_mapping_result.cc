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
std::vector<std::shared_ptr<raptor::RegionBase>> GraphMappingResult::CollectRegions(bool one_hit_per_target) const {

    int32_t num_paths = static_cast<int32_t>(paths().size());

    std::vector<std::pair<int64_t, size_t>> path_scores;
    for (size_t path_id = 0; path_id < paths().size(); ++path_id) {
        const auto& path_aln = paths()[path_id];
        path_scores.emplace_back(std::make_pair(path_aln->score(), path_id));
    }
    std::sort(path_scores.begin(), path_scores.end());
    std::reverse(path_scores.begin(), path_scores.end());

    std::vector<std::shared_ptr<raptor::RegionBase>> ret;
    // Used for filtering multiple hits to the same target.
    std::unordered_map<std::string, int32_t> query_target_pairs;
    // Secondary and SecondarySupplementary alignments.
    for (int64_t path_id = 0; path_id < static_cast<int64_t>(path_scores.size()); ++path_id) {
        auto& curr_path = paths()[std::get<1>(path_scores[path_id])];
        auto merged_path = LocalPathTools::MergeImplicitEdges(curr_path);
        if (merged_path == nullptr) {
            continue;
        }
        int32_t num_segments = static_cast<int32_t>(merged_path->nodes().size());
        // Annotate all chunks as either Secondary or SecondarySupplementary.
        for (size_t aln_id = 0; aln_id < merged_path->nodes().size(); ++aln_id) {
            auto& aln = merged_path->nodes()[aln_id]->data();
            aln->SetRegionPriority(path_id);
            aln->SetRegionIsSupplementary(aln_id > 0);

            if (one_hit_per_target) {
                std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
                if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
                    continue;
                }
                query_target_pairs[pair_name] = 1;
            }
            aln->path_id(static_cast<int32_t>(path_id));
            aln->num_paths(num_paths);
            aln->segment_id(static_cast<int32_t>(aln_id));
            aln->num_segments(num_segments);
            ret.emplace_back(aln);
        }
    }

    raptor::RelabelSupplementary(ret);

//    std::cerr << "[GraphMappingResult::CollectRegions]:\n";
//    for (size_t i = 0; i < ret.size(); ++i) {
//        std::cerr << "[" << i << "] " << ret[i]->WriteAsCSV('\t') << "\n";
//    }

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
                                                int32_t min_map_len,
                                                int32_t min_mapq,
                                                int32_t min_score_diff_margin,
                                                bool just_sort) {

  std::vector<std::shared_ptr<raptor::LocalPath>> sorted;

  // First, filter nullptr paths.
  for (size_t i = 0; i < paths.size(); ++i) {
      if (paths[i] == nullptr) {
          continue;
      }
      sorted.emplace_back(paths[i]);
  }

  // Then, sort.
  std::sort(sorted.begin(), sorted.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );

  if (sorted.size() == 0) {
    return sorted;
  }

  /*
   * If only sorting is needed, just return.
  */
  if (just_sort) {
    return sorted;
  }

  int32_t best_score = sorted[0]->score();

  auto similar_score_count = CountSimilarMappings_(sorted, min_score_diff_margin);

  /*
   * Filter the scores according to the parameters.
  */
  std::vector<std::shared_ptr<raptor::LocalPath>> ret;
  ret.reserve(sorted.size());

  for (size_t i = 0; i < sorted.size(); i++) {
    if (bestn > 0 && i >= bestn) {
      break;
    }

    if (max_fraction_diff >= 0.0 && i > 0) {
      int32_t score = sorted[i]->score();
      int32_t score_diff = abs(best_score - score);

      double fraction_diff = std::abs((double) (best_score - score)) / std::abs((double) best_score);

      if (score_diff > 0 && score_diff > min_score_diff_margin && fraction_diff > max_fraction_diff) {
        break;
      }

    }

    if (sorted[i]->nodes().size() > 0 && min_map_len > 0) {
        int32_t map_len = sorted[i]->nodes().back()->data()->QueryEnd() - sorted[i]->nodes().front()->data()->QueryStart();   // TODO: Handling query boundaries needs to be handled better.
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

void GraphMappingResult::Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, bool just_sort) {
  paths_ = GenerateFiltered(paths_, bestn, max_fraction_diff, min_map_len, min_mapq, index_->k() * 3, just_sort);
}



}
