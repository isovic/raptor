#include <containers/mapping_result/aligned_mapping_result.h>
#include <aligner/aligner_util.hpp>
#include <utility/revcmp.hpp>
#include <math.h>
#include <containers/mapping_result/mapping_result_common.h>

namespace raptor {

std::shared_ptr<raptor::AlignedMappingResult> createRaptorAlignmentResult(const mindex::SequencePtr& qseq,
                                                                       mindex::IndexPtr index) {
    return std::shared_ptr<raptor::AlignedMappingResult>(new raptor::AlignedMappingResult(qseq, index));
}

AlignedMappingResult::AlignedMappingResult(const mindex::SequencePtr& qseq, mindex::IndexPtr _index)
    : qseq_(qseq),
      index_(_index),
      rv_(raptor::MapperReturnValueBase::NotRunYet),
      path_alignments_() {}

// Interface implementation.
std::vector<std::shared_ptr<raptor::RegionBase>> AlignedMappingResult::CollectRegions(bool one_hit_per_target) const {
    // // Sort paths by score.
    // auto sorted = paths_;
    // std::sort(sorted.begin(), sorted.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );

    std::vector<std::shared_ptr<raptor::RegionBase>> ret;
    std::vector<std::pair<std::shared_ptr<raptor::RegionBase>, int32_t>> secondary;

    // Used for filtering multiple hits to the same target.
    std::unordered_map<std::string, int8_t> query_target_pairs;

    // Collect the PRIMARY and SECONDARY alignments separately.
    for (size_t i = 0; i < path_alignments_.size(); ++i) {
        const auto& path = path_alignments_[i];
        if (path == nullptr) {
            continue;
        }
        // Annotate all chunks as either Secondary or SecondarySupplementary.
        for (size_t aln_id = 0; aln_id < path->alns().size(); ++aln_id) {
            const auto& aln = path->alns()[aln_id];
            if (one_hit_per_target) {
                std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
                if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
                    continue;
                }
                query_target_pairs[pair_name] = 1;
            }
            // Priority of 0 are the primary alignments, and < 0 should be filtered.
            auto priority = aln->GetRegionPriority();
            if (priority == 0) {
                ret.emplace_back(aln);
            } else if (priority > 0) {
                secondary.emplace_back(std::make_pair(aln, path->path_score()));
            }
        }
    }

    // For now, output absolutely any alignment. Filtering should have happened earlier.
    for (const auto& vals: secondary) {
        ret.emplace_back(std::get<0>(vals));
    }

    return ret;
}

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
                                                const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths,
                                                int32_t bestn,
                                                double max_fraction_diff,
                                                int32_t min_map_len,
                                                double min_idt
                                                ) {

    /*
     * This function keeps any path which has at least one region within the bestn
     * priority. Also, any path within the max_fraction_diff from the best scoring
     * path is kept even if it's beyond bestn priority.
     *
    */

    std::vector<std::shared_ptr<raptor::PathAlignment>> ret;

    if (paths.empty()) {
        return ret;
    }

    // Find the best score.
    int64_t best_score = std::numeric_limits<int32_t>::lowest();
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        if (path == nullptr) { continue; }
        best_score = std::max(best_score, path->path_score());
    }
    // Keep any score above this fraction difference even if it's above bestn regions.
    int64_t min_score = best_score * (1.0 - max_fraction_diff);
    bool ret_contains_primary = false;

    // Filter the bestn alternative paths.
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        if (path == nullptr) { continue; }
        if (path->alns().empty()) { continue; }

        // Keep path if any region has priority less than bestn.
        bool keep_path = false;
        bool contains_primary = false;
        for (size_t aln_id = 0; aln_id < path->alns().size(); ++aln_id) {
            auto& aln = path->alns()[aln_id];
            auto priority = aln->GetRegionPriority();
            if (priority == 0) {
                contains_primary = true;
            }
            // If bestn == 1, then no secondary alignments should explicitly be taken.
            if (bestn <= 0 || (bestn == 1 && priority == 0) || (bestn > 1 && priority >= 0 && priority < bestn)) {
                keep_path = true;
            }
        }

        // Keep any score above this fraction difference, but only if bestn != 1.
        // If bestn == 1, then only primary alignments should be kept.
        if (path->path_score() >= min_score && (bestn != 1 || (bestn == 1 && contains_primary))) {
            keep_path = true;
        }

        // Filter paths that are too short.
        int32_t qspan = path->alns().front()->QueryEnd() - path->alns().front()->QueryStart();
        if (qspan < min_map_len) {
            keep_path = false;
        }

        if (path->entire_alignment() != nullptr &&
            path->entire_alignment()->op_counts().identity_min < min_idt) {
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

void AlignedMappingResult::Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len,
                                   double min_idt, bool just_sort) {
    path_alignments_ = GenerateFiltered(path_alignments_, bestn, max_fraction_diff, min_map_len, min_idt);
}

}  // namespace raptor
