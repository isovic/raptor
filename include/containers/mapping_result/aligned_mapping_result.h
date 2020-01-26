/*
 * raptor_alignment_results.h
 *
 *  Created on: Jan 03, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_RAPTOR_ALIGNMENT_RESULT_H_
#define SRC_CONTAINERS_RAPTOR_ALIGNMENT_RESULT_H_

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>
#include <containers/target_hits.hpp>
#include <raptor/yield_index.h>
#include <types/typedefs.h>
#include <containers/path_alignment.h>
#include <utility/mapq.hpp>
#include <sequences/sequence.h>
#include <containers/mapping_result/mapping_result_base.h>

namespace raptor {

class AlignedMappingResult;

typedef std::shared_ptr<raptor::AlignedMappingResult> RaptorAlignmentResultPtr;

/*
 * Using a shared_ptr instead of unique_ptr because this structure might
 * get juggled around and moved from vector to vector.
 * Setting it as a shared_ptr reduced the wrong usage at the expense of
 * a few bytes of space.
*/
std::shared_ptr<raptor::AlignedMappingResult> createRaptorAlignmentResult(const mindex::SequencePtr& qseq, mindex::IndexPtr index);

class AlignedMappingResult : public raptor::MappingResultBase {
public:
    friend std::shared_ptr<raptor::AlignedMappingResult> createRaptorAlignmentResult(const mindex::SequencePtr& qseq, mindex::IndexPtr index);

    ~AlignedMappingResult() = default;

    // Interface implementation.
    std::vector<std::shared_ptr<raptor::RegionBase>> CollectRegions(bool one_hit_per_target, bool do_relabel_sec_supp) const;
    int64_t QueryId() const;
    int64_t QueryLen() const;
    std::string QueryHeader() const;
    const mindex::IndexPtr Index() const;
    MapperReturnValueBase ReturnValue() const;
    const std::unordered_map<std::string, double>& Timings() const;

    // Implementation-specific interfaces.
    raptor::MapperReturnValueBase rv() const { return rv_; }

    const mindex::IndexPtr index() const { return index_; }

    std::vector<std::shared_ptr<raptor::PathAlignment>>& path_alignments() {
        return path_alignments_;
    }

    const std::vector<std::shared_ptr<raptor::PathAlignment>>& path_alignments() const {
        return path_alignments_;
    }

    void path_alignments(const std::vector<std::shared_ptr<raptor::PathAlignment>>& new_alns) {
        path_alignments_ = new_alns;

        // Ensure paths are sorted by path_score.
        std::sort(path_alignments_.begin(), path_alignments_.end(),
                [](const std::shared_ptr<raptor::PathAlignment>& a,
                    const std::shared_ptr<raptor::PathAlignment>& b) {
                    return a->path_score() > b->path_score();
                });

        all_similar_scores_ = CountSimilarMappings_(path_alignments_);
        fraction_query_covered_ = CalcBestFractionQueryCovered_(path_alignments_, qseq_->data().size());
    }

    std::vector<std::shared_ptr<raptor::PathAlignment>> GenerateFiltered(
                                                    const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths,
                                                    int32_t bestn,
                                                    double max_fraction_diff,
                                                    int32_t min_map_len,
                                                    double min_idt);

    void Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, double min_idt, bool just_sort);

    /*
    * @brief Calculates the mapping quality for the alignments, based on the number
    * of objects in the path_alignments_ vector.
    */
    inline int32_t CalcMapq() const {
        int32_t mapq = raptor::CalcMapqFromScoreCounts(all_similar_scores_);

        // int32_t augment = 0;
        // if (mapq > 3 && path_alignments_.size() > 0) {
        //   auto& best_path = path_alignments_[0];
        //   // augment = -10 * log10(best_path->entire_alignment()->op_counts().identity / 100.0);
        //   double qspan = (best_path->alns().back()->QueryEnd() - best_path->alns().front()->QueryStart());
        //   augment = -10 * log10(((double) best_path->path_score()) / qspan);
        //   // printf ("qspan = %f\npath_score = %d\n", qspan, best_path->path_score());
        // }
        // // int32_t augment = (fraction_query_covered_ == 1.0) ? (-256) :
        // //                   ((int32_t) std::round(log10(100 * (1.0 - fraction_query_covered_))));
        // mapq = (mapq > 3) ? (mapq - augment) : mapq;
        // mapq = std::min(mapq, 60);
        // mapq = std::max(0, mapq);

        return mapq;
    }

    /*
    * Sets the return value of the mapping method.
    */
    void SetReturnValue(MapperReturnValueBase _rv) {
        rv_ = _rv;
    }

private:
    AlignedMappingResult(const AlignedMappingResult&) = delete;
    AlignedMappingResult& operator=(const AlignedMappingResult&) = delete;
    AlignedMappingResult(const mindex::SequencePtr& qseq, mindex::IndexPtr _index); // Prevent creating objects in any way other than RAII.

    static std::vector<int32_t> CountSimilarMappings_(const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths);
    static double CalcBestFractionQueryCovered_(const std::vector<std::shared_ptr<raptor::PathAlignment>>& paths, int64_t qseq_len);

    const mindex::SequencePtr& qseq_;
    mindex::IndexPtr index_;
    raptor::MapperReturnValueBase rv_;
    std::vector<std::shared_ptr<raptor::PathAlignment>> path_alignments_;
    std::vector<int32_t> all_similar_scores_;         // Distribution of similar scores. Needed to augment the mapq.
    double fraction_query_covered_;                   // Fraction of the query covered by the best scoring path. Needed to augment the mapq.
    // Time measurements, useful for debugging.
    std::unordered_map<std::string, double> timings_;
};

}

#endif
