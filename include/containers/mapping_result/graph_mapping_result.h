/*
 * graph_mapping_result.h
 *
 *  Created on: Dec 18, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_GRAPH_MAPPING_RESULT_H_
#define SRC_CONTAINERS_GRAPH_MAPPING_RESULT_H_

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>
#include <containers/target_hits.hpp>
#include <raptor/index_factory.h>
#include <types/typedefs.h>
#include <graph/anchor_graph.h>
#include <graph/local_path.h>
#include <utility/mapq.hpp>
#include <containers/mapping_result/mapping_result_base.h>

namespace raptor {

class GraphMappingResult;

/*
 * Using a shared_ptr instead of unique_ptr because this structure might
 * get juggled around and moved from vector to vector.
 * Setting it as a shared_ptr reduced the wrong usage at the expense of
 * a few bytes of space.
*/
std::shared_ptr<raptor::GraphMappingResult> createGraphMappingResult(int64_t qseq_len, int64_t qseq_id, std::string qseq_header, mindex::IndexPtr index);

class GraphMappingResult : public raptor::MappingResultBase {
public:
    friend std::shared_ptr<raptor::GraphMappingResult> createGraphMappingResult(int64_t qseq_len, int64_t qseq_id, std::string qseq_header, mindex::IndexPtr index);

    ~GraphMappingResult() = default;

    // Interface implementation.
    std::vector<std::shared_ptr<raptor::RegionBase>> CollectRegions(bool one_hit_per_target) const;
    int64_t QueryId() const;
    int64_t QueryLen() const;
    std::string QueryHeader() const;
    const mindex::IndexPtr Index() const;
    MapperReturnValueBase ReturnValue() const;
    const std::unordered_map<std::string, double>& Timings() const;

    // Custom methods.
    /*
        * @brief Copies the vector of chains internally, sorts them, and selects the best
        * ones to return. The best scores are selected by the following criteria:
        *   - If just_sort is true, all other filters are ignored.
        *   - If bestn > 0 the bestn number of best scores will be returned and the rest ignored.
        *   - If max_fraction_diff >= 0.0, then all chains with scores within this fraction from the best one will be output and the rest ignored.
        *   - If min_map_len > 0, only chains which span more than min_map_len will be output.
        *  To summarize, any filter with value < 0 will turn that filter off.
        *
        * The min_score_diff_margin is a margin allowed between two scores to consider them almost equal. Used to overcome the numerical rounding and errors.
    */
    static std::vector<std::shared_ptr<raptor::LocalPath>> GenerateFiltered(
                                const std::vector<std::shared_ptr<raptor::LocalPath>>& chains,
                                int32_t bestn, double max_fraction_diff,
                                int32_t min_map_len, int32_t min_mapq, int32_t min_score_diff_margin, bool just_sort);

    void Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, bool just_sort);

    /*
    * @brief Calculates the mapping quality for the alignments, based on the number
    * of objects in the path_alignments_ vector.
    */
    // inline int32_t CalcMapq() const {
    //     int32_t mapq = raptor::CalcMapqFromScoreCounts(all_similar_scores_);
    //     int32_t augment = (fraction_query_covered_ == 1.0) ? (-256) :
    //                       ((int32_t) std::round(log10(100 * (1.0 - fraction_query_covered_))));
    //     mapq = (mapq > 3) ? (mapq - augment) : mapq;
    //     mapq = std::min(mapq, 60);
    //     mapq = std::max(0, mapq);
    //     return mapq;
    // }

    inline int32_t CalcMapq() const {
        int32_t mapq = raptor::CalcMapqFromScoreCounts(all_similar_scores_);

        // int32_t augment = 0;
        // if (mapq > 3 && paths_.size() > 0) {
        //     auto& best_path = paths_[0];
        //     double qspan = (best_path->nodes().back()->data()->QueryEnd() - best_path->nodes().front()->data()->QueryStart());
        //     augment = -10 * log10(((double) best_path->score()) / qspan) / index_->w();
        // }
        // mapq = (mapq > 3) ? (mapq - augment) : mapq;
        // mapq = std::min(mapq, 60);
        // mapq = std::max(0, mapq);

        return mapq;
    }

    int64_t qseq_len() const { return qseq_len_; }
    int64_t qseq_id() const { return qseq_id_; }
    raptor::MapperReturnValueBase rv() const { return rv_; }

    void paths(const std::vector<std::shared_ptr<raptor::LocalPath>>& new_paths) {
        paths_ = new_paths;
        // Ensure that the paths are sorted. This is the only place where paths can be set.
        std::sort(paths_.begin(), paths_.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );
        all_similar_scores_ = CountSimilarMappings_(paths_, index_->params()->k * 3);
        fraction_query_covered_ = CalcBestFractionQueryCovered_(paths_, qseq_len_);
    }

    const std::vector<std::shared_ptr<raptor::LocalPath>>& paths() const {
        return paths_;
    }

    /*
    * Sets the return value of the mapping method.
    */
    void SetReturnValue(MapperReturnValueBase _rv) {
        rv_ = _rv;
    }

    std::unordered_map<std::string, double>& timings() {
        return timings_;
    }

    const std::unordered_map<std::string, double>& timings() const {
        return timings_;
    }

    void timings(const std::unordered_map<std::string, double>& _timings) {
        timings_ = _timings;
    }

private:
    GraphMappingResult(int64_t _qseq_id, int64_t _qseq_len, std::string _qseq_header, mindex::IndexPtr _index);
    GraphMappingResult(const GraphMappingResult&) = delete;
    GraphMappingResult& operator=(const GraphMappingResult&) = delete;

    static std::vector<int32_t> CountSimilarMappings_(const std::vector<std::shared_ptr<raptor::LocalPath>>& paths, int32_t min_score_diff_margin);
    static double CalcBestFractionQueryCovered_(const std::vector<std::shared_ptr<raptor::LocalPath>>& paths, int64_t qseq_len);

    // Interface-specific values.
    int64_t qseq_id_;
    int64_t qseq_len_;
    std::string qseq_header_;
    mindex::IndexPtr index_;
    raptor::MapperReturnValueBase rv_;
    std::unordered_map<std::string, double> timings_;

    // Implementation-specific values.
    std::vector<std::shared_ptr<raptor::LocalPath>> paths_;
    std::vector<int32_t> all_similar_scores_;           // Distribution of similar scores. Needed to augment the mapq.
    double fraction_query_covered_;                     // Fraction of the query covered by the best scoring path. Needed to augment the mapq.
};

}

#endif
