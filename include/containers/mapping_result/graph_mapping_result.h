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
#include <raptor/yield_index.h>
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
        * ones to return.
    */
    static std::vector<std::shared_ptr<raptor::LocalPath>> GenerateFiltered(
                                const std::vector<std::shared_ptr<raptor::LocalPath>>& paths,
                                int32_t bestn,
                                double max_fraction_diff,
                                int32_t min_score_diff_margin,
                                int32_t min_map_len);

    void Filter(int32_t bestn, double max_fraction_diff, int32_t min_map_len, int32_t min_mapq, bool just_sort);

    int64_t qseq_len() const { return qseq_len_; }
    int64_t qseq_id() const { return qseq_id_; }
    raptor::MapperReturnValueBase rv() const { return rv_; }

    void paths(const std::vector<std::shared_ptr<raptor::LocalPath>>& new_paths) {
        paths_ = new_paths;
        // Ensure that the paths are sorted. This is the only place where paths can be set.
        std::sort(paths_.begin(), paths_.end(), [](const std::shared_ptr<raptor::LocalPath>& a, const std::shared_ptr<raptor::LocalPath>& b){ return a->score() > b->score(); } );
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
};

}

#endif
