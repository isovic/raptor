/*
 * params_mapper.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PARAMS_SOVA_MAPPER_H_
#define SRC_RAPTOR_PARAMS_SOVA_MAPPER_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace raptor {
namespace sova {

class ParamsSovaMapper {
   public:
    friend std::shared_ptr<raptor::sova::ParamsSovaMapper> createParamsSovaMapper();
    ~ParamsSovaMapper() = default;

    int64_t verbose_level = 0;
    int64_t debug_qid = -1;
    std::string debug_qname = "";

    int32_t min_qlen = 50;
    bool overlap_skip_self_hits = false;
    bool overlap_single_arc = false;
    int32_t seed_join_dist = 5000;
    int32_t min_num_seeds = 3;
    int32_t min_cov_bases = 30;
    int32_t chain_min_span = 0;
    bool ref_and_reads_path_same = false;          // Applicable for overlapping. If true, and "is_overlapper == true", then any overlap with t_id >= q_id will be ignored. Only overlaps with t_id < q_id will be output. Although it may seem redundant, there might be use cases when we want a sanity check that the mapper will find perfect mappings of a set onto itself. Or map only to reverse complement, in which case inputs are the same.
    int32_t chain_bandwidth = 100;
    double align_bandwidth = 0.01;
    double align_max_diff = 0.03;

   private:
    ParamsSovaMapper() = default;
    ParamsSovaMapper(const ParamsSovaMapper&) = delete;
    ParamsSovaMapper& operator=(const ParamsSovaMapper&) = delete;
};

inline std::shared_ptr<raptor::sova::ParamsSovaMapper> createParamsSovaMapper() {
    return std::shared_ptr<raptor::sova::ParamsSovaMapper>(new raptor::sova::ParamsSovaMapper);
}

}
}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_MAPPER_H_ */
