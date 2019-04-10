/*
 * params_mapper.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PARAMS_MAPPER_H_
#define SRC_RAPTOR_PARAMS_MAPPER_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace raptor {

class ParamsMapper;

std::shared_ptr<raptor::ParamsMapper> createParamsMapper();

class ParamsMapper {
   public:
    friend std::shared_ptr<raptor::ParamsMapper> createParamsMapper();
    ~ParamsMapper() = default;

    int64_t verbose_level;
    int64_t debug_qid;
    std::string debug_qname;

    int32_t min_qlen;
    bool overlap_skip_self_hits;
    bool overlap_single_arc;
    bool is_rna = false;
    int32_t diag_margin;
    int32_t seed_join_dist;
    int32_t min_num_seeds;
    int32_t min_cov_bases;
    int32_t min_dp_score;
    bool score_anchors;
    bool graph_score_anchors;
    int32_t chain_max_skip;
    int32_t chain_max_predecessors;
    int32_t chain_max_dist;
    int32_t chain_max_bandwidth;
    double chain_penalty_gap;
    double chain_penalty_match;
    int32_t graph_max_path_edges;
    bool add_symmetric_arcs;
    int32_t max_hits;
    bool graph_chain_inversions;
    bool ref_and_reads_path_same;          // Applicable for overlapping. If true, and "is_overlapper == true", then any overlap with t_id >= q_id will be ignored. Only overlaps with t_id < q_id will be output. Although it may seem redundant, there might be use cases when we want a sanity check that the mapper will find perfect mappings of a set onto itself. Or map only to reverse complement, in which case inputs are the same.
    double graph_allowed_score_diff_frac;
    int32_t flank_ext_len;

   private:
    ParamsMapper();
    ParamsMapper(const ParamsMapper&) = delete;
    ParamsMapper& operator=(const ParamsMapper&) = delete;
};

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_MAPPER_H_ */
