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

// class ParamsMapper;
// inline std::shared_ptr<raptor::ParamsMapper> createParamsMapper();

class ParamsMapper {
   public:
    friend std::shared_ptr<raptor::ParamsMapper> createParamsMapper();
    ~ParamsMapper() = default;

    int64_t verbose_level = 0;
    int64_t debug_qid = -1;
    std::string debug_qname = "";

    int32_t min_qlen = 50;
    bool overlap_skip_self_hits = false;
    bool overlap_single_arc = false;
    bool is_rna = false;
    int32_t diag_margin = 500;
    int32_t seed_join_dist = 5000;
    int32_t min_num_seeds = 3;
    int32_t min_cov_bases = 30;
    int32_t min_dp_score = 60;
    bool score_anchors = false;
    bool graph_score_anchors = false;
    int32_t chain_max_skip = 50;
    int32_t chain_max_predecessors = 500;
    int32_t chain_max_dist = 10000;
    int32_t chain_max_bandwidth = 3000;
    int32_t chain_min_span = 0;
    double chain_penalty_gap = 1.0;
    double chain_penalty_match = 0.20;
    int32_t graph_max_path_edges = 2;
    bool add_symmetric_arcs = false;
    int32_t max_hits = -1;
    bool graph_chain_inversions = false;
    bool ref_and_reads_path_same = false;          // Applicable for overlapping. If true, and "is_overlapper == true", then any overlap with t_id >= q_id will be ignored. Only overlaps with t_id < q_id will be output. Although it may seem redundant, there might be use cases when we want a sanity check that the mapper will find perfect mappings of a set onto itself. Or map only to reverse complement, in which case inputs are the same.
    double graph_allowed_score_diff_frac = 0.95;
    int32_t flank_ext_len = 200;
    bool no_graph_mapping = false;
    bool relabel_secondary_supp = true;
    double min_secondary_to_primary_ratio = 0.80;
    int32_t allowed_suppl_overlap = 0;

   private:
    ParamsMapper() = default;
    ParamsMapper(const ParamsMapper&) = delete;
    ParamsMapper& operator=(const ParamsMapper&) = delete;
};

inline std::shared_ptr<raptor::ParamsMapper> createParamsMapper() {
    return std::shared_ptr<raptor::ParamsMapper>(new raptor::ParamsMapper);
}

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_MAPPER_H_ */
