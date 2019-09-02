/*
 * params_mapper.cc
 *
 *  Created on: Jan 19, 2018
 *      Author: Ivan Sovic
 */

#include <params/params_mapper.h>

namespace raptor {

std::shared_ptr<raptor::ParamsMapper> createParamsMapper() {
    return std::shared_ptr<raptor::ParamsMapper>(new raptor::ParamsMapper);
}

ParamsMapper::ParamsMapper()
    : verbose_level(0),
      debug_qid(-1),
      debug_qname(""),
      min_qlen(50),
      overlap_skip_self_hits(false),
      overlap_single_arc(false),
      is_rna(false),
      diag_margin(500),
      seed_join_dist(5000),
      min_num_seeds(3),
      min_cov_bases(30),
      min_dp_score(60),
      score_anchors(false),
      graph_score_anchors(false),
      chain_max_skip(50),
      chain_max_predecessors(500),
      chain_max_dist(10000),
      chain_max_bandwidth(3000),
      chain_penalty_gap(1.0),
      chain_penalty_match(0.20),
      graph_max_path_edges(2),
      max_hits(-1),
      ref_and_reads_path_same(false),
      graph_allowed_score_diff_frac(0.95),
      flank_ext_len(200),
      no_graph_mapping(false)
{}

}  // namespace raptor