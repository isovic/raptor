/*
 * dp_chain.h
 *
 *  Created on: Dec 16, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_DP_CHAIN_H_
#define SRC_RAPTOR_DP_CHAIN_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <sequences/sequence_file.h>
#include <containers/target_hits.hpp>
#include <containers/region/region_mapped.h>
#include <params/params_mapper.h>
#include <index/minimizer_index.h>

namespace raptor {

std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>> ChainHits(
    const mindex::IndexPtr index,
    int64_t qseq_abs_id,
    int64_t qseq_len,
    int32_t chain_max_skip,
    int32_t chain_max_predecessors,
    int32_t seed_join_dist,
    int32_t diag_margin,
    int32_t min_num_seeds,
    int32_t min_cov_bases, int32_t min_dp_score, int32_t k,
    const std::vector<mindex::SeedHitPacked>& hits);

// std::vector<mindex::SeedHitPacked> ChainTargetHits2(
//     const SingleSequence& qseq,
//     const std::vector<mindex::SeedHitPacked>& hits,
//     const std::shared_ptr<raptor::ParamsMapper> params, int32_t k);

} /* namespace raptor */

#endif /* SRC_RAPTOR_RAPTOR_H_ */
