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
#include <types/typedefs.h>

namespace raptor {

/*
 * Performs a dynamic programming based chaining on a set of hits.
 * Converts the packed seed hits into the unpacked form.
 * Does not require the IndexPtr for input, but it will mark target sequence
 * lengths as 0 because of this.
 * However, sometimes it's useful to simply chain the hits, without extra
 * mapping environment information.
*/
std::vector<raptor::ChainPtr> ChainHits(
    const std::vector<mindex::SeedHitPacked>& hits,
    int64_t qseq_abs_id,
    int64_t qseq_len,
    int32_t chain_max_skip,
    int32_t chain_max_predecessors,
    int32_t seed_join_dist,
    int32_t diag_margin,
    int32_t min_num_seeds,
    int32_t min_cov_bases, int32_t min_dp_score, int32_t k);

/*
 * Performs a dynamic programming based chaining on a set of hits.
 * Converts the packed seed hits into the unpacked form.
 * For each chain, a separate TargetHits object will be created, and each
 * TargetHits has it's own MappingEnv describing the target and query
 * information (IDs, lengths and strands).
*/
std::vector<raptor::ChainPtr> ChainHits(
    const std::vector<mindex::SeedHitPacked>& hits,
    const mindex::IndexPtr index,       // Optional. If nullptr, then the target sequence length won't be initialized in the MappingEnv.
    int64_t qseq_abs_id,
    int64_t qseq_len,
    int32_t chain_max_skip,
    int32_t chain_max_predecessors,
    int32_t seed_join_dist,
    int32_t diag_margin,
    int32_t min_num_seeds,
    int32_t min_cov_bases, int32_t min_dp_score, int32_t k);

} /* namespace raptor */

#endif /* SRC_RAPTOR_RAPTOR_H_ */
