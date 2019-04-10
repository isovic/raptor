/*
 * mapper_tools.h
 *
 *  Created on: Jan 27, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_MAPPER_TOOLS_H_
#define SRC_RAPTOR_MAPPER_TOOLS_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <types/typedefs.h>
#include <index/minimizer_index.h>

namespace raptor {
namespace mapper {

void CalcHitCoverage(const std::vector<mindex::MinimizerHitPacked>& hits, int32_t seed_len,
                    int32_t hits_begin, int32_t hits_end, int32_t& cov_bases_q,
                    int32_t& cov_bases_r);

std::vector<std::shared_ptr<raptor::TargetAnchorType>> MakeAnchors(
    const std::vector<std::shared_ptr<raptor::TargetHits<mindex::MinimizerHitPacked>>>& target_hits);


std::vector<std::shared_ptr<raptor::TargetHits<mindex::MinimizerHitPacked>>> GroupTargetSeedHits(
                    std::vector<mindex::MinimizerHitPacked> seed_hits,  // Copy.
                    int32_t k,
                    int32_t qid,
                    int32_t qlen);

}
}

#endif