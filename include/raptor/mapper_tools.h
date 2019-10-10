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

void CalcHitCoverage(const std::vector<mindex::SeedHitPacked>& hits, int32_t seed_len,
                    int32_t hits_begin, int32_t hits_end, int32_t& cov_bases_q,
                    int32_t& cov_bases_r);

std::vector<std::shared_ptr<raptor::TargetAnchorType>> MakeAnchors(
    const std::vector<raptor::ChainPtr>& target_hits);


std::vector<raptor::ChainPtr> GroupTargetSeedHits(
                    std::vector<mindex::SeedHitPacked> seed_hits,  // Copy.
                    int32_t k,
                    int32_t qid,
                    int32_t qlen);


/*
* Uses alignment (Edlib edit distance) to estimate the number of matches in anchors.
* This is then used for DP chaining. I used the term "estimate" because Edlib returns
* the edit distance instead of the matches we need, so I subtracted the edit distance
* from the query length. This produces a pessimistic score because indels should
* not be counted. In this way, this actually provides an estimate of accuracy (but not
* normalized).
*/
int32_t CalcMatchRate(const mindex::SequencePtr& qseq, mindex::IndexPtr index,
                    const std::shared_ptr<raptor::RegionMapped>& anchor);

}
}

#endif