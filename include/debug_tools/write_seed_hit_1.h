/*
 * write_seed_hit_1.h
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_DEBUG_TOOLS_WRITE_SEED_HIT_1_H_
#define SRC_DEBUG_TOOLS_WRITE_SEED_HIT_1_H_

#include <cstdint>
#include <string>
#include <vector>
#include <index/seed_hit.hpp>
#include <containers/target_hits.hpp>
#include <types/typedefs.h>

namespace raptor {
void WriteSeedHits(const std::string& out_path,
                   const std::vector<mindex::SeedHitPacked>& seed_hits, int32_t seed_len,
                   const std::string& qname, int64_t qlen,
                   const std::string& rname, int64_t rlen);

void WriteTargetHits(const std::string& out_path,
                     const std::vector<raptor::ChainPtr>& target_hits,
                     int32_t seed_len,
                     const std::string& qname, int64_t qlen,
                     const std::string& rname, int64_t rlen);

}

#endif
