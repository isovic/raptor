/*
 * typedefs.h
 *
 *  Created on: Sep 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_TYPES_TYPEDEFS_
#define SRC_TYPES_TYPEDEFS_

#include <memory>
#include <vector>
#include <containers/target_hits.hpp>
#include <containers/region/region_mapped.h>
#include <aligner/aligner_base.h>
// #include <index/minimizer_index.h>
#include <index/minimizer_hit.hpp>
#include <graph/segment_graph.h>

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

namespace raptor {

typedef std::shared_ptr<raptor::AlignerBase> AlignerPtr;
typedef std::shared_ptr<raptor::SegmentGraph> GraphPtr;
typedef std::shared_ptr<raptor::RegionMapped> AnchorPtr;
typedef raptor::TargetHits<raptor::AnchorPtr> TargetAnchorType;
typedef std::shared_ptr<raptor::TargetHits<raptor::AnchorPtr>> TargetAnchorPtr;
typedef std::vector<TargetAnchorPtr> TargetAnchorPtrVector;

typedef raptor::TargetHits<mindex::SeedHitPacked> ChainType;
typedef std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>> ChainPtr;

}

#endif
