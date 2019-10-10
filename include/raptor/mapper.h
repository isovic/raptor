/*
 * mapper.h
 *
 *  Created on: Sep 2, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_MAPPER_H_
#define SRC_RAPTOR_MAPPER_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <sequences/sequence_file.h>
#include <containers/range.h>
#include <containers/target_hits.hpp>
#include <containers/region/region_mapped.h>
#include <containers/mapping_result/linear_mapping_result.h>
#include <raptor/yield_index.h>
#include <types/typedefs.h>
#include <params/params_mapper.h>

namespace raptor {

class Mapper;

std::unique_ptr<raptor::Mapper> createMapper(const mindex::IndexPtr index,
                                         const std::shared_ptr<raptor::ParamsMapper> params);

/* Mapper object is intended for mapping and alignment
   of one particular read to the reference.
 */
class Mapper {
   public:
    friend std::unique_ptr<raptor::Mapper> createMapper(
        const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsMapper> params);
    friend std::unique_ptr<raptor::Mapper> createMapper(
        const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsMapper> params);
    ~Mapper();

    std::shared_ptr<raptor::LinearMappingResult> Map(const mindex::SequencePtr& qseq);

   private:
    Mapper(const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsMapper> params);
    Mapper(const Mapper&) = delete;
    Mapper& operator=(const Mapper&) = delete;


    std::vector<std::shared_ptr<raptor::TargetAnchorType>> DeduplicateAnchors_(
        const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& anchors,
        double max_allowed_overlap_frac);

#ifdef USE_LIS_FILTER
    std::vector<raptor::ChainPtr> LISFilterAndGroupByTarget_(
        std::vector<mindex::SeedHitPacked>& seed_hits, indid_t q_id, ind_t q_len,
        ind_t diag_margin, ind_t min_cov_bases, int32_t min_num_hits);
#endif

    template <class T>
    std::vector<std::shared_ptr<raptor::TargetHits<T>>> FilterForOverlapping_(
        const std::vector<std::shared_ptr<raptor::TargetHits<T>>>& target_hits, const mindex::SequencePtr& qseq);

    const mindex::IndexPtr index_;
    const std::shared_ptr<raptor::ParamsMapper> params_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_MAPPER_H_ */
