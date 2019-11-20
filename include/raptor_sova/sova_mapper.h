/*
 * sova_mapper.h
 *
 *  Created on: Nov 13, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_SOVA_MAPPER_H_
#define SRC_RAPTOR_SOVA_MAPPER_H_

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
#include <raptor_sova/params_sova_mapper.h>

namespace raptor {
namespace sova {

class SovaMapper;

std::unique_ptr<raptor::sova::SovaMapper> createSovaMapper(const mindex::IndexPtr index,
                                         const std::shared_ptr<raptor::sova::ParamsSovaMapper> params);

class SovaMapper {
   public:
    friend std::unique_ptr<raptor::sova::SovaMapper> createSovaMapper(
        const mindex::IndexPtr index, const std::shared_ptr<raptor::sova::ParamsSovaMapper> params);
    friend std::unique_ptr<raptor::sova::SovaMapper> createSovaMapper(
        const mindex::IndexPtr index, const std::shared_ptr<raptor::sova::ParamsSovaMapper> params);
    ~SovaMapper();

    std::shared_ptr<raptor::LinearMappingResult> Map(const mindex::SequencePtr& qseq);

   private:
    SovaMapper(const mindex::IndexPtr index, const std::shared_ptr<raptor::sova::ParamsSovaMapper> params);
    SovaMapper(const SovaMapper&) = delete;
    SovaMapper& operator=(const SovaMapper&) = delete;


    std::vector<std::shared_ptr<raptor::TargetAnchorType>> DeduplicateAnchors_(
        const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& anchors,
        double max_allowed_overlap_frac);

    std::vector<raptor::ChainPtr> LISFilterAndGroupByTarget_(
        std::vector<mindex::SeedHitPacked>& seed_hits, indid_t q_id, ind_t q_len,
        ind_t diag_margin, ind_t min_cov_bases, int32_t min_num_hits);

    template <class T>
    std::vector<std::shared_ptr<raptor::TargetHits<T>>> FilterForOverlapping_(
        const std::vector<std::shared_ptr<raptor::TargetHits<T>>>& target_hits, const mindex::SequencePtr& qseq);

    const mindex::IndexPtr index_;
    const std::shared_ptr<raptor::sova::ParamsSovaMapper> params_;
};

}
} /* namespace raptor */

#endif /* SRC_RAPTOR_MAPPER_H_ */
