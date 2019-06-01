/*
 * raptor_aligner.h
 *
 *  Created on: Sep 2, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_ALIGNER_H_
#define SRC_RAPTOR_ALIGNER_H_

#include <stdint.h>
#include <memory>
#include <sequences/sequence_file.h>
#include <aligner/aligner_base.h>
#include <containers/mapping_result/aligned_mapping_result.h>
#include <raptor/index_factory.h>
#include <types/typedefs.h>
#include <graph/local_path.h>
#include <params/params_aligner.h>

namespace raptor {

class RaptorAligner;

/*
 * Aligners can be specified for three types of alignment:
 *  - Aligning an anchor. Anchors are expected to have smaller indels. The "aligner" parameter is
 * used for this.
 *  - Aligning a gap. Gaps between anchors in a chain can have larger gaps. A piecewise alignment is
 * preferable.
 *  - Extending alignment beyond chain ends. A SW based aligner is needed for this task.
 */
std::unique_ptr<raptor::RaptorAligner> createRaptorAligner(
    const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsAligner> params,
    std::shared_ptr<raptor::AlignerBase> aligner, std::shared_ptr<raptor::AlignerBase> aligner_gap,
    std::shared_ptr<raptor::AlignerBase> aligner_ext);

/* Mapper object is intended for mapping and alignment
   of one particular read to the reference.
 */
class RaptorAligner {
   public:
    friend std::unique_ptr<raptor::RaptorAligner> createRaptorAligner(
        const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsAligner> params,
        std::shared_ptr<raptor::AlignerBase> aligner, std::shared_ptr<raptor::AlignerBase> aligner_gap,
        std::shared_ptr<raptor::AlignerBase> aligner_ext);
    ~RaptorAligner();

    std::shared_ptr<raptor::AlignedMappingResult> AlignPaths(
        const mindex::SequencePtr& qseq, const std::vector<std::shared_ptr<raptor::LocalPath>>& paths);

   private:
    RaptorAligner(const mindex::IndexPtr _index, const std::shared_ptr<raptor::ParamsAligner> _params,
                  std::shared_ptr<raptor::AlignerBase> _aligner,
                  std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                  std::shared_ptr<raptor::AlignerBase> _aligner_ext);

    RaptorAligner(const RaptorAligner&) = delete;

    RaptorAligner& operator=(const RaptorAligner&) = delete;

    const mindex::IndexPtr index_;
    const std::shared_ptr<raptor::ParamsAligner> params_;
    const std::shared_ptr<raptor::AlignerBase> aligner_;
    const std::shared_ptr<raptor::AlignerBase> aligner_gap_;
    const std::shared_ptr<raptor::AlignerBase> aligner_ext_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_MAPPER_H_ */
