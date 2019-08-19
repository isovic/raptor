/*
 * aligner_base.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_BASE_H_
#define SRC_ALIGNER_ALIGNER_BASE_H_

#include <memory>
#include <vector>
#include <aligner/aligner_containers.h>
#include <aligner/pairwise_penalties.h>

namespace raptor {

class AlignerBase {
   public:
    virtual ~AlignerBase() {}

    // virtual AlignmentReturnValue Align(const char* q, int64_t qlen, const char* t, int64_t tlen,
    // AlignmentType type) = 0;      // Selects the alignment mode based on a parameter.

    virtual std::shared_ptr<raptor::AlignmentResult> Global(
        const char* qseq, int64_t qlen, const char* tseq,
        int64_t tlen) = 0;  // Global alignment mode.

    virtual std::shared_ptr<raptor::AlignmentResult> Local(
        const char* qseq, int64_t qlen, const char* tseq,
        int64_t tlen) = 0;  // Local alignment mode.

    virtual std::shared_ptr<raptor::AlignmentResult> Semiglobal(
        const char* qseq, int64_t qlen, const char* tseq,
        int64_t tlen) = 0;  // Semiglobal alignment mode.

    virtual std::shared_ptr<raptor::AlignmentResult> Extend(
        const char* qseq, int64_t qlen, const char* tseq,
        int64_t tlen) = 0;  // Extend alignment mode. Does not necessarily
                            //  produce CIGAR,but generate max alignment coords
    virtual raptor::AlignmentOptions GetAlignmentOptions() = 0;
};

} /* namespace raptor */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */
