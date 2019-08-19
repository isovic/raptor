/*
 * aligner_base.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_KSW2_DOUBLE_H_
#define SRC_ALIGNER_ALIGNER_KSW2_DOUBLE_H_

#include <memory>
#include <vector>
#include <aligner/aligner_base.h>
#include <aligner/aligner_containers.h>
#include <aligner/pairwise_penalties.h>
#include <aligner/aligner_util.hpp>
#include <ksw2/ksw2.h>

namespace raptor {

class AlignerKSW2Double;

std::shared_ptr<AlignerBase> createAlignerKSW2Double(const raptor::AlignmentOptions& opt);

class AlignerKSW2Double : public AlignerBase {
   public:
    friend std::shared_ptr<AlignerBase> createAlignerKSW2Double(
        const raptor::AlignmentOptions& opt);

    ~AlignerKSW2Double();

    std::shared_ptr<raptor::AlignmentResult> Global(const char* qseq, int64_t qlen,
                                                    const char* tseq,
                                                    int64_t tlen);  // Global alignment mode.

    std::shared_ptr<raptor::AlignmentResult> Local(const char* qseq, int64_t qlen, const char* tseq,
                                                   int64_t tlen);  // Local alignment mode.

    std::shared_ptr<raptor::AlignmentResult> Semiglobal(
        const char* qseq, int64_t qlen, const char* tseq,
        int64_t tlen);  // Semiglobal alignment mode.

    std::shared_ptr<raptor::AlignmentResult> Extend(const char* qseq, int64_t qlen,
                                                    const char* tseq, int64_t tlen);

    raptor::AlignmentOptions GetAlignmentOptions() {
        return opt_;
    }

   protected:
    AlignerKSW2Double(
        const raptor::AlignmentOptions& opt);  // We don't want users attempting to instantiate
                                               // manually, even though the class is virtual.

   private:
    AlignerKSW2Double(const AlignerKSW2Double&) = delete;              // No copying.
    AlignerKSW2Double& operator=(const AlignerKSW2Double&) = delete;   // No copying.
    AlignerKSW2Double(AlignerKSW2Double&&) = delete;                   // No move constructor.
    AlignerKSW2Double& operator=(const AlignerKSW2Double&&) = delete;  // No copying.

    void KSW2GlobalAlnWrapper_(void* km, const int8_t* qseq_, int qlen, const int8_t* tseq_,
                               int tlen, int8_t m, const int8_t* mat, int8_t q, int8_t e, int8_t q2,
                               int8_t e2, int w, int zdrop, int end_bonus, int flag,
                               ksw_extz_t* ez);

    const raptor::AlignmentOptions& opt_;
};

} /* namespace raptor */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */