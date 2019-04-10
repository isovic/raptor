/*
 * aligner_edlib.h
 *
 *  Created on: Jan 7, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_EDLIB_H_
#define SRC_ALIGNER_ALIGNER_EDLIB_H_

#include <memory>
#include <vector>
#include <aligner/aligner_base.h>
#include <aligner/aligner_containers.h>
#include <aligner/pairwise_penalties.h>
#include <aligner/aligner_util.hpp>
#include <lib/edlib.h>

namespace raptor {

class AlignerEdlib;

std::shared_ptr<AlignerBase> createAlignerEdlib(const raptor::AlignmentOptions& opt);

class AlignerEdlib : public AlignerBase {
   public:
    friend std::shared_ptr<AlignerBase> createAlignerEdlib(const raptor::AlignmentOptions& opt);

    ~AlignerEdlib();

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

    static std::shared_ptr<raptor::AlignmentResult> RunEdlib(const char* qseq, int64_t qlen,
                                                             const char* tseq, int64_t tlen,
                                                             const raptor::AlignmentOptions& opt,
                                                             EdlibAlignMode edlib_mode,
                                                             EdlibAlignTask edlib_task);

   protected:
    AlignerEdlib(
        const raptor::AlignmentOptions& opt);  // We don't want users attempting to instantiate
                                               // manually, even though the class is virtual.

   private:
    AlignerEdlib(const AlignerEdlib&) = delete;              // No copying.
    AlignerEdlib& operator=(const AlignerEdlib&) = delete;   // No copying.
    AlignerEdlib(AlignerEdlib&&) = delete;                   // No move constructor.
    AlignerEdlib& operator=(const AlignerEdlib&&) = delete;  // No copying.

    const raptor::AlignmentOptions& opt_;
};

} /* namespace raptor */

#endif /* SRC_ALIGNER_ALIGNER_BASE_H_ */