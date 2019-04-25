/*
 * difflib_edlib.h
 *
 *  Created on: Apr 25, 2019
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_DIFFLIB_EDLIB_H_
#define SRC_ALIGNER_DIFFLIB_EDLIB_H_

#include <memory>
#include <vector>
#include <aligner/difflib_base.h>

namespace raptor {

class DifflibEdlib;

std::shared_ptr<DifflibBase> createDifflibEdlib();

class DifflibEdlib : public DifflibBase {
   public:
    friend std::shared_ptr<DifflibBase> createDifflibEdlib();

    ~DifflibEdlib();

    int64_t CalcDiffs(const char* qseq, int64_t qstart, int64_t qend, int64_t qlen,
                        const char* tseq, bool trev, int64_t tstart, int64_t tend, int64_t tlen) const;

   private:
    DifflibEdlib();
    DifflibEdlib(const DifflibEdlib&) = delete;              // No copying.
    DifflibEdlib& operator=(const DifflibEdlib&) = delete;   // No copying.
    DifflibEdlib(DifflibEdlib&&) = delete;                   // No move constructor.
    DifflibEdlib& operator=(const DifflibEdlib&&) = delete;  // No copying.
};

} /* namespace raptor */

#endif /* SRC_ALIGNER_DIFFLIB_EDLIB_H_ */
