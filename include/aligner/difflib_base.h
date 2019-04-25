/*
 * difflib_base.h
 *
 *  Created on: Apr 25, 2019
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_DIFFLIB_BASE_H_
#define SRC_ALIGNER_DIFFLIB_BASE_H_

#include <memory>
#include <vector>

namespace raptor {

class DifflibBase {
   public:
    virtual ~DifflibBase() {}

    virtual int64_t CalcDiffs(const char* qseq, int64_t qstart, int64_t qend, int64_t qlen,
                                const char* tseq, bool trev, int64_t tstart, int64_t tend, int64_t tlen) const = 0;
};

} /* namespace raptor */

#endif /* SRC_ALIGNER_DIFFLIB_BASE_H_ */
