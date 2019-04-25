#include <aligner/difflib_edlib.h>

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <lib/edlib.h>
#include <utility/revcmp.hpp>

namespace raptor {

std::shared_ptr<DifflibBase> createDifflibEdlib() {
    return std::shared_ptr<DifflibBase>(new DifflibEdlib());
}

DifflibEdlib::DifflibEdlib() = default;

DifflibEdlib::~DifflibEdlib() = default;

int64_t DifflibEdlib::CalcDiffs(const char* qseq, int64_t qstart, int64_t qend, int64_t qlen,
                                const char* tseq, bool trev, int64_t tstart, int64_t tend, int64_t tlen) const {

    int64_t diffs = 0;

    EdlibAlignTask task = EDLIB_TASK_DISTANCE;
    EdlibAlignMode edlib_mode = EDLIB_MODE_NW;

    if (trev) {
        std::string target = raptor::ReverseComplement((const char*)(tseq + tstart), (tend - tstart));
        EdlibAlignResult result =
            edlibAlign((const char*)(qseq + qstart), (qend - qstart),
                        target.c_str(), target.size(),
                        edlibNewAlignConfig(-1, edlib_mode, task));
        diffs = result.editDistance;
        edlibFreeAlignResult(result);

    } else {
        EdlibAlignResult result =
            edlibAlign((const char*)(qseq + qstart), (qend - qstart),
                        (tseq + tstart), (tend - tstart),
                        edlibNewAlignConfig(-1, edlib_mode, task));
        diffs = result.editDistance;
        edlibFreeAlignResult(result);

    }

    return diffs;
}

}  // namespace raptor
