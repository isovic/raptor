#ifndef SRC_RAPTOR_SES_BANDED_H_
#define SRC_RAPTOR_SES_BANDED_H_

#include <cstdint>
#include <string>
#include <vector>
#include <limits>

namespace raptor {
namespace ses {

const int32_t MINUS_INF = std::numeric_limits<int32_t>::min() + 10000;  // We need a margin to avoid overflows.

class SesResults {
  public:
    int32_t last_q = 0;
    int32_t last_t = 0;
    int32_t last_score = MINUS_INF;
    int32_t max_q = 0;
    int32_t max_t = 0;
    int32_t max_score = MINUS_INF;
    int32_t max_score_diffs = MINUS_INF;
    int32_t diffs = 0;
    bool valid = false;
};

/*
 * These routines implement a banded diff calculation.
 * They are fast, but do not offer traceback.
*/
int32_t BandedSESDistance(const std::string& q, const std::string& t, double maxd_frac, double bandw_frac);
int32_t BandedSESDistance(const char* q, size_t qlen, const char* t, size_t tlen, double maxd_frac, double bandw_frac);
SesResults BandedSESDistanceAdvanced(const char* q, size_t qlen, const char* t, size_t tlen, double maxd_frac, double bandw_frac, int32_t match_score, int32_t indel_penalty);

}
}

#endif
