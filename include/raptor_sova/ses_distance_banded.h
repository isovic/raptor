#ifndef SRC_RAPTOR_SES_BANDED_H_
#define SRC_RAPTOR_SES_BANDED_H_

#include <cstdint>
#include <string>
#include <vector>
#include <limits>

namespace raptor {
namespace ses {

const int32_t MINUS_INF = std::numeric_limits<int32_t>::min() + 10000;  // We need a margin to avoid overflows.

/*
 * These routines implement a banded diff calculation.
 * They are fast, but do not offer traceback.
*/
int32_t AutoBandedSESDistance(const std::string& q, const std::string& t, float max_q_err);
int32_t AutoBandedSESDistance(const std::string& q, const std::string& t);
int32_t AutoBandedSESDistance(const char* q, size_t qlen, const char* t, size_t tlen);
int32_t AutoBandedSESDistance(const char* q, size_t qlen, const char* t, size_t tlen, float max_err);
int32_t BandedSESDistance(const std::string& q, const std::string& t, int32_t band_w);
int32_t BandedSESDistance(const char* q, size_t qlen, const char* t, size_t tlen, int32_t band_w);

}
}

#endif
