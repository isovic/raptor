/*
 * types.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_CONSTS_H_
#define SRC_ALIGNER_CONSTS_H_

#include <stdint.h>
#include <limits>

namespace raptor {

static constexpr int64_t LARGE_NEGATIVE_INT64 = std::numeric_limits<int64_t>::min() + 10000;
}

#endif
