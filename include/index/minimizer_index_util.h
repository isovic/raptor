/*
 * minimizer_index_util.h
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX_UTIL_H_
#define SRC_MINIMIZER_INDEX_UTIL_H_

#include <string>
#include <cstdint>
#include <index/lookups.h>
#include <index/minimizer_index_types.h>

namespace mindex {

std::string MinimizerKeyToString(minkey_t seed, int32_t k);
}

#endif
