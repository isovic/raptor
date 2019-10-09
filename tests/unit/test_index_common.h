#ifndef TESTS_UNIT_INDEX_COMMON_H_
#define TESTS_UNIT_INDEX_COMMON_H_

#include <index/index_params.h>
#include <index/index_factory.h>
#include <index/seed.hpp>
#include <vector>

namespace raptor {
namespace unit {

void VerboseSeeds(const std::vector<mindex128_t>& results, const std::vector<mindex128_t>& expected);

}
}

#endif
