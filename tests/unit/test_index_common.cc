#include <tests/unit/test_index_common.h>

#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#define TEST_DEBUG_VERBOSE_

namespace raptor {
namespace unit {

void VerboseSeeds(const std::vector<mindex128_t>& results, const std::vector<mindex128_t>& expected) {
#ifdef TEST_DEBUG_VERBOSE_
    std::cerr << "Result:" << std::endl;
    for (int32_t i = 0; i < results.size(); i++) {
        auto minimizer = mindex::Seed(results[i]);
        std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    }
    std::cerr << std::endl;

    std::cerr << "Expected:" << std::endl;
    for (int32_t i = 0; i < expected.size(); i++) {
        auto minimizer = mindex::Seed(expected[i]);
        std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    }
    std::cerr << std::endl;
#endif
}

}
}
