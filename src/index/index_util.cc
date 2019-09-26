#include <index/index_util.h>

#include <algorithm>
#include <sstream>

namespace mindex {

std::string MinimizerKeyToString(minkey_t seed, int32_t k) {
    std::stringstream ss;
    for (int32_t i = 0; i < k; i++) {
        ss << (char)twobit_to_nuc[seed & 0x03];
        seed = seed >> 2;
    }
    std::string seed_str;
    seed_str = ss.str();
    std::reverse(seed_str.begin(), seed_str.end());
    return seed_str;
}

}  // namespace mindex
