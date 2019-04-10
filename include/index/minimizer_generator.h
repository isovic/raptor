#ifndef MINIMIZER_GENERATOR_2_H_
#define MINIMIZER_GENERATOR_2_H_

#include <cstdint>
#include <deque>
#include <memory>
#include <vector>
#include <tuple>
#include <index/minimizer_index_types.h>
#include <index/minimizer.hpp>

namespace mindex {

class MinimizerGenerator;

std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int32_t window_len);

class MinimizerGenerator {
   public:
    friend std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int32_t window_len);
    ~MinimizerGenerator();

    /* Pushes the new seed into the queue and returns one or more elements with the minimum key
     * value within the given window. Returns 0 if the generator has enough seeds to fill a window,
     * which means that seed_out is output. Returns 1 if the queue has not yet been filled up to
     * window_len size.
     */
    int Yield(const mindex::Minimizer& seed_in, std::vector<mindex::Minimizer>& seeds_out);

   private:
    MinimizerGenerator(int32_t window_len);

    int32_t window_len_;
    size_t num_added_keys_;
    std::deque<std::pair<mindex::Minimizer, int64_t>> q_;
};

}  // namespace mindex

#endif
