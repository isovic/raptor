#include <index/minimizer_generator.h>

namespace mindex {

std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int32_t window_len) {
    return std::unique_ptr<MinimizerGenerator>(new MinimizerGenerator(window_len));
}

MinimizerGenerator::MinimizerGenerator(int32_t window_len)
    : window_len_(window_len), num_added_keys_(0) {}

MinimizerGenerator::~MinimizerGenerator() {}

int MinimizerGenerator::Yield(const mindex::Minimizer& seed_in,
                              std::vector<mindex::Minimizer>& seeds_out) {
    if (num_added_keys_ < window_len_) {
        // Remove smaller elements if any.
        while ((!q_.empty()) && seed_in.key <= std::get<0>(q_.back()).key) {
            q_.pop_back();
        }

    } else {
        // Return one or more elements of the window
        // with the minimal key value.
        seeds_out.clear();
        if (q_.size() > 0) {
            minkey_t key = std::get<0>(q_.front()).key;
            auto it_next = q_.begin();
            while (it_next != q_.end() && std::get<0>(*(it_next)).key == key) {
                seeds_out.emplace_back(std::get<0>(*(it_next)));
                ++it_next;
            }
        }

        // Remove the elements which are out of this window.
        while ((!q_.empty()) && std::get<1>(q_.front()) <= (num_added_keys_ - window_len_)) {
            q_.pop_front();
        }
        // Remove smaller elements if any.
        while ((!q_.empty()) && seed_in.key <= std::get<0>(q_.back()).key) {
            q_.pop_back();
        }
    }

    q_.push_back(std::make_pair(seed_in, num_added_keys_));
    num_added_keys_ += 1;

    if (num_added_keys_ <= window_len_) {
        return 2;
    }

    return 0;
}

}  // namespace mindex
