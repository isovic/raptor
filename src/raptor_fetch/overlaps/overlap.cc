#include <raptor_fetch/overlaps/overlap.h>
#include <string>

namespace raptor {

std::unique_ptr<raptor::Overlap> createOverlap() {
    return std::unique_ptr<raptor::Overlap>(new raptor::Overlap());
}

std::unique_ptr<raptor::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len) {
    return std::unique_ptr<raptor::Overlap>(new raptor::Overlap(
                                            a_name, b_name, score, identity,
                                            a_rev, a_start, a_end, a_len,
                                            b_rev, b_start, b_end, b_len));
}

std::unique_ptr<raptor::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
            coord_t clip_a_start, coord_t clip_a_end,
            coord_t clip_b_start, coord_t clip_b_end,
            const std::string& type) {

    return std::unique_ptr<raptor::Overlap>(new raptor::Overlap(
                                            a_name, b_name, score, identity,
                                            a_rev, a_start, a_end, a_len,
                                            b_rev, b_start, b_end, b_len,
                                            clip_a_start, clip_a_end,
                                            clip_b_start, clip_b_end,
                                            type));
}


}
