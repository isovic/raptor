#include <raptor_fetch/overlaps/overlap_compact.h>

namespace raptor {

std::unique_ptr<raptor::OverlapCompact> createOverlapCompact() {
    return std::unique_ptr<raptor::OverlapCompact>(new raptor::OverlapCompact());
}

std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(
            seqid_t a_id, seqid_t b_id,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end,
            bool b_rev, coord_t b_start, coord_t b_end) {
    return std::unique_ptr<raptor::OverlapCompact>(new raptor::OverlapCompact(
                                            a_id, b_id, score, identity,
                                            a_rev, a_start, a_end,
                                            b_rev, b_start, b_end));
}

std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(const std::unique_ptr<raptor::OverlapCompact>& ovl) {
    return std::unique_ptr<raptor::OverlapCompact>(new raptor::OverlapCompact(ovl));
}

}
