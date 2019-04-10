#include <aligner/alignment_result.h>
#include <sstream>
#include <string>

namespace raptor {

std::shared_ptr<raptor::AlignmentResult> createAlignmentResult() {
    return std::shared_ptr<raptor::AlignmentResult>(new raptor::AlignmentResult());
}

AlignmentResult::AlignmentResult()
    : score_(0),
      edit_dist_(0),
      position_(),
      max_score_(raptor::LARGE_NEGATIVE_INT64),
      max_q_pos_(-1),
      max_t_pos_(-1),
      final_band_(-1),
      status_(AlignmentReturnValue::AlignmentNotPerformed) {}

AlignmentResult::~AlignmentResult() {}

std::string AlignmentResult::Verbose() const {
    std::ostringstream oss;

    oss << "q = [" << position().qstart << ", " << position().qend << "]"
        << ", t = [" << position().tstart << ", " << position().tend << "]"
        << ", max_q_pos = " << max_q_pos() << ", max_t_pos = " << max_t_pos()
        << ", max_score = " << max_score() << ", score = " << score()
        << ", edit_dist = " << edit_dist() << ", cigar = ";

    for (auto& op : cigar()) {
        oss << op.count << op.op;
    }

    return oss.str();
}

}  // namespace raptor
