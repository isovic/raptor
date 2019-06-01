#include <raptor/raptor_aligner.h>
#include <raptor/path_aligner.h>
#include <sequences/sequence.h>

namespace raptor {

std::unique_ptr<raptor::RaptorAligner> createRaptorAligner(
    const mindex::IndexPtr index, const std::shared_ptr<raptor::ParamsAligner> params,
    std::shared_ptr<raptor::AlignerBase> aligner, std::shared_ptr<raptor::AlignerBase> aligner_gap,
    std::shared_ptr<raptor::AlignerBase> aligner_ext) {
    return std::unique_ptr<raptor::RaptorAligner>(
        new raptor::RaptorAligner(index, params, aligner, aligner_gap, aligner_ext));
}

RaptorAligner::RaptorAligner(const mindex::IndexPtr _index,
                             const std::shared_ptr<raptor::ParamsAligner> _params,
                             std::shared_ptr<raptor::AlignerBase> _aligner,
                             std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                             std::shared_ptr<raptor::AlignerBase> _aligner_ext)
    : index_(_index),
      params_(_params),
      aligner_(_aligner),
      aligner_gap_(_aligner_gap),
      aligner_ext_(_aligner_ext) {}

RaptorAligner::~RaptorAligner() {}

std::shared_ptr<raptor::AlignedMappingResult> RaptorAligner::AlignPaths(
    const mindex::SequencePtr& qseq, const std::vector<std::shared_ptr<raptor::LocalPath>>& paths) {
    auto result = raptor::createRaptorAlignmentResult(qseq, index_);

    raptor::MapperReturnValueBase status = raptor::MapperReturnValueBase::Failed;

    auto path_aligner = raptor::createPathAligner(index_, aligner_, aligner_gap_, aligner_ext_);
    status = raptor::MapperReturnValueBase::OK;

    // Calculate and accumulate the alignments.
    std::vector<std::shared_ptr<raptor::PathAlignment>> alns;
    for (const auto& path : paths) {
        auto aln = path_aligner->Align(qseq, path, !params_->no_extend_alignment, params_);
        alns.emplace_back(aln);
    }

    // Update the alignment IDs (tracking of primary, secondary and supplementary alignments).
    int32_t num_paths = static_cast<int32_t>(alns.size());
    for (int32_t path_id = 0; path_id < num_paths; ++path_id) {
        auto path_aln = alns[path_id];
        if (path_aln == nullptr) {
            continue;
        }
        int32_t num_segments = static_cast<int32_t>(path_aln->alns().size());
        for (int32_t seg_id = 0; seg_id < path_aln->alns().size(); ++seg_id) {
            auto& aln = path_aln->alns()[seg_id];
            aln->path_id(static_cast<int32_t>(path_id));
            aln->num_paths(num_paths);
            aln->segment_id(static_cast<int32_t>(seg_id));
            aln->num_segments(num_segments);
        }
    }

    result->path_alignments(alns);
    result->SetReturnValue(status);

    return result;
}

}  // namespace raptor
