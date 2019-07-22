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
    int32_t num_paths = paths.size();
    std::vector<std::shared_ptr<raptor::PathAlignment>> alns;
    for (int32_t path_id = 0; path_id < static_cast<int32_t>(paths.size()); ++path_id) {
        const auto& path = paths[path_id];
        auto aln = path_aligner->Align(qseq, path, path_id, num_paths, !params_->no_extend_alignment, params_);
        alns.emplace_back(aln);
    }

    result->path_alignments(alns);
    result->SetReturnValue(status);

    return result;
}

}  // namespace raptor
