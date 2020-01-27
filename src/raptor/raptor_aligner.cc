#include <raptor/raptor_aligner.h>
#include <raptor/path_aligner.h>
#include <sequences/sequence.h>
#include <containers/mapping_result/mapping_result_common.h>

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

    LabelSupplementaryAndSecondary_(alns, params_->relabel_secondary_supp,
                                    params_->min_secondary_to_primary_ratio,
                                    params_->allowed_suppl_overlap);

    result->path_alignments(alns);
    result->SetReturnValue(status);

    return result;
}

void RaptorAligner::LabelSupplementaryAndSecondary_(
        std::vector<std::shared_ptr<raptor::PathAlignment>>& paths,
        bool do_relabel_sec_supp, double min_sec_to_prim_ratio,
        int32_t allowed_suppl_overlap) {

    int32_t num_paths = static_cast<int32_t>(paths.size());

    // Sort the paths by score, but preserve the original order.
    // Each path can have multiple aligned regions.
    std::vector<std::pair<int64_t, size_t>> path_scores;
    for (size_t path_id = 0; path_id < paths.size(); ++path_id) {
        const auto& path_aln = paths[path_id];
        path_scores.emplace_back(std::make_pair(path_aln->path_score(), path_id));
    }
    std::sort(path_scores.begin(), path_scores.end());
    std::reverse(path_scores.begin(), path_scores.end());

    // Secondary and SecondarySupplementary alignments.
    std::vector<std::shared_ptr<raptor::RegionBase>> sorted_regions;
    for (int64_t i = 0; i < static_cast<int64_t>(path_scores.size()); ++i) {
        auto path_score = std::get<0>(path_scores[i]);
        auto path_id = std::get<1>(path_scores[i]);
        auto& curr_path = paths[path_id];
        int32_t num_segments = static_cast<int32_t>(curr_path->alns().size());
        // Annotate all chunks as either Secondary or SecondarySupplementary.
        for (size_t aln_id = 0; aln_id < curr_path->alns().size(); ++aln_id) {
            auto& aln = curr_path->alns()[aln_id];
            aln->SetRegionPriority(i);
            aln->SetRegionIsSupplementary(aln_id > 0);
            aln->SetAltRegionCount(1);
            aln->SetMappingQuality(255);
            aln->path_id(static_cast<int32_t>(path_id));
            aln->num_paths(num_paths);
            aln->segment_id(static_cast<int32_t>(aln_id));
            aln->num_segments(num_segments);
            sorted_regions.emplace_back(aln);
        }
    }

    if (do_relabel_sec_supp) {
        // The grace distance for mapq scaling is an arbitrary value.
        raptor::RelabelSupplementary(sorted_regions, min_sec_to_prim_ratio, (index_->params()->k + index_->params()->w), allowed_suppl_overlap);
    }
}

}  // namespace raptor
