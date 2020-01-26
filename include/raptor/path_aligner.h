/*
 * path_aligner.h
 *
 *  Created on: Dec 29, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PATH_ALIGNER_H_
#define SRC_RAPTOR_PATH_ALIGNER_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <types/typedefs.h>
#include <sequences/sequence.h>
#include <raptor/graph_mapper.h>
#include <containers/region/region_aligned.h>
#include <containers/path_alignment.h>
#include <graph/local_path.h>
#include <params/params_aligner.h>

namespace raptor {

class AlignmentRegion {
   public:
    // AlignmentRegion() : start(0), end(0), env(nullptr) { }
    // AlignmentRegion(int64_t _start, int64_t _end, std::shared_ptr<raptor::MappingEnv> _env) :
    // start(_start), end(_end), env(_env) { }
    ~AlignmentRegion() = default;
    AlignmentRegion(int64_t _ref_start_pos, int64_t _ref_end_pos, int64_t _spl_start_pos,
                    int64_t _spl_end_pos, std::shared_ptr<raptor::MappingEnv> _env)
        : ref_start_pos(_ref_start_pos),
          ref_end_pos(_ref_end_pos),
          spl_start_pos(_spl_start_pos),
          spl_end_pos(_spl_end_pos),
          env(_env) {}

    int64_t ref_start_pos;  // Start position on the target where the region is from.
    int64_t ref_end_pos;    // End position on the target where the region is from.
    int64_t spl_start_pos;  // Start position in the spliced reference string.
    int64_t spl_end_pos;    // End position in the spliced reference string.
    std::shared_ptr<raptor::MappingEnv> env;
};

class PathAligner;

std::shared_ptr<raptor::PathAligner> createPathAligner(const mindex::IndexPtr _index,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner_ext);

class PathAligner {
   public:
    friend std::shared_ptr<raptor::PathAligner> createPathAligner(
        const mindex::IndexPtr _index, std::shared_ptr<raptor::AlignerBase> _aligner,
        std::shared_ptr<raptor::AlignerBase> _aligner_gap,
        std::shared_ptr<raptor::AlignerBase> _aligner_ext);

    std::shared_ptr<raptor::PathAlignment> Align(const mindex::SequencePtr& qseq,
                                             const std::shared_ptr<raptor::LocalPath> path,
                                            int32_t path_id, int32_t num_paths,
                                             bool use_extend_alignment,
                                             const std::shared_ptr<raptor::ParamsAligner> params);

    static std::shared_ptr<raptor::LocalPath> FlankExtend(const mindex::IndexPtr& index,
                                            const mindex::SequencePtr& qseq,
                                            std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                                            int32_t flank_ext_len,
                                            const std::shared_ptr<raptor::LocalPath>& path);

    ~PathAligner();

   private:
    PathAligner(const PathAligner&) = delete;
    PathAligner& operator=(const PathAligner&) = delete;

    PathAligner(const mindex::IndexPtr _index, std::shared_ptr<raptor::AlignerBase> _aligner,
                std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                std::shared_ptr<raptor::AlignerBase> _aligner_ext);

    bool ComposeImplicitSeq_(const mindex::SequencePtr& qseq, std::shared_ptr<raptor::LocalPath> path,
                             std::shared_ptr<raptor::AnchorGraphEdge> local_edge,
                             int32_t implicit_first_node, int32_t implicit_last_node,
                             std::string& node_seq, AlignmentRegion& node_seq_region) const;

    bool ComposeEdgeSeq_(const mindex::SequencePtr& qseq, std::shared_ptr<raptor::AnchorGraphEdge> edge,
                         std::string& edge_seq,
                         std::vector<AlignmentRegion>& edge_seq_regions) const;

    bool ComposeTargetSequence_(const mindex::SequencePtr& qseq, std::shared_ptr<raptor::LocalPath> path,
                                std::string& ref_seq,
                                std::vector<AlignmentRegion>& ref_seq_regions) const;
    std::vector<std::shared_ptr<raptor::RegionAligned>> SplitAlignment_(
        const mindex::SequencePtr& qseq, const std::shared_ptr<raptor::LocalPath> path,
        const std::string& ref_seq, const std::vector<AlignmentRegion>& ref_seq_regions,
        const std::shared_ptr<raptor::AlignmentResult> entire_alignment,
        double min_idt);

    static std::shared_ptr<raptor::RegionMapped> ExtendAlignmentBack_(
                    const mindex::IndexPtr index,
                    const mindex::SequencePtr& qseq,
                    const std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                    int32_t max_flank_len,  // If < 0, unrestricted.
                    const std::shared_ptr<raptor::RegionMapped> anchor);
    static std::shared_ptr<raptor::RegionMapped> ExtendAlignmentFront_(
                    const mindex::IndexPtr index,
                    const mindex::SequencePtr& qseq,
                    const std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                    int32_t max_flank_len,  // If < 0, unrestricted.
                    const std::shared_ptr<raptor::RegionMapped> anchor);

    const mindex::IndexPtr index_;
    const std::shared_ptr<raptor::AlignerBase> aligner_;
    const std::shared_ptr<raptor::AlignerBase> aligner_gap_;
    const std::shared_ptr<raptor::AlignerBase> aligner_ext_;
};

}  // namespace raptor

#endif
