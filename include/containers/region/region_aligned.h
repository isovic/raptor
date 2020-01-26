/*
 * aligned_region.h
 *
 *  Created on: Dec 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTANIERS_ALIGNED_REGION_H_
#define SRC_CONTANIERS_ALIGNED_REGION_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <containers/region/region_base.h>
#include <containers/mapping_env.h>
#include <aligner/alignment_result.h>
#include <log/log_tools.h>
#include <containers/region/region_type.h>

namespace raptor {

class RegionAligned;

typedef std::shared_ptr<raptor::RegionAligned> AlignedRegionPtr;

std::shared_ptr<raptor::RegionAligned> createAlignedRegion(std::shared_ptr<raptor::MappingEnv> env,
                                                       std::shared_ptr<raptor::AlignmentResult> aln,
                                                       int32_t _path_id, int32_t _num_paths,
                                                       int32_t _segment_id, int32_t _num_segments);

/*
 * Wraps the generic raptor::Alignment object in the context of a
 * concrete region where the alignment was performed.
 * This includes the MappingEnv object which specifies the query and
 * target IDs, their orientations, lengths and the beginning of the
 * target in the index.
 * RegionAligned inherits from RegionBase to ensure the interface to
 * provide the start and end locations of the alignment, yield a CIGAR object,
 * and all other values needed for outputting to any mapping/alignment format.
*/
class RegionAligned : public raptor::RegionBase {
 public:
    friend std::shared_ptr<raptor::RegionAligned> createAlignedRegion(std::shared_ptr<raptor::MappingEnv> env,
                                                                  std::shared_ptr<raptor::AlignmentResult> aln,
                                                                  int32_t _path_id, int32_t _num_paths,
                                                                  int32_t _segment_id, int32_t _num_segments);

    ~RegionAligned() = default;

    std::string WriteAsCSV(const char separator) const;

    /*
     * Implementation of base class methods.
    */
    int32_t QueryID() const {
        return env_->q_id;
    }
    bool QueryRev() const {
        return env_->q_rev;
    }
    int32_t QueryStart() const {
        return aln_->position().qstart;
    }
    int32_t QueryEnd() const {
        return aln_->position().qend;
    }
    int32_t QueryLen() const {
        return env_->q_len;
    }
    int32_t QuerySpan() const {
        return (aln_->position().qend - aln_->position().qstart);
    }
    int32_t TargetID() const {
        return env_->t_id;
    }
    bool TargetRev() const {
        return env_->t_rev;
    }
    int32_t TargetStart() const {
        return aln_->position().tstart;
    }
    int32_t TargetEnd() const {
        return aln_->position().tend;
    }
    int32_t TargetLen() const {
        return env_->t_len;
    }
    int32_t TargetSpan() const {
        return (aln_->position().tend - aln_->position().tstart);
    }
    int32_t TargetIndexStart() const {
        return env_->index_t_start;
    }
    int32_t TargetFwdStart() const {
        int32_t rstart = (TargetRev()) ? (TargetLen() - TargetEnd()) : TargetStart();
        return rstart;
    }
    int32_t TargetFwdEnd() const {
        int32_t rend = (TargetRev()) ? (TargetLen()- TargetStart()) : TargetEnd();
        return rend;
    }
    int32_t Score() const {
        return aln_->score();
    }
    int32_t NumSeeds() const {
        return aln_->op_counts().eq;
    }
    int32_t CoveredBasesQuery() const {
        return (QueryEnd() - QueryStart());
    }
    int32_t CoveredBasesTarget() const {
        return (TargetEnd() - TargetStart());
    }
    int32_t EditDistance() const {
        return aln_->edit_dist();
    }
    std::vector<raptor::CigarOp> Cigar() const {
        return aln_->cigar();
    }
    int32_t MatchBases() const {
        return aln_->op_counts().eq;
    }
    int32_t MismatchBases() const {
        return aln_->op_counts().x;
    }
    int32_t InsertionBases() const {
        return aln_->op_counts().i;
    }
    int32_t DeletionBases() const {
        return aln_->op_counts().d;
    }

    int32_t PathId() const {
        return path_id_;
    }
    int32_t PathsNum() const {
        return num_paths_;
    }
    int32_t SegmentId() const {
        return segment_id_;
    }
    int32_t SegmentsNum() const {
        return num_segments_;
    }
    bool IsPrimary() const {
        return (region_priority_ == 0 && region_is_supplementary_ == false);
    }
    bool IsSecondary() const {
        return (region_priority_ > 0);
    }
    bool IsSupplementary() const {
        return region_is_supplementary_;
    }

    const std::unordered_map<std::string, raptor::SamTag>& ExtraTags() const {
        return extra_tags_;
    }

    void SetRegionPriority(int32_t val) {
        region_priority_ = val;
    }
    void SetRegionIsSupplementary(bool val) {
        region_is_supplementary_ = val;
    }
    int32_t GetRegionPriority() const {
        return region_priority_;
    }
    bool GetRegionIsSupplementary() const {
        return region_is_supplementary_;
    }
    raptor::RegionType GetRegionType() const {
        if (region_priority_ == 0 && region_is_supplementary_ == false) {
            return raptor::RegionType::Primary;
        } else if (region_priority_ == 0 && region_is_supplementary_ == true) {
            return raptor::RegionType::PrimarySupplementary;
        } else if (region_priority_ > 0 && region_is_supplementary_ == false) {
            return raptor::RegionType::Secondary;
        } else if (region_priority_ > 0 && region_is_supplementary_ == true) {
            return raptor::RegionType::SecondarySupplementary;
        }
        return raptor::RegionType::Undefined;
    }
    int32_t GetAltRegionCount() const {
        // Number of alternative regions (secondary mappings) covering the same coordinates (either in query or target coords). Not counting this particular region.
        return alt_region_count_;
    }

    ////////////////
    /// Setters. ///
    ////////////////
    void SetCoveredBasesQuery(int32_t val) {
        LOG_ALL("Warning: Base-level setters not implemented in RegionAligned.\n");
    }
    void SetCoveredBasesTarget(int32_t val) {
        LOG_ALL("Warning: Base-level setters not implemented in RegionAligned.\n");
    }
    void SetEditDistance(int32_t val) {
        LOG_ALL("Warning: Base-level setters not implemented in RegionAligned.\n");
    }
    void SetScore(int32_t val) {
        LOG_ALL("Warning: Base-level setters not implemented in RegionAligned.\n");
    }
    void AddTag(const raptor::SamTag& val) {
        extra_tags_[val.name] = val;
    }
    void SetAltRegionCount(int32_t val) {
        alt_region_count_ = val;
    }

    /*
     * Getters.
    */
    const std::shared_ptr<raptor::MappingEnv> env() const { return env_; }
    const std::shared_ptr<raptor::AlignmentResult> aln() const { return aln_; }
    int32_t path_id() const { return path_id_; }
    int32_t num_paths() const { return num_paths_; }
    int32_t segment_id() const { return segment_id_; }
    int32_t num_segments() const { return num_segments_; }

    /*
     * Setters.
    */
    void env(std::shared_ptr<raptor::MappingEnv> env) { env_ = env; }
    void aln(std::shared_ptr<raptor::AlignmentResult> _aln) { aln_ = _aln; }
    void path_id(int32_t _path_id) { path_id_ = _path_id; }
    void num_paths(int32_t _num_paths) { num_paths_ = _num_paths; }
    void segment_id(int32_t _segment_id) { segment_id_ = _segment_id; }
    void num_segments(int32_t _num_segments) { num_segments_ = _num_segments; }

 private:
    RegionAligned(std::shared_ptr<raptor::MappingEnv> env)
                            :   env_(env),
                                aln_(nullptr),
                                path_id_(-1), num_paths_(-1),
                                segment_id_(-1), num_segments_(-1),
                                region_priority_(0), region_is_supplementary_(false),
                                alt_region_count_(0),
                                extra_tags_{}
    {
    }

    RegionAligned(std::shared_ptr<raptor::MappingEnv> env,
                  std::shared_ptr<raptor::AlignmentResult> aln,
                  int32_t _path_id, int32_t _num_paths,
                  int32_t _segment_id, int32_t _num_segments)
                            :   env_(env),
                                aln_(aln),
                                path_id_(_path_id), num_paths_(_num_paths),
                                segment_id_(_segment_id), num_segments_(_num_segments),
                                region_priority_(0), region_is_supplementary_(false),
                                alt_region_count_(0),
                                extra_tags_{}
    {
    }

    std::shared_ptr<raptor::MappingEnv> env_;
    std::shared_ptr<raptor::AlignmentResult> aln_;
    int32_t path_id_;
    int32_t num_paths_;
    int32_t segment_id_;
    int32_t num_segments_;

    int32_t region_priority_;           // Priority 0 is a primary alignment, and > 0 secondary. There can be more than 1 regions of priority "0" but only one is primary, others are supplementary.
    bool region_is_supplementary_;
    int32_t alt_region_count_;

    // If there is any additional data which needs to be available for output, it
    // can be encoded here.
    std::unordered_map<std::string, raptor::SamTag> extra_tags_;
};

}

#endif
