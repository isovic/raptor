/*
 * anchor.h
 *
 *  Created on: Sep 21, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_ANCHOR_H_
#define SRC_ANCHOR_H_

#include <memory>
#include <sstream>
#include <string>
#include <containers/region/region_base.h>
#include <containers/mapping_env.h>
#include <containers/region/region_type.h>

namespace raptor {

class RegionMapped;

std::shared_ptr<raptor::RegionMapped> createRegionMapped(
                int32_t _id,
                int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
                std::shared_ptr<raptor::MappingEnv> env,
                int32_t _qstart, int32_t _qend,
                int32_t _tstart, int32_t _tend,
                int32_t _cov_bases_q, int32_t _cov_bases_t,
                int32_t _num_seeds, int32_t _edit_dist,
                int32_t _score,
                int32_t _path_id, int32_t _num_paths,
                int32_t _segment_id, int32_t _num_segments);

class RegionMapped : public raptor::RegionBase {
public:

    friend std::shared_ptr<raptor::RegionMapped> createRegionMapped(
                    int32_t _id,
                    int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
                    std::shared_ptr<raptor::MappingEnv> env,
                    int32_t _qstart, int32_t _qend,
                    int32_t _tstart, int32_t _tend,
                    int32_t _cov_bases_q, int32_t _cov_bases_t,
                    int32_t _num_seeds, int32_t _edit_dist,
                    int32_t _score,
                    int32_t _path_id, int32_t _num_paths,
                    int32_t _segment_id, int32_t _num_segments);

    std::string Verbose() const;

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
        return qstart_;
    }
    int32_t QueryEnd() const {
        return qend_;
    }
    int32_t QueryLen() const {
        return env_->q_len;
    }
    int32_t QuerySpan() const {
        return (qend_ - qstart_);
    }
    int32_t TargetID() const {
        return env_->t_id;
    }
    bool TargetRev() const {
        return env_->t_rev;
    }
    int32_t TargetStart() const {
        return tstart_;
    }
    int32_t TargetEnd() const {
        return tend_;
    }
    int32_t TargetLen() const {
        return env_->t_len;
    }
    int32_t TargetSpan() const {
        return (tend_ - tstart_);
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
        return score_;
    }
    int32_t NumSeeds() const {
        return num_seeds_;
    }
    int32_t CoveredBasesQuery() const {
        return cov_bases_q_;
    }
    int32_t CoveredBasesTarget() const {
        return cov_bases_t_;
    }
    int32_t EditDistance() const {
        return edit_dist_;
    }
    std::vector<raptor::CigarOp> Cigar() const {
        return std::vector<raptor::CigarOp>();
    }
    int32_t MatchBases() const {
        return 0;
    }
    int32_t MismatchBases() const {
        return 0;
    }
    int32_t InsertionBases() const {
        return 0;
    }
    int32_t DeletionBases() const {
        return 0;
    }
    virtual int32_t MappingQuality() const {
        return mapq_;
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

    /*
    * Returns the score if it was initialized e.g. via alignment,
    * otherwise the query span as a score substitute.
    */
    inline int32_t ScoreOrSpan() const {
        return (score_ >= 0) ? score_ : QuerySpan();
    }

    ////////////////
    /// Setters. ///
    ////////////////
    void SetCoveredBasesQuery(int32_t val) {
        cov_bases_q_ = val;
    }
    void SetCoveredBasesTarget(int32_t val) {
        cov_bases_t_ = val;
    }
    void SetEditDistance(int32_t val) {
        edit_dist_ = val;
    }
    void SetScore(int32_t val) {
        score_ = val;
    }
    void AddTag(const raptor::SamTag& val) {
        extra_tags_[val.name] = val;
    }
    void SetAltRegionCount(int32_t val) {
        alt_region_count_ = val;
    }
    void SetMappingQuality(int32_t val) {
        mapq_ = val;
    }

    /*
     * Getters.
    */
    const std::shared_ptr<raptor::MappingEnv> env() { return env_; }
    int32_t id() const { return id_; }
    int32_t target_hits_id() const { return target_hits_id_; }
    int32_t qstart() const { return qstart_; }
    int32_t qend() const { return qend_; }
    int32_t tstart() const { return tstart_; }
    int32_t tend() const { return tend_; }
    int32_t cov_bases_q() const { return cov_bases_q_; }
    int32_t cov_bases_t() const { return cov_bases_t_; }
    int32_t num_seeds() const { return num_seeds_; }
    int32_t edit_dist() const { return edit_dist_; }
    int32_t score() const { return score_; }
    int32_t path_id() const { return path_id_; }
    int32_t num_paths() const { return num_paths_; }
    int32_t segment_id() const { return segment_id_; }
    int32_t num_segments() const { return num_segments_; }

    /*
     * Setters.
    */
    void env(std::shared_ptr<raptor::MappingEnv> env) { env_ = env; }
    void qstart(int32_t _qstart) { qstart_ = _qstart; }
    void qend(int32_t _qend) { qend_ = _qend; }
    void tstart(int32_t _tstart) { tstart_ = _tstart; }
    void tend(int32_t _tend) { tend_ = _tend; }
    void cov_bases_q(int32_t _cov_bases_q) { cov_bases_q_ = _cov_bases_q; }
    void cov_bases_t(int32_t _cov_bases_t) { cov_bases_t_ = _cov_bases_t; }
    void num_seeds(int32_t _num_seeds) { num_seeds_ = _num_seeds; }
    void edit_dist(int32_t _edit_dist) { edit_dist_ = _edit_dist; }
    void score(int32_t _score) { score_ = _score; }
    void path_id(int32_t _path_id) { path_id_ = _path_id; }
    void num_paths(int32_t _num_paths) { num_paths_ = _num_paths; }
    void segment_id(int32_t _segment_id) { segment_id_ = _segment_id; }
    void num_segments(int32_t _num_segments) { num_segments_ = _num_segments; }

private:
    RegionMapped(
            int32_t _id,
            int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
            std::shared_ptr<raptor::MappingEnv> _env,
            int32_t _qstart, int32_t _qend,
            int32_t _tstart, int32_t _tend);
    RegionMapped(
            int32_t _id,
            int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
            std::shared_ptr<raptor::MappingEnv> _env,
            int32_t _qstart, int32_t _qend,
            int32_t _tstart, int32_t _tend,
            int32_t _cov_bases_q, int32_t _cov_bases_t,
            int32_t _num_seeds, int32_t _edit_dist,
            int32_t _score,
            int32_t _path_id, int32_t _num_paths,
            int32_t _segment_id, int32_t _num_segments);

    // Allow copy constructors in this case, to enable copying of all data.

    int32_t id_;
    int32_t target_hits_id_;
    std::shared_ptr<raptor::MappingEnv> env_;
    int32_t qstart_, qend_;
    int32_t tstart_, tend_;
    int32_t cov_bases_q_, cov_bases_t_;
    int32_t num_seeds_;
    int32_t edit_dist_;
    int32_t score_;    // Score is < 0 if not initialized. Score can represent an arbitrary value such as the number of matching bases.
    int32_t path_id_;
    int32_t num_paths_;
    int32_t segment_id_;
    int32_t num_segments_;

    int32_t region_priority_;           // Priority 0 is a primary alignment, and > 0 secondary. There can be more than 1 regions of priority "0" but only one is primary, others are supplementary.
    bool region_is_supplementary_;
    int32_t alt_region_count_;
    int32_t mapq_;

    // If there is any additional data which needs to be available for output, it
    // can be encoded here.
    std::unordered_map<std::string, raptor::SamTag> extra_tags_;

};

}

#endif
