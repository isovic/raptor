#include <containers/region/region_mapped.h>

namespace raptor {

std::shared_ptr<raptor::RegionMapped> createRegionMapped(
                int32_t _id,                // ID of this current anchor.
                int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
                std::shared_ptr<raptor::MappingEnv> _env,
                int32_t _qstart, int32_t _qend,
                int32_t _tstart, int32_t _tend,
                int32_t _cov_bases_q, int32_t _cov_bases_t,
                int32_t _num_seeds, int32_t _score,
                int32_t _path_id, int32_t _num_paths,
                int32_t _segment_id, int32_t _num_segments) {

    return std::shared_ptr<raptor::RegionMapped>(new raptor::RegionMapped(
                                            _id,
                                            _target_hits_id,
                                            _env,
                                            _qstart, _qend,
                                            _tstart, _tend,
                                            _cov_bases_q, _cov_bases_t,
                                            _num_seeds, _score,
                                            _path_id, _num_paths,
                                            _segment_id, _num_segments));

}

RegionMapped::RegionMapped(
        int32_t _id,                // ID of this current anchor.
        int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
        std::shared_ptr<raptor::MappingEnv> _env,
        int32_t _qstart, int32_t _qend,
        int32_t _tstart, int32_t _tend)
        :   id_(_id),
            target_hits_id_(_target_hits_id),
            env_(_env),
            qstart_(_qstart), qend_(_qend),
            tstart_(_tstart), tend_(_tend),
            cov_bases_q_(0), cov_bases_t_(0),
            num_seeds_(0), score_(-1),
            path_id_(-1), num_paths_(-1),
            segment_id_(-1), num_segments_(-1) {

}

RegionMapped::RegionMapped(
        int32_t _id,                // ID of this current anchor.
        int32_t _target_hits_id,    // ID of the TargetHits object which holds all seeds used to construct this anchor.
        std::shared_ptr<raptor::MappingEnv> _env,
        int32_t _qstart, int32_t _qend,
        int32_t _tstart, int32_t _tend,
        int32_t _cov_bases_q, int32_t _cov_bases_t,
        int32_t _num_seeds, int32_t _score,
        int32_t _path_id, int32_t _num_paths,
        int32_t _segment_id, int32_t _num_segments)
        :   id_(_id),
            target_hits_id_(_target_hits_id),
            env_(_env),
            qstart_(_qstart), qend_(_qend),
            tstart_(_tstart), tend_(_tend),
            cov_bases_q_(_cov_bases_q), cov_bases_t_(_cov_bases_t),
            num_seeds_(_num_seeds), score_(_score),
            path_id_(_path_id), num_paths_(_num_paths),
            segment_id_(_segment_id), num_segments_(_num_segments)  {

}

std::string RegionMapped::Verbose() const {
    std::ostringstream ss;
    ss << "(RegionMapped) "
        << "q = [id:'" << QueryID() << "', s:" << QueryStart() << ", e:" << QueryEnd() << ", len:" << QueryLen() << ", cov:" << CoveredBasesQuery() << "], "
        << "t = [id:'" << TargetID() << "', s:" << TargetStart() << ", e:" << TargetEnd() << ", len:" << TargetLen() << ", cov:" << CoveredBasesTarget()
        << ", rev:" << TargetRev()
        // << ", ind_s:" << TargetIndexStart()
        << "]"
        << ", seeds = " << NumSeeds()
        << ", score = " << Score()
        << ", diag = " << (TargetStart() - QueryStart())
        << ", id = " << id() << ", path_id = " << path_id_ << "/" << num_paths_ << ", seg_id = " << segment_id_ << "/" << num_segments_;

    return ss.str();
}

std::string RegionMapped::WriteAsCSV(const char separator) const {
    std::ostringstream ss;
    ss  << QueryID() << separator
        << TargetID() <<  separator
        << Score() <<  separator
        << NumSeeds() <<  separator
        << 0 <<  separator
        << QueryStart() <<  separator
        << QueryEnd() <<  separator
        << QueryLen() <<  separator
        << ((TargetRev()) ? "1" : "0") <<  separator
        << TargetStart() <<  separator
        << TargetEnd() <<  separator
        << TargetLen() << separator
        << CoveredBasesQuery() << separator
        << CoveredBasesTarget() << separator
        << NumSeeds();
    return ss.str();
}

}
