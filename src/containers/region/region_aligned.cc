#include <containers/region/region_aligned.h>

namespace raptor {

std::shared_ptr<raptor::RegionAligned> createAlignedRegion(std::shared_ptr<raptor::MappingEnv> env,
                                                       std::shared_ptr<raptor::AlignmentResult> aln,
                                                       int32_t path_id, int32_t num_paths,
                                                       int32_t segment_id, int32_t num_segments) {
    return std::shared_ptr<raptor::RegionAligned>(new raptor::RegionAligned(env, aln, path_id, num_paths, segment_id, num_segments));
}

std::string RegionAligned::WriteAsCSV(const char separator) const {
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
