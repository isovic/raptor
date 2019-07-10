#ifndef _IPA_INCLUDE_OVERLAP_FILTER_H_
#define _IPA_INCLUDE_OVERLAP_FILTER_H_

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <raptor_fetch/overlaps/overlap_file.h>
#include <raptor_fetch/overlaps/overlap_compact.h>

namespace raptor {

class QueryRegion {
  public:
	QueryRegion() : id(0), start(0), end(0), len(0) {}
	QueryRegion(int32_t _id, int32_t _start, int32_t _end, int32_t _len) : id(_id), start(_start), end(_end), len(_len) {}
	int32_t id;
	int32_t start;
	int32_t end;
	int32_t len;
};

class CoverageEvent {
  public:
	CoverageEvent() : pos(0), is_start(0), source_id(0), cov(0), fwd_closing(0) {}
	CoverageEvent(int32_t _pos, bool _is_start, int32_t _source_id, int32_t _cov, int32_t _fwd_closing)
		: pos(_pos), is_start(_is_start), source_id(_source_id), cov(_cov), fwd_closing(_fwd_closing) {}
	int32_t pos;
	bool is_start;
	int32_t source_id;
	int32_t cov;
	int32_t fwd_closing;
};

enum class OverlapType {
	Overlap,		// Generic overlap, perhaps not valid, but can't be called internal or containment clearly.
	Dovetail5,		// Nice and clean dovetail overlap on the 5' end of the A-read.
	Dovetail3,		// Nice and clean dovetail overlap on the 3' end of the A-read.
	Internal,
	ContainedA,
	ContainedB,
	Sketchy			// It's not completely clear what type of overlap this is, so best not to filter it out.
};

inline std::string OverlapTypeToString(OverlapType ot) {
	std::string ret = "unknown";
	switch(ot) {
		case OverlapType::Overlap:
			ret = "overlap";
			break;
		case OverlapType::Dovetail5:
			ret = "dovetail5";
			break;
		case OverlapType::Dovetail3:
			ret = "dovetail3";
			break;
		case OverlapType::Internal:
			ret = "internal";
			break;
		case OverlapType::ContainedA:
			ret = "contained_A";
			break;
		case OverlapType::ContainedB:
			ret = "contained_B";
			break;
		case OverlapType::Sketchy:
			ret = "sketchy";
			break;
		default:
			ret = "unknown";
			break;
	}
	return ret;
}

raptor::OverlapCompactPtr GetDovetailOverlap(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 int32_t max_allowed_clip, int32_t min_edge_len);



OverlapType DetermineOverlapTypeNew(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 // If min_lhang or min_hang is > max_valid_aln_clip, it's called an internal overlap.
								 int32_t internal_max_valid_aln_clip,
								 // Used for checking internal overlaps.
								 float internal_min_frac_cov,
								 // Perform coordinate offsetting only if the distance is <= offset_max_dist. Offsetting is performed before calling contained and dovetail overlaps.
								 int32_t offset_max_dist,
                                 // Similar to internal_max_valid_aln_clip, but should be smaller. If unaligned portion on the shorter end is <= internal_max_valid_aln_clip,
                                 // this value will be used to offset the coordinates and create fake dovetail overlaps.
                                 // The fake dovetail overlap will be used to call contained reads or dovetail overlaps.
								 int32_t dovetail_max_valid_aln_clip,
 								 // For containment. If e.g. A is end-to-end contained within B, it still has to be this far from one of the edges of B, just in case.
								 int32_t contained_min_dist_from_edge
								 );



raptor::OverlapPtr CalculateClippedOverlap(const raptor::OverlapCompactPtr& ovl,
								 			const raptor::QueryDataPtr& qda,
								 			const raptor::QueryDataPtr& qdb);

bool FindClippingRegions(raptor::OverlapFilePtr& ovl_file, size_t batch_start, size_t bach_end,
                         int32_t fwd_dist,
						 int32_t min_user_cov_threshold, bool clip_ends_on_coverage,
						 std::vector<QueryRegion>& ret_regions);

int32_t FilterRegions(const std::vector<QueryRegion>& in_regions, int32_t min_span, int32_t bestn, std::vector<QueryRegion>& out_regions);


OverlapType DetermineOverlapType(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 int32_t max_valid_aln_clip,							// Used for calling contained reads, should be smaller tan max_internal_hang, but not very high to not call SV regions falsely.
								 int32_t contained_min_dist_from_edge,						// For containment. If e.g. A is end-to-end contained within B, it still has to be this far from one of the edges of B, just in case.
								 int32_t internal_min_hang, float min_frac_cov		// Used for checking internal overlaps.
								 );

std::string FormatOverlapForFalcon(const raptor::OverlapCompactPtr& ovl, const raptor::QueryDataPtr& qda, const raptor::QueryDataPtr& qdb);
std::string FormatOverlapToCage11(const raptor::OverlapCompactPtr& ovl, const raptor::QueryDataPtr& qda, const raptor::QueryDataPtr& qdb);

void FindQueryClips(int32_t fwd_dist, int32_t min_cov, int32_t min_len, raptor::OverlapFilePtr& overlap_file);
void ApplyClips(raptor::OverlapFilePtr& overlap_file);

}

#endif
