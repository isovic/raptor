#include <index/index_params.h>

namespace mindex {

std::shared_ptr<mindex::IndexParams> createIndexParams() {
    return std::shared_ptr<mindex::IndexParams>(new mindex::IndexParams());
}

std::shared_ptr<mindex::IndexParams> createIndexParams(int32_t _k, int32_t _w, bool _hp_supp,
                                                       int32_t _max_hp_len, double _freq_percentile,
                                                       int32_t _min_occ_cutoff,
                                                       bool _index_only_fwd_strand,
                                                       int64_t _min_tlen) {
    return std::shared_ptr<mindex::IndexParams>(new mindex::IndexParams(
        _k, _w, _hp_supp, _max_hp_len, _freq_percentile, _min_occ_cutoff, _index_only_fwd_strand, _min_tlen));
}

IndexParams::IndexParams()
    : k(15),
      w(5),
      homopolymer_suppression(false),
      max_homopolymer_len(10),
      freq_percentil(0.002),
      min_occ_cutoff(-1),
      is_region_specified(false),
      region_rname(),
      region_rstart(0),
      region_rend(-1),
      index_only_fwd_strand(false),
      min_tlen(0) {}

IndexParams::IndexParams(int32_t _k, int32_t _w, bool _hp_supp, int32_t _max_hp_len,
                         double _freq_percentile, int32_t _min_occ_cutoff,
                         bool _index_only_fwd_strand, int64_t _min_tlen)
    : k(_k),
      w(_w),
      homopolymer_suppression(_hp_supp),
      max_homopolymer_len(_max_hp_len),
      freq_percentil(_freq_percentile),
      min_occ_cutoff(_min_occ_cutoff),
      is_region_specified(false),
      region_rname(),
      region_rstart(0),
      region_rend(-1),
      index_only_fwd_strand(_index_only_fwd_strand),
      min_tlen(_min_tlen) {}

IndexParams::~IndexParams() {}

}  // namespace mindex
