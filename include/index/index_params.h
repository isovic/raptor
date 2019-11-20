/*
 * index_params.h
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX_2_INDEX_PARAMS
#define SRC_MINIMIZER_INDEX_2_INDEX_PARAMS

#include <cstdint>
#include <vector>
#include <string>
#include <memory>

namespace mindex {

class IndexParams {
   public:
    friend std::shared_ptr<mindex::IndexParams> createIndexParams();
    friend std::shared_ptr<mindex::IndexParams> createIndexParams(
        int32_t _k, int32_t _w, bool _hp_supp, int32_t _max_hp_len, double _freq_percentile,
        int32_t _min_occ_cutoff, bool _index_only_fwd_strand, int64_t _min_tlen);
    ~IndexParams() = default;

    int32_t k = 15;
    int32_t w = 5;
    bool homopolymer_suppression = false;
    int32_t max_homopolymer_len = 10;
    double freq_percentil = 0.002;
    int32_t min_occ_cutoff = -1;
    bool is_region_specified = false;
    std::string region_rname;
    int64_t region_rstart = 0;
    int64_t region_rend = -1;
    bool index_only_fwd_strand = false;
    int64_t min_tlen = 0;

   private:
    IndexParams() = default;
    IndexParams(int32_t _k, int32_t _w, bool _hp_supp, int32_t _max_hp_len, double _freq_percentile,
                int32_t _min_occ_cutoff, bool _index_only_fwd_strand, int64_t _min_tlen);
    IndexParams(const IndexParams&) = delete;
    IndexParams& operator=(const IndexParams&) = delete;
};

inline std::shared_ptr<mindex::IndexParams> createIndexParams() {
    return std::shared_ptr<mindex::IndexParams>(new mindex::IndexParams());
}

inline std::shared_ptr<mindex::IndexParams> createIndexParams(int32_t _k, int32_t _w, bool _hp_supp,
                                                       int32_t _max_hp_len, double _freq_percentile,
                                                       int32_t _min_occ_cutoff,
                                                       bool _index_only_fwd_strand,
                                                       int64_t _min_tlen) {
    return std::shared_ptr<mindex::IndexParams>(new mindex::IndexParams(
        _k, _w, _hp_supp, _max_hp_len, _freq_percentile, _min_occ_cutoff, _index_only_fwd_strand, _min_tlen));
}


}  // namespace mindex

#endif
