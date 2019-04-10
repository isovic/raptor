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

class IndexParams;

std::shared_ptr<mindex::IndexParams> createIndexParams();
std::shared_ptr<mindex::IndexParams> createIndexParams(int32_t _k, int32_t _w, bool _hp_supp,
                                                       int32_t _max_hp_len, double _freq_percentile,
                                                       int32_t _min_occ_cutoff,
                                                       bool _index_only_fwd_strand,
                                                       int64_t min_tlen);

class IndexParams {
   public:
    friend std::shared_ptr<mindex::IndexParams> createIndexParams();
    friend std::shared_ptr<mindex::IndexParams> createIndexParams(
        int32_t _k, int32_t _w, bool _hp_supp, int32_t _max_hp_len, double _freq_percentile,
        int32_t _min_occ_cutoff, bool _index_only_fwd_strand, int64_t _min_tlen);
    ~IndexParams();

    int32_t k;
    int32_t w;
    bool homopolymer_suppression;
    int32_t max_homopolymer_len;
    double freq_percentil;
    int32_t min_occ_cutoff;

    bool is_region_specified;
    std::string region_rname;
    int64_t region_rstart;
    int64_t region_rend;

    bool index_only_fwd_strand;

    int64_t min_tlen;

   private:
    IndexParams();
    IndexParams(int32_t _k, int32_t _w, bool _hp_supp, int32_t _max_hp_len, double _freq_percentile,
                int32_t _min_occ_cutoff, bool _index_only_fwd_strand, int64_t _min_tlen);
};

}  // namespace mindex

#endif
