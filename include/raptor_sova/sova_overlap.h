/*
 * sova_overlap.h
 *
 *  Created on: Nov 21, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_SOVA_OVERLAP_H_
#define SRC_RAPTOR_SOVA_OVERLAP_H_

#include <cstdint>
#include <memory>
#include <string>

namespace raptor {
namespace sova {

class Overlap;

using coord_t = int32_t;
using OverlapPtr = std::unique_ptr<raptor::sova::Overlap>;

class Overlap {
public:
    std::string a_name;
    std::string b_name;
    float score = 0.0f;
    float identity = 0.0f;
    bool a_rev = false;
    coord_t a_start = 0;
    coord_t a_end = 0;
    coord_t a_len = 0;
    bool b_rev = false;
    coord_t b_start = 0;
    coord_t b_end = 0;
    coord_t b_len = 0;
    std::string type;
    int32_t a_id = -1;
    int32_t b_id = -1;
    int32_t edit_dist = -1;
    int32_t num_seeds = 0;

public:
    friend std::unique_ptr<raptor::sova::Overlap> createOverlap();
    friend std::unique_ptr<raptor::sova::Overlap> createOverlap(
                const std::string& a_name, const std::string& b_name,
                float score, float identity,
                bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
                bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
                const std::string& type,
                int32_t a_id, int32_t b_id,
                int32_t edit_dist, int32_t num_seeds);
    friend std::unique_ptr<raptor::sova::Overlap> createOverlap(const std::unique_ptr<raptor::sova::Overlap>& ovl);

    Overlap() = default;
    ~Overlap() = default;

    Overlap(const std::string& _a_name, const std::string& _b_name,
            float _score, float _identity,
            bool _a_rev, coord_t _a_start, coord_t _a_end, coord_t _a_len,
            bool _b_rev, coord_t _b_start, coord_t _b_end, coord_t _b_len,
            const std::string& _type, int32_t _a_id, int32_t _b_id,
            int32_t _edit_dist, int32_t _num_seeds)
        :
            a_name(_a_name), b_name(_b_name), score(_score), identity(_identity),
            a_rev(_a_rev), a_start(_a_start), a_end(_a_end), a_len(_a_len),
            b_rev(_b_rev), b_start(_b_start), b_end(_b_end), b_len(_b_len),
            type(_type), a_id(_a_id), b_id(_b_id),
            edit_dist(_edit_dist), num_seeds(_num_seeds)
    { }

public:
    int32_t ASpan() const { return (a_end - a_start); }
    int32_t BSpan() const { return (b_end - b_start); }

};

inline std::unique_ptr<raptor::sova::Overlap> createOverlap() {
    return std::unique_ptr<raptor::sova::Overlap>(new raptor::sova::Overlap());
}

inline std::unique_ptr<raptor::sova::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
            const std::string& type, int32_t a_id, int32_t b_id,
            int32_t edit_dist, int32_t num_seeds) {
    return std::unique_ptr<raptor::sova::Overlap>(new raptor::sova::Overlap(
                                            a_name, b_name, score, identity,
                                            a_rev, a_start, a_end, a_len,
                                            b_rev, b_start, b_end, b_len, type,
                                            a_id, b_id,
                                            edit_dist, num_seeds));
}

inline std::unique_ptr<raptor::sova::Overlap> createOverlap(
            const std::unique_ptr<raptor::sova::Overlap>& ovl) {
    return std::unique_ptr<raptor::sova::Overlap>(new raptor::sova::Overlap(
                                            ovl->a_name, ovl->b_name, ovl->score, ovl->identity,
                                            ovl->a_rev, ovl->a_start, ovl->a_end, ovl->a_len,
                                            ovl->b_rev, ovl->b_start, ovl->b_end, ovl->b_len, ovl->type,
                                            ovl->a_id, ovl->b_id,
                                            ovl->edit_dist, ovl->num_seeds));
}

}
} /* namespace raptor */

#endif /* SRC_RAPTOR_MAPPER_H_ */
