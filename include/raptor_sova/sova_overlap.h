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
    float score;
    float identity;
    bool a_rev;
    coord_t a_start;
    coord_t a_end;
    coord_t a_len;
    bool b_rev;
    coord_t b_start;
    coord_t b_end;
    coord_t b_len;
    std::string type;

public:
    friend std::unique_ptr<raptor::sova::Overlap> createOverlap();
    friend std::unique_ptr<raptor::sova::Overlap> createOverlap(
                const std::string& a_name, const std::string& b_name,
                float score, float identity,
                bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
                bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
                const std::string& type);
    ~Overlap() = default;

    Overlap()
        :
            a_name(0), b_name(0), score(0), identity(0),
            a_rev(0), a_start(0), a_end(0), a_len(0),
            b_rev(0), b_start(0), b_end(0), b_len(0),
            type()
    { }

    Overlap(const std::string& _a_name, const std::string& _b_name,
            float _score, float _identity,
            bool _a_rev, coord_t _a_start, coord_t _a_end, coord_t _a_len,
            bool _b_rev, coord_t _b_start, coord_t _b_end, coord_t _b_len,
            const std::string& _type)
        :
            a_name(_a_name), b_name(_b_name), score(_score), identity(_identity),
            a_rev(_a_rev), a_start(_a_start), a_end(_a_end), a_len(_a_len),
            b_rev(_b_rev), b_start(_b_start), b_end(_b_end), b_len(_b_len),
            type(_type)
    { }
};

inline std::unique_ptr<raptor::sova::Overlap> createOverlap() {
    return std::unique_ptr<raptor::sova::Overlap>(new raptor::sova::Overlap());
}

inline std::unique_ptr<raptor::sova::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
            const std::string& type) {
    return std::unique_ptr<raptor::sova::Overlap>(new raptor::sova::Overlap(
                                            a_name, b_name, score, identity,
                                            a_rev, a_start, a_end, a_len,
                                            b_rev, b_start, b_end, b_len, type));
}


}
} /* namespace raptor */

#endif /* SRC_RAPTOR_MAPPER_H_ */
