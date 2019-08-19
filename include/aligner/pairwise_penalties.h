/*
 * pairwise_penalties.h
 *
 *  Created on: Jul 1, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_ALIGNER_PAIRWISE_PENALTIES_H_
#define SRC_ALIGNER_PAIRWISE_PENALTIES_H_

#include <cstdint>
#include <vector>
#include <string>
#include <sstream>
#include <vector>

namespace raptor {

/* Regular alignment penalties for a single piece Gotoh alignment.
 */
class Penalties {
   public:
    Penalties() : match(5), mismatch(-4), gapopen(-8), gapext(-6) {}
    Penalties(int32_t _match, int32_t _mismatch, int32_t _gapopen, int32_t _gapext)
        : match(_match), mismatch(_mismatch), gapopen(_gapopen), gapext(_gapext) {}
    int32_t match, mismatch, gapopen, gapext;
};

/* A helper class for a linear function. Used for piecewise Gotoh alignment.
 */
class AffinePiece {
   public:
    AffinePiece() : ext(-6.0), open(-8.0) {}
    AffinePiece(float _ext, float _open) : ext(_ext), open(_open) {}

    inline float calc(int32_t k) const { return (ext * (k) + open); }

    float ext, open;  // Line equation parameters: w(k) = ext * k + open.
};

/* Penalties for a multiple affine function alignment.
 */
class PiecewisePenalties {
   public:
    PiecewisePenalties()
        : match(5), mismatch(-4), w(std::vector<AffinePiece>{AffinePiece(-6.0, -8.0)}) {}

    PiecewisePenalties(int32_t _match, int32_t _mismatch, const std::vector<AffinePiece>& _w)
        : match(_match), mismatch(_mismatch), w(_w) {}

    std::string Verbose() {
        std::stringstream ss;
        ss << "match = " << match << ", mismatch = " << mismatch << "";
        for (size_t l = 0; l < w.size(); l++) {
            ss << ", w[" << l << "] = {ext = " << w[l].ext << ", open = " << w[l].open << "}";
        }
        ss << "\n";
        return ss.str();
    }

    float match, mismatch;
    std::vector<AffinePiece> w;
};

}  // namespace raptor

#endif