/*
 * minimizer.hpp
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX2_MINIMIZER_H_
#define SRC_MINIMIZER_INDEX2_MINIMIZER_H_

#include <index/minimizer_index_types.h>
#include <sstream>

namespace mindex {

static const int8_t MINIMIZER_FLAG_DEFAULT_FWD = 0x00;
static const int8_t MINIMIZER_FLAG_IS_REV = 0x01;

static const mindex128_t MINIMIZER_CODED_REV_BIT = (((mindex128_t)1) << 32);

class Minimizer {
   public:
    Minimizer() : key(0), seq_id(0), pos(0), flag(0) {}

    Minimizer(minkey_t _key, indid_t _seq_id, ind_t _pos, bool _is_rev)
        : key(_key),
          seq_id(_seq_id),
          pos(_pos),
          flag((_is_rev) ? MINIMIZER_FLAG_IS_REV : MINIMIZER_FLAG_DEFAULT_FWD) {}
    Minimizer(const mindex128_t& coded_keypos)
        : key(DecodeKey(coded_keypos)),
          seq_id(DecodeSeqId(coded_keypos)),
          pos(DecodePos(coded_keypos)),
          flag(DecodeIsRev(coded_keypos)) {}

    bool is_rev() { return (flag & MINIMIZER_FLAG_IS_REV); }

    inline mindex128_t to_128t() {
        mindex128_t ret = ((mindex128_t)key) << 64;
        ret |= ((((mindex128_t)(seq_id)) << 1) & MINIMIZER_32bit_MASK) << 32;
        ret |= ((mindex128_t)(pos)) & MINIMIZER_32bit_MASK;

        if (flag & MINIMIZER_FLAG_IS_REV) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }

        return ret;
    }

    static inline mindex128_t Encode(minkey_t _key, indid_t _seq_id, ind_t _pos, bool _is_rev) {
        mindex128_t ret = ((mindex128_t)_key) << 64;
        ret |= (((mindex128_t)(_seq_id & MINIMIZER_32bit_MASK)) << 33);
        ret |= ((mindex128_t)(_pos & MINIMIZER_32bit_MASK));
        if (_is_rev) {
            ret |= MINIMIZER_CODED_REV_BIT;
        }
        return ret;
    }

    static inline minkey_t DecodeKey(const mindex128_t& seed) {
        return (minkey_t)((seed >> 64) & MINIMIZER_64bit_MASK);
    }

    static inline ind_t DecodePos(const mindex128_t& seed) { return ((ind_t)(seed & MINIMIZER_32bit_MASK)); }

    static inline indid_t DecodeSeqId(const mindex128_t& seed) {
        return (indid_t)((seed >> 33) & MINIMIZER_32bit_MASK);
    }

    /*
     * Unlike DecodeSeqId, this method keeps the info about
     * of the strand still encoded in the return value.
     * The LSB is 0 for fwd, and 1 for rev.
     */
    static inline indid_t DecodeSeqIdWithRev(const mindex128_t& seed) {
        return (indid_t)((seed >> 32) & MINIMIZER_32bit_MASK);
    }

    static inline bool DecodeIsRev(const mindex128_t& seed) { return (seed & MINIMIZER_CODED_REV_BIT); }

    std::string Verbose() const {
        std::stringstream ss;
        ss << "pos = " << pos << ", seq_id = " << seq_id << ", flag = " << (int32_t)flag
           << ", key = " << key;
        return ss.str();
    }

    minkey_t key;
    indid_t seq_id;
    ind_t pos;
    int8_t flag;
};

}  // namespace mindex

#endif