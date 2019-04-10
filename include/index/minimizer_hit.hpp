/*
 * minimizer_hit.hpp
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX2_MINIMIZER_HIT_H_
#define SRC_MINIMIZER_INDEX2_MINIMIZER_HIT_H_

#include <cstdint>
#include <string>
#include <sstream>
#include <index/minimizer_index_types.h>

namespace mindex {

#define MINIMIZER_HIT_TANDEM_FLAG   (0x01)
constexpr mindex128_t MINIMIZER_HIT_TARGET_REV = ((mindex128_t) (0x01)) << 96;

class MinimizerHitPacked;
class MinimizerHitUnpacked;

class MinimizerHitPacked {
   public:
    MinimizerHitPacked()
            : data_(0x0) {
    }

    MinimizerHitPacked(indid_t _t_id, bool _t_rev, ind_t _t_pos,
                        ind_t _q_mask, ind_t _q_pos)
            : data_(PackTo128t(_t_id, _t_rev, _t_pos, _q_mask, _q_pos)) {
    }

    MinimizerHitPacked(mindex128_t _data)
        : data_(_data) {
    }

    MinimizerHitPacked(const MinimizerHitPacked& op)
        : data_(op.data()) {
    }

    bool operator==(const MinimizerHitPacked& op) const {
        return this->data_ == op.data_;
    }

    std::string Verbose() const {
        std::stringstream ss;
        ss << "t_id = " << TargetId() << ", t_rev = " << (TargetRev() ? "true" : "false") << ", t_pos = " << TargetPos() << ", q_mask = " << QueryMask()
           << ", q_pos = " << QueryPos();
        return ss.str();
    }

    std::string WriteAsCSV(char separator) const {
        std::stringstream ss;
        ss << TargetId() << separator << (TargetRev() ? 1 : 0) << separator << TargetPos()
            << separator << QueryMask() << separator << QueryPos();
        return ss.str();
    }

    bool operator<(const mindex::MinimizerHitPacked& op) const {
        return this->data_ < op.data_;
    }

    bool operator>(const mindex::MinimizerHitPacked& op) const {
        return this->data_ > op.data_;
    }

    inline indid_t TargetId() const {
        return DecodeTargetId(data_);
    }

    inline indid_t TargetRev() const {
        return DecodeTargetRev(data_);
    }

    inline ind_t QueryMask() const {
        return DecodeQueryMask(data_);
    }

    inline ind_t TargetPos() const {
        return DecodeTargetPos(data_);
    }

    inline ind_t QueryPos() const {
        return DecodeQueryPos(data_);
    }

    // inline uint64_t TargetQueryPos() const { return DecodeTargetQueryPos(data_); }

    inline mindex128_t data() const {
        return data_;
    }

    static inline mindex128_t PackTo128t(indid_t _t_id, bool _t_rev, ind_t _t_pos, ind_t _q_mask, ind_t _q_pos) {
        mindex128_t ret = ((mindex128_t)_t_id) << 97;
        if (_t_rev) {
            ret |= MINIMIZER_HIT_TARGET_REV;
        }
        ret |= (((mindex128_t)_t_pos) & MINIMIZER_32bit_MASK_FULL) << 64;
        ret |= (((mindex128_t)_q_mask) & MINIMIZER_32bit_MASK_FULL) << 32;
        ret |= (((mindex128_t)_q_pos) & MINIMIZER_32bit_MASK_FULL);
        return ret;
    }

    static inline indid_t DecodeTargetId(mindex128_t packed) {
        // The bit 96 encodes the reverse complement flag.
        return (indid_t)(packed >> 97);
    }

    static inline bool DecodeTargetRev(mindex128_t packed) {
        return static_cast<bool>(packed & MINIMIZER_HIT_TARGET_REV);
    }

    static inline ind_t DecodeTargetPos(mindex128_t packed) {
        return (ind_t)((packed >> 64) & MINIMIZER_32bit_MASK_FULL);
    }

    static inline ind_t DecodeQueryMask(mindex128_t packed) {
        return (ind_t)(((packed >> 32) & MINIMIZER_32bit_MASK_FULL));
    }

    static inline ind_t DecodeQueryPos(mindex128_t packed) {
        return (ind_t)(packed & MINIMIZER_32bit_MASK_FULL);
    }

    // static inline uint64_t DecodeTargetQueryPos(mindex128_t packed) {
    //     return ((uint64_t)(packed & MINIMIZER_64bit_MASK));
    // }

   private:
    mindex128_t data_;
};

class MinimizerHitUnpacked {
   public:
    MinimizerHitUnpacked()
        : t_id_(0), t_rev_(false), t_pos_(0),
          q_mask_(0), q_pos_(0) {
    }

    MinimizerHitUnpacked(indid_t _t_id, bool _t_rev, ind_t _t_pos, ind_t _q_pos)
        : t_id_(_t_id),
          t_rev_(_t_rev),
          t_pos_(_t_pos),
          q_mask_(0),  // Query mask.
          q_pos_(_q_pos) {
    }

    MinimizerHitUnpacked(const mindex128_t& packed)
        : t_id_(MinimizerHitPacked::DecodeTargetId(packed)),
          t_rev_(MinimizerHitPacked::DecodeTargetRev(packed)),
          t_pos_(MinimizerHitPacked::DecodeTargetPos(packed)),
          q_mask_(MinimizerHitPacked::DecodeQueryMask(packed)),
          q_pos_(MinimizerHitPacked::DecodeQueryPos(packed)) {
    }

    MinimizerHitUnpacked(const MinimizerHitPacked& packed)
        : t_id_(packed.TargetId()),
          t_rev_(packed.TargetRev()),
          t_pos_(packed.TargetPos()),
          q_mask_(packed.QueryMask()),
          q_pos_(packed.QueryPos()) {
    }

    ~MinimizerHitUnpacked() = default;

    inline mindex128_t to_128t() {
        return MinimizerHitPacked::PackTo128t(t_id_, t_rev_, t_pos_, q_mask_, q_pos_);
    }

    inline indid_t TargetId() const {
        return t_id_;
    }

    inline bool TargetRev() const {
        return t_rev_;
    }

    inline ind_t QueryMask() const {
        return q_mask_;
    }

    inline ind_t TargetPos() const {
        return t_pos_;
    }

    inline ind_t QueryPos() const {
        return q_pos_;
    }

  private:
    indid_t t_id_;
    bool t_rev_;
    ind_t t_pos_;
    ind_t q_mask_;
    ind_t q_pos_;
};

inline std::vector<MinimizerHitUnpacked> UnpackMinimizerHitVector(const std::vector<MinimizerHitPacked>& packed) {
    std::vector<MinimizerHitUnpacked> unpacked;
    unpacked.reserve(packed.size());
    for (auto& hit: packed) {
        unpacked.emplace_back(MinimizerHitUnpacked(hit));
    }
    return unpacked;
}

}  // namespace mindex

#endif
