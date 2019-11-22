/*
 * seed_hit_diag.hpp
 *
 *  Created on: Nov 14, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX2_SEED_HIT_DIAGONAL_H_
#define SRC_MINIMIZER_INDEX2_SEED_HIT_DIAGONAL_H_

#include <cstdint>
#include <string>
#include <sstream>
#include <vector>
#include <index/index_types.h>
#include <index/seed_hit.hpp>

namespace mindex {

#define MINIMIZER_HIT_TANDEM_FLAG   (0x01)
constexpr mindex128_t SEED_HIT_DIAG_TARGET_REV_BIT = ((mindex128_t) (0x01)) << 96;
const mindex128_t MASK_SEED_HIT_DIAG_0 = (static_cast<mindex128_t>(0x0FFFFFFFF));
const mindex128_t MASK_SEED_HIT_DIAG_1 = (static_cast<mindex128_t>(0x0FFFFFFFF) << 32);
const mindex128_t MASK_SEED_HIT_DIAG_2 = (static_cast<mindex128_t>(0x0FFFFFFFF) << 64);
const mindex128_t MASK_SEED_HIT_DIAG_3 = (static_cast<mindex128_t>(0x0FFFFFFFF) << 96);
const mindex128_t MASK_SEED_HIT_DIAG_TARGET_ID = (static_cast<mindex128_t>(0x0FFFFFFFF) << 97);
const mindex128_t MASK_INV_SEED_HIT_DIAG_0 = MASK_SEED_HIT_DIAG_3 | MASK_SEED_HIT_DIAG_2 | MASK_SEED_HIT_DIAG_1;
const mindex128_t MASK_INV_SEED_HIT_DIAG_1 = MASK_SEED_HIT_DIAG_3 | MASK_SEED_HIT_DIAG_2 | MASK_SEED_HIT_DIAG_0;
const mindex128_t MASK_INV_SEED_HIT_DIAG_2 = MASK_SEED_HIT_DIAG_3 | MASK_SEED_HIT_DIAG_1 | MASK_SEED_HIT_DIAG_0;
const mindex128_t MASK_INV_SEED_HIT_DIAG_3 = MASK_SEED_HIT_DIAG_2 | MASK_SEED_HIT_DIAG_1 | MASK_SEED_HIT_DIAG_0;
const mindex128_t MASK_INV_SEED_HIT_TARGET_ID = SEED_HIT_DIAG_TARGET_REV_BIT | MASK_SEED_HIT_DIAG_2 | MASK_SEED_HIT_DIAG_1 | MASK_SEED_HIT_DIAG_0;
const mindex128_t MASK_INV_SEED_HIT_TARGET_REV = MASK_SEED_HIT_DIAG_TARGET_ID | MASK_SEED_HIT_DIAG_2 | MASK_SEED_HIT_DIAG_1 | MASK_SEED_HIT_DIAG_0;

const mindex128_t MASK_SEED_HIT_DIAG_LOWER_64 = MASK_SEED_HIT_DIAG_0 | MASK_SEED_HIT_DIAG_1;

class SeedHitDiagPacked;
class SeedHitDiagUnpacked;

class SeedHitDiagPacked {
   public:
    SeedHitDiagPacked()
            : data_(0x0) {
    }

    SeedHitDiagPacked(indid_t _t_id, bool _t_rev, ind_t _t_pos,
                        ind_t _q_pos)
            : data_(PackTo128t(_t_id, _t_rev, _t_pos, _q_pos)) {
    }

    SeedHitDiagPacked(mindex128_t _data)
        : data_(_data) {
    }

    SeedHitDiagPacked(const SeedHitPacked& shp)
        : data_(RepackFromSeedHit(shp)) {
    }

    SeedHitDiagPacked(const SeedHitDiagPacked& op)
        : data_(op.data()) {
    }

    bool operator==(const SeedHitDiagPacked& op) const {
        return this->data_ == op.data_;
    }

    std::string Verbose() const {
        std::stringstream ss;
        ss << "t_id = " << TargetId() << ", t_rev = " << (TargetRev() ? "true" : "false")
            << ", diag = " << Diagonal()
            << ", t_pos = " << TargetPos() << ", q_pos = " << QueryPos();
        return ss.str();
    }

    std::string WriteAsCSV(char separator) const {
        std::stringstream ss;
        ss << TargetId() << separator << (TargetRev() ? 1 : 0)
            << separator << Diagonal()
            << separator << TargetPos() << separator << QueryPos();
        return ss.str();
    }

    bool operator<(const mindex::SeedHitDiagPacked& op) const {
        return this->data_ < op.data_;
    }

    bool operator>(const mindex::SeedHitDiagPacked& op) const {
        return this->data_ > op.data_;
    }

    inline indid_t TargetId() const {
        return DecodeTargetId(data_);
    }

    inline indid_t TargetRev() const {
        return DecodeTargetRev(data_);
    }

    inline ind_t Diagonal() const {
        return DecodeDiagonal(data_);
    }

    inline ind_t TargetPos() const {
        return DecodeTargetPos(data_);
    }

    inline ind_t QueryPos() const {
        return DecodeQueryPos(data_);
    }

    inline indid_t TargetIdRev() const {
        return DecodeTargetIdRev(data_);
    }

    inline mindex128_t data() const {
        return data_;
    }

    inline void TargetId(int32_t val) {
        data_ &= MASK_INV_SEED_HIT_TARGET_ID;
        data_ |= ((mindex128_t)val) << 97;
    }

    inline void TargetRev(bool val) {
        data_ &= MASK_INV_SEED_HIT_TARGET_REV;
        data_ |= ((mindex128_t)val) << 96;
    }

    inline void Diagonal(int32_t val) {
        data_ &= MASK_INV_SEED_HIT_DIAG_2;
        data_ |= (((mindex128_t)val) & MINIMIZER_32bit_MASK_FULL) << 64;
    }

    inline void TargetPos(int32_t val) {
        data_ &= MASK_INV_SEED_HIT_DIAG_1;
        data_ |= ((mindex128_t)val) << 32;
    }

    inline void QueryPos(int32_t val) {
        data_ &= MASK_INV_SEED_HIT_DIAG_0;
        data_ |= ((mindex128_t)val);
    }

    static inline mindex128_t PackTo128t(indid_t _t_id, bool _t_rev, ind_t _t_pos, ind_t _q_pos) {
        mindex128_t ret = ((mindex128_t)_t_id) << 97;
        if (_t_rev) {
            ret |= SEED_HIT_DIAG_TARGET_REV_BIT;
        }
        int32_t diag = _t_pos - _q_pos;
        ret |= (((mindex128_t)diag) & MINIMIZER_32bit_MASK_FULL) << 64;
        ret |= (((mindex128_t)_t_pos) & MINIMIZER_32bit_MASK_FULL) << 32;
        ret |= (((mindex128_t)_q_pos) & MINIMIZER_32bit_MASK_FULL);
        return ret;
    }

    static inline indid_t DecodeTargetId(mindex128_t packed) {
        // The bit 96 encodes the reverse complement flag.
        return (indid_t)(packed >> 97);
    }

    static inline bool DecodeTargetRev(mindex128_t packed) {
        return static_cast<bool>(packed & SEED_HIT_DIAG_TARGET_REV_BIT);
    }

    static inline ind_t DecodeDiagonal(mindex128_t packed) {
        return (ind_t)(((packed >> 64) & MINIMIZER_32bit_MASK_FULL));
    }

    static inline ind_t DecodeTargetPos(mindex128_t packed) {
        return (ind_t)((packed >> 32) & MINIMIZER_32bit_MASK_FULL);
    }

    static inline ind_t DecodeQueryPos(mindex128_t packed) {
        return (ind_t)(packed & MINIMIZER_32bit_MASK_FULL);
    }

    static inline indid_t DecodeTargetIdRev(mindex128_t packed) {
        // The bit 96 encodes the reverse complement flag.
        return (indid_t)(packed >> 96);
    }

    static inline mindex128_t RepackFromSeedHit(const SeedHitPacked& shp) {
        SeedHitDiagPacked shdp(shp.TargetId(), shp.TargetRev(),
                                shp.TargetPos(), shp.QueryPos());
        return shdp.data();
    }

    static inline void EncodeDiagonal(mindex128_t& packed, int32_t val) {
        packed &= MASK_INV_SEED_HIT_DIAG_2;
        packed |= (((mindex128_t)val) & MINIMIZER_32bit_MASK_FULL) << 64;
    }


   private:
    mindex128_t data_;
};

class SeedHitDiagUnpacked {
   public:
    SeedHitDiagUnpacked()
        : t_id_(0), t_rev_(false),
          diag_(0), t_pos_(0), q_pos_(0) {
    }

    SeedHitDiagUnpacked(indid_t _t_id, bool _t_rev, ind_t _t_pos, ind_t _q_pos)
        : t_id_(_t_id),
          t_rev_(_t_rev),
          diag_(_t_pos - _q_pos),
          t_pos_(_t_pos),
          q_pos_(_q_pos) {
    }

    SeedHitDiagUnpacked(indid_t _t_id, bool _t_rev, ind_t _diag, ind_t _t_pos,
                        ind_t _q_pos)
        : t_id_(_t_id),
          t_rev_(_t_rev),
          diag_(_diag),
          t_pos_(_t_pos),
          q_pos_(_q_pos) {
    }


    SeedHitDiagUnpacked(const mindex128_t& packed)
        : t_id_(SeedHitDiagPacked::DecodeTargetId(packed)),
          t_rev_(SeedHitDiagPacked::DecodeTargetRev(packed)),
          diag_(SeedHitDiagPacked::DecodeDiagonal(packed)),
          t_pos_(SeedHitDiagPacked::DecodeTargetPos(packed)),
          q_pos_(SeedHitDiagPacked::DecodeQueryPos(packed)) {
    }

    SeedHitDiagUnpacked(const SeedHitDiagPacked& packed)
        : t_id_(packed.TargetId()),
          t_rev_(packed.TargetRev()),
          diag_(packed.Diagonal()),
          t_pos_(packed.TargetPos()),
          q_pos_(packed.QueryPos()) {
    }

    ~SeedHitDiagUnpacked() = default;

    inline mindex128_t to_128t() {
        return SeedHitDiagPacked::PackTo128t(t_id_, t_rev_, t_pos_, q_pos_);
    }

    inline indid_t TargetId() const {
        return t_id_;
    }

    inline ind_t Diagonal() const {
        return diag_;
    }

    inline bool TargetRev() const {
        return t_rev_;
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
    ind_t diag_;
    ind_t t_pos_;
    ind_t q_pos_;
};

inline std::vector<SeedHitDiagUnpacked> UnpacSeedHitDiagVector(const std::vector<SeedHitDiagPacked>& packed) {
    std::vector<SeedHitDiagUnpacked> unpacked;
    unpacked.reserve(packed.size());
    for (auto& hit: packed) {
        unpacked.emplace_back(SeedHitDiagUnpacked(hit));
    }
    return unpacked;
}

inline std::vector<SeedHitDiagPacked> ConvertToSeedHitDiagPackedVector(const std::vector<SeedHitPacked>& in_shp) {
    std::vector<SeedHitDiagPacked> ret;
    ret.reserve(in_shp.size());
    for (auto& hit: in_shp) {
        ret.emplace_back(SeedHitDiagPacked(hit));
    }
    return ret;
}

inline std::vector<mindex128_t> ConvertToSeedHitDiagPacked128Vector(const std::vector<SeedHitPacked>& in_shp) {
    std::vector<mindex128_t> ret;
    ret.reserve(in_shp.size());
    for (auto& hit: in_shp) {
        ret.emplace_back(SeedHitDiagPacked::RepackFromSeedHit(hit));
    }
    return ret;
}

}  // namespace mindex

#endif
