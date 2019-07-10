#ifndef _IPA_INCLUDE_OVERLAP_COMPACT_H_
#define _IPA_INCLUDE_OVERLAP_COMPACT_H_

#include <cstdint>
#include <memory>

namespace raptor {

class OverlapCompact;

using seqid_t = int32_t;
using coord_t = int32_t;

constexpr uint32_t FLAG_OVERLAP_LOW_LEN     = ((static_cast<uint32_t>(1)) << 0);
constexpr uint32_t FLAG_OVERLAP_LOW_SCORE   = ((static_cast<uint32_t>(1)) << 1);
constexpr uint32_t FLAG_OVERLAP_LOW_IDT     = ((static_cast<uint32_t>(1)) << 2);
constexpr uint32_t FLAG_OVERLAP_CONTAINED_A = ((static_cast<uint32_t>(1)) << 3);
constexpr uint32_t FLAG_OVERLAP_CONTAINED_B = ((static_cast<uint32_t>(1)) << 4);
constexpr uint32_t FLAG_OVERLAP_FILTERED_A  = ((static_cast<uint32_t>(1)) << 5);
constexpr uint32_t FLAG_OVERLAP_FILTERED_B  = ((static_cast<uint32_t>(1)) << 6);
constexpr uint32_t FLAG_OVERLAP_LOW_SPAN    = ((static_cast<uint32_t>(1)) << 7);
constexpr uint32_t FLAG_OVERLAP_INTERNAL    = ((static_cast<uint32_t>(1)) << 8);
constexpr uint32_t FLAG_OVERLAP_NOT_DOVETAIL    = ((static_cast<uint32_t>(1)) << 9);

constexpr uint32_t FLAG_OVERLAP_TYPE_5_PRIME    = ((static_cast<uint32_t>(1)) << 0);
constexpr uint32_t FLAG_OVERLAP_TYPE_3_PRIME    = ((static_cast<uint32_t>(1)) << 1);

// constexpr uint32_t FLAG_OVERLAP_IS_CLIPPED_TO_REGION    = ((static_cast<uint32_t>(1)) << 9);    // If true, the overlap coordinates are already offset by the query_data->clip_start, so that the coords are relative to the clipped region.
// constexpr uint32_t FLAG_OVERLAP_ORIENT_IN_STRAND        = ((static_cast<uint32_t>(1)) << 10);   // By default, overlap coordinates are in the fwd strand always. If this flag is true, then the overlap coords are in the strand of alignment.

using OverlapCompactPtr = std::unique_ptr<raptor::OverlapCompact>;

std::unique_ptr<raptor::OverlapCompact> createOverlapCompact();
std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(
            seqid_t Aid, seqid_t Bid,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end,
            bool b_rev, coord_t b_start, coord_t b_end);
std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(const std::unique_ptr<raptor::OverlapCompact>& ovl);

class OverlapCompact {
  public:
    friend std::unique_ptr<raptor::OverlapCompact> createOverlapCompact();
    friend std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(
                seqid_t Aid, seqid_t Bid,
                float score, float identity,
                bool a_rev, coord_t a_start, coord_t a_end,
                bool b_rev, coord_t b_start, coord_t b_end);
    friend std::unique_ptr<raptor::OverlapCompact> createOverlapCompact(const std::unique_ptr<raptor::OverlapCompact>& ovl);
    ~OverlapCompact() = default;

    bool IsFiltered() const {
        // This can evolve in time.
        return (flag_ != 0);
    }

    // Setters.
    void a_id(seqid_t val) { a_id_ = val; }
    void b_id(seqid_t val) { b_id_ = val; }
    void score(float score) { score_ = score; }
    void identity(float val) { identity_ = val; }
    void a_rev(bool val) { a_rev_ = val; }
    void a_start(seqid_t val) { a_start_ = val; }
    void a_end(seqid_t val) { a_end_ = val; }
    void b_rev(bool val) { b_rev_ = val; }
    void b_start(seqid_t val) { b_start_ = val; }
    void b_end(seqid_t val) { b_end_ = val; }
    void flag(uint32_t val) { flag_ = val; }

    // Getters.
    seqid_t a_id() const { return a_id_; }
    seqid_t b_id() const { return b_id_; }
    float score() const { return score_; }
    float identity() const { return identity_; }
    bool a_rev() const { return a_rev_; }
    coord_t a_start() const { return a_start_; }
    coord_t a_end() const { return a_end_; }
    bool b_rev() const { return b_rev_; }
    coord_t b_start() const { return b_start_; }
    coord_t b_end() const { return b_end_; }
    uint32_t flag() const { return flag_; }

    seqid_t a_span() const { return (a_end_ - a_start_); }
    seqid_t b_span() const { return (b_end_ - b_start_); }

    // Flag getters.
    bool GetType5Prime() const {
        return type_ & FLAG_OVERLAP_TYPE_5_PRIME;
    }
    bool GetType3Prime() const {
        return type_ & FLAG_OVERLAP_TYPE_3_PRIME;
    }

    bool GetFlag(uint32_t flag) const {
        return flag_ & flag;
    }

    // Flag setters.
    void SetType5Prime(bool val) {
        SetType_(FLAG_OVERLAP_TYPE_5_PRIME, val);
    }
    void SetType3Prime(bool val) {
        SetType_(FLAG_OVERLAP_TYPE_3_PRIME, val);
    }

    void SetFlagLowLen(bool val) {
        SetFlag_(FLAG_OVERLAP_LOW_LEN, val);
    }
    void SetFlagLowSpan(bool val) {
        SetFlag_(FLAG_OVERLAP_LOW_SPAN, val);
    }
    void SetFlagLowScore(bool val) {
        SetFlag_(FLAG_OVERLAP_LOW_SCORE, val);
    }
    void SetFlagLowIdentity(bool val) {
        SetFlag_(FLAG_OVERLAP_LOW_IDT, val);
    }
    void SetFlagContainedA(bool val) {
        SetFlag_(FLAG_OVERLAP_CONTAINED_A, val);
    }
    void SetFlagContainedB(bool val) {
        SetFlag_(FLAG_OVERLAP_CONTAINED_B, val);
    }
    void SetFlagFilteredA(bool val) {
        SetFlag_(FLAG_OVERLAP_FILTERED_A, val);
    }
    void SetFlagFilteredB(bool val) {
        SetFlag_(FLAG_OVERLAP_FILTERED_B, val);
    }
    void SetFlagInternal(bool val) {
        SetFlag_(FLAG_OVERLAP_INTERNAL, val);
    }
    void SetFlagNotDovetail(bool val) {
        SetFlag_(FLAG_OVERLAP_NOT_DOVETAIL, val);
    }

    // void SetFlagClippedToRegion(bool val) {
    //     SetFlag_(FLAG_OVERLAP_IS_CLIPPED_TO_REGION, val);
    // }
    // void SetFlagOverlapOrientedInStrand(bool val) {
    //     SetFlag_(FLAG_OVERLAP_ORIENT_IN_STRAND, val);
    // }

  private:
    OverlapCompact()
        :
            a_id_(0), b_id_(0), score_(0.0f), identity_(0.0f),
            a_rev_(false), a_start_(0), a_end_(0),
            b_rev_(false), b_start_(0), b_end_(0),
            flag_(0), type_(0)
    { }


    OverlapCompact(seqid_t a_id, seqid_t b_id,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end,
            bool b_rev, coord_t b_start, coord_t b_end)
        :
            a_id_(a_id), b_id_(b_id), score_(score), identity_(identity),
            a_rev_(a_rev), a_start_(a_start), a_end_(a_end),
            b_rev_(b_rev), b_start_(b_start), b_end_(b_end),
            flag_(0),
            type_(0)
    { }

    OverlapCompact(const std::unique_ptr<raptor::OverlapCompact>& ovl)
        :
            a_id_(ovl->a_id_), b_id_(ovl->b_id_), score_(ovl->score_), identity_(ovl->identity_),
            a_rev_(ovl->a_rev_), a_start_(ovl->a_start_), a_end_(ovl->a_end_),
            b_rev_(ovl->b_rev_), b_start_(ovl->b_start_), b_end_(ovl->b_end_),
            flag_(ovl->flag_),
            type_(ovl->type_)
    { }

    void SetFlag_(uint32_t flag, bool val) {
        if (val) {  flag_ |= flag; }
        else {      flag_ &= ~flag; }
    }

    void SetType_(uint32_t flag, bool val) {
        if (val) {  type_ |= flag; }
        else {      type_ &= ~flag; }
    }

    seqid_t a_id_;
    seqid_t b_id_;

    float score_;
    float identity_;

    bool a_rev_;
    coord_t a_start_;
    coord_t a_end_;

    bool b_rev_;
    coord_t b_start_;
    coord_t b_end_;

    uint32_t flag_;

    uint32_t type_;
};

}

#endif
