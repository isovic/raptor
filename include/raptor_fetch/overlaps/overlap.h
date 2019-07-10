#ifndef _IPA_INCLUDE_OVERLAP_H_
#define _IPA_INCLUDE_OVERLAP_H_

#include <cstdint>
#include <memory>
#include <string>

namespace raptor {

class Overlap;

using seqid_t = int32_t;
using coord_t = int32_t;

using OverlapPtr = std::unique_ptr<raptor::Overlap>;

std::unique_ptr<raptor::Overlap> createOverlap();
std::unique_ptr<raptor::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len);
std::unique_ptr<raptor::Overlap> createOverlap(
            const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
            coord_t clip_a_start, coord_t clip_a_end,
            coord_t clip_b_start, coord_t clip_b_end,
            const std::string& type);

class Overlap {
  public:
    friend std::unique_ptr<raptor::Overlap> createOverlap();
    friend std::unique_ptr<raptor::Overlap> createOverlap(
                const std::string& a_name, const std::string& b_name,
                float score, float identity,
                bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
                bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len);
    friend std::unique_ptr<raptor::Overlap> createOverlap(
                const std::string& a_name, const std::string& b_name,
                float score, float identity,
                bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
                bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
                coord_t clip_a_start, coord_t clip_a_end,
                coord_t clip_b_start, coord_t clip_b_end,
                const std::string& type);
    ~Overlap() = default;

    void FlipBCoords() {
		std::swap(b_start_, b_end_);
		b_start_ = b_len_ - b_start_;
		b_end_ = b_len_ - b_end_;
    }

    // Getters.
    const std::string& a_name() const { return a_name_; }
    const std::string& b_name() const { return b_name_; }
    float score() const { return score_; }
    float identity() const { return identity_; }
    bool a_rev() const { return a_rev_; }
    coord_t a_start() const { return a_start_; }
    coord_t a_end() const { return a_end_; }
    coord_t a_len() const { return a_len_; }
    bool b_rev() const { return b_rev_; }
    coord_t b_start() const { return b_start_; }
    coord_t b_end() const { return b_end_; }
    coord_t b_len() const { return b_len_; }
    coord_t clip_a_start() const { return clip_a_start_; }
    coord_t clip_a_end() const { return clip_a_end_; }
    coord_t clip_b_start() const { return clip_b_start_; }
    coord_t clip_b_end() const { return clip_b_end_; }
    seqid_t a_id() const { return a_id_; }
    seqid_t b_id() const { return b_id_; }

    void a_id(seqid_t _a_id) { a_id_ = _a_id; }
    void b_id(seqid_t _b_id) { b_id_ = _b_id; }

    seqid_t a_span() const { return (a_end_ - a_start_); }
    seqid_t b_span() const { return (b_end_ - b_start_); }

  private:
    Overlap()
        :
            a_name_(0), b_name_(0), score_(0), identity_(0),
            a_rev_(0), a_start_(0), a_end_(0), a_len_(0),
            b_rev_(0), b_start_(0), b_end_(0), b_len_(0),
            clip_a_start_(0), clip_a_end_(0),
            clip_b_start_(0), clip_b_end_(0),
            type_(),
            a_id_(0), b_id_(0)
    { }

    Overlap(const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len)
        :
            a_name_(a_name), b_name_(b_name), score_(score), identity_(identity),
            a_rev_(a_rev), a_start_(a_start), a_end_(a_end), a_len_(a_len),
            b_rev_(b_rev), b_start_(b_start), b_end_(b_end), b_len_(b_len),
            clip_a_start_(0), clip_a_end_(a_len),
            clip_b_start_(0), clip_b_end_(b_len),
            type_(),
            a_id_(0), b_id_(0)
    { }

    Overlap(const std::string& a_name, const std::string& b_name,
            float score, float identity,
            bool a_rev, coord_t a_start, coord_t a_end, coord_t a_len,
            bool b_rev, coord_t b_start, coord_t b_end, coord_t b_len,
            coord_t clip_a_start, coord_t clip_a_end,
            coord_t clip_b_start, coord_t clip_b_end,
            const std::string& type)
        :
            a_name_(a_name), b_name_(b_name), score_(score), identity_(identity),
            a_rev_(a_rev), a_start_(a_start), a_end_(a_end), a_len_(a_len),
            b_rev_(b_rev), b_start_(b_start), b_end_(b_end), b_len_(b_len),
            clip_a_start_(clip_a_start), clip_a_end_(clip_a_end),
            clip_b_start_(clip_b_start), clip_b_end_(clip_b_end),
            type_(),
            a_id_(0), b_id_(0)
    { }

    std::string a_name_;
    std::string b_name_;

    float score_;
    float identity_;

    bool a_rev_;
    coord_t a_start_;
    coord_t a_end_;
    coord_t a_len_;

    bool b_rev_;
    coord_t b_start_;
    coord_t b_end_;
    coord_t b_len_;

    coord_t clip_a_start_;
    coord_t clip_a_end_;
    coord_t clip_b_start_;
    coord_t clip_b_end_;

    std::string type_;

    seqid_t a_id_;
    seqid_t b_id_;
};

}

#endif
