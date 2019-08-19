/*
 * cigar.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_CIGAR_H_
#define SRC_ALIGNER_CIGAR_H_

#include <algorithm>
#include <cstdint>
#include <string>
#include <vector>
// #include <aligner/aligner_util.hpp>

namespace raptor {

#define IsCigarOp(x) \
    ((x) == 'M' || (x) == '=' || (x) == 'X' || (x) == 'I' || (x) == 'D' || (x) == 'S' || (x) == 'H')
#define IsCigarMatch(x) ((x) == 'M' || (x) == '=' || (x) == 'X')
#define IsCigarIns(x) ((x) == 'I')
#define IsCigarDel(x) ((x) == 'D')
#define IsCigarSoft(x) ((x) == 'S')
#define IsCigarHard(x) ((x) == 'H')
#define IsCigarRef(x) ((x) == 'M' || (x) == '=' || (x) == 'X' || (x) == 'D' || (x) == 'N')
#define IsCigarQuery(x) ((x) == 'M' || (x) == '=' || (x) == 'X' || (x) == 'I' || (x) == 'S')

class CigarOp;
class CigarVector;

std::vector<raptor::CigarOp> CigarStringToVector(const std::string& cigar_str);

int64_t QueryLengthFromCigar(const std::vector<raptor::CigarOp>& split_cigar,
                             bool count_clipping_ops);

int64_t ReferenceLengthFromCigar(const std::vector<raptor::CigarOp>& split_cigar);

std::vector<raptor::CigarOp> ConvertBasicToExtCIGAR(
    const char* qseq, int64_t qlen, const char* tseq, int64_t tlen,
    const std::vector<raptor::CigarOp>& basic_cigar);

int64_t EditDistFromExtCIGAR(const std::vector<raptor::CigarOp>& extended_cigar);

int64_t MatchesFromExtCIGAR(const std::vector<raptor::CigarOp>& extended_cigar);

std::vector<raptor::CigarOp> ExtractCigarBetweenQueryCoords(
    const std::vector<raptor::CigarOp>& cigar, int64_t qstart, int64_t qend);

std::string CigarToString(const std::vector<raptor::CigarOp>& cigar, bool skip_clipping_ops);

std::vector<raptor::CigarOp> AlignmentArrayToCigar(const unsigned char* aln, int aln_len);

std::vector<int8_t> CigarToAlignmentArray(const std::vector<raptor::CigarOp>& cigar);

int64_t ScoreCigarAlignment(const std::vector<raptor::CigarOp>& cigar,
                                        int32_t match, int32_t mismatch,
                                        int32_t gap_open, int32_t gap_ext);

/** @brief A container for a single CIGAR operation.
 *
 */
class CigarOp {
   public:
    CigarOp() : op(0), count(0) {}

    CigarOp(char _op, int32_t _count) : op(_op), count(_count) {}

    ~CigarOp() {}

    char op;
    int32_t count;
};

class CigarVector {
   public:
    CigarVector() = default;
    ~CigarVector() = default;

    const std::vector<raptor::CigarOp>& vec() const { return vec_; }

    void Add(char op, int32_t count) { Add(raptor::CigarOp(op, count)); }

    void Add(const CigarOp& new_op) {
        if (vec_.size() > 0 && vec_.back().op == new_op.op) {
            vec_.back().count += new_op.count;
        } else {
            vec_.emplace_back(new_op);
        }
    }

    void Add(const CigarVector& source) {
        if (source.vec().size() == 0) {
        } else if (vec_.size() == 0) {
            vec_.insert(vec_.end(), source.vec().begin(), source.vec().end());

        } else if (source.vec().front().op == vec_.back().op) {
            vec_.back().count += source.vec().front().count;
            vec_.insert(vec_.end(), source.vec().begin() + 1, source.vec().end());

        } else {
            vec_.insert(vec_.end(), source.vec().begin(), source.vec().end());
        }
    }

    std::string ToString(bool skip_clipping_ops = false) const {
        return raptor::CigarToString(vec_, skip_clipping_ops);
    }

    void FromString(const std::string& cigar_str) { vec_ = raptor::CigarStringToVector(cigar_str); }

    void Reverse() { std::reverse(vec_.begin(), vec_.end()); }

   private:
    std::vector<raptor::CigarOp> vec_;
};

}  // namespace raptor

#endif
