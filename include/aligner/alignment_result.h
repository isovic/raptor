/*
 * alignment_result.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNMENT_RESULT_H_
#define SRC_ALIGNER_ALIGNMENT_RESULT_H_

#include <stdint.h>
#include <vector>
#include <memory>
#include <aligner/consts.h>
#include <aligner/alignment_position.h>
#include <aligner/cigar.h>

namespace raptor {

class AlignmentResult;

enum class AlignmentReturnValue {  // Strongly typed enum, C++11 feature.
    OK,                            // Everything went ok.
    Suboptimal,                    // Alignment stepped out of defined band. Result is not optimal.
    InvalidOptions,                // In case parameters of values are invalid.
    QlenIsZero,
    TlenIsZero,
    WrongEditDist,
    AlignmentNotPerformed,  // A default value for an alignment which wasn't performed.
    AlignerFailure,
    NotImplementedYet  // For features in development.
};

class CigarOpCounts {
   public:
    CigarOpCounts() = default;
    int64_t eq = 0;  // Num matching bases (equal operations).
    int64_t x = 0;   // Mismatches.
    int64_t i = 0;   // Insertions.
    int64_t d = 0;   // Deletions.
    int64_t iu = 0;  // Insertions, where each streak of CIGAR operations is counted as 1 (unique event).
    int64_t du = 0;  // Deletions, where each streak of CIGAR operations is counted as 1 (unique event).
    int64_t n = 0;   // Num 'N' bases.
    int64_t s = 0;   // Soft clipped bases.
    int64_t h = 0;   // Hard clipped bases.
    int64_t aligned_qlen = 0;  // Length of the alignment on the query sequence.
    int64_t aligned_tlen = 0;  // Length of the alignment on the target sequence.
    int64_t num_ops = 0;
    double identity_q = 0.0;
    double identity_t = 0.0;
    double identity_min = 0.0;
    double error_rate = 0.0;
    double error_rate_u = 0.0;
};

std::shared_ptr<raptor::AlignmentResult> createAlignmentResult();

class AlignmentResult {
   public:
    friend std::shared_ptr<raptor::AlignmentResult> createAlignmentResult();
    ~AlignmentResult();

    std::string Verbose() const;

    /*
     * Getters.
    */
    int64_t score() const { return score_; }
    int64_t edit_dist() const { return edit_dist_; }
    int64_t max_score() const { return max_score_; }
    int64_t max_q_pos() const { return max_q_pos_; }
    int64_t max_t_pos() const { return max_t_pos_; }
    int64_t final_band() const { return final_band_; }
    const raptor::AlignmentPosition& position() const { return position_; }
    AlignmentReturnValue status() const { return status_; }
    const std::vector<raptor::CigarOp>& cigar() const { return cigar_; }
    const CigarOpCounts& op_counts() const { return op_counts_; }

    /*
     * Setters.
    */
    void score(int64_t _score) { score_ = _score; }
    void edit_dist(int64_t _edit_dist) { edit_dist_ = _edit_dist; }
    void max_score(int64_t _max_score) { max_score_ = _max_score; }
    void max_q_pos(int64_t _max_q_pos) { max_q_pos_ = _max_q_pos; }
    void max_t_pos(int64_t _max_t_pos) { max_t_pos_ = _max_t_pos; }
    void final_band(int64_t _final_band) { final_band_ = _final_band; }
    void position(const raptor::AlignmentPosition& _position) { position_ = _position; }
    void status(AlignmentReturnValue _status) { status_ = _status; }
    void cigar(const std::vector<raptor::CigarOp>& _cigar) {
        cigar_ = _cigar;
        op_counts_ = CigarOpCounts();
        for (size_t i = 0; i < _cigar.size(); i++) {
            op_counts_.num_ops += _cigar[i].count;
            if (_cigar[i].op == '=') {
                op_counts_.eq += _cigar[i].count;
                op_counts_.aligned_qlen += _cigar[i].count;
                op_counts_.aligned_tlen += _cigar[i].count;
            } else if (_cigar[i].op == 'X') {
                op_counts_.x += _cigar[i].count;
                op_counts_.aligned_qlen += _cigar[i].count;
                op_counts_.aligned_tlen += _cigar[i].count;
            } else if (_cigar[i].op == 'I') {
                op_counts_.iu += 1;
                op_counts_.i += _cigar[i].count;
                op_counts_.aligned_qlen += _cigar[i].count;
            } else if (_cigar[i].op == 'D') {
                op_counts_.du += 1;
                op_counts_.d += _cigar[i].count;
                op_counts_.aligned_tlen += _cigar[i].count;
            } else if (_cigar[i].op == 'N') {
                op_counts_.n += _cigar[i].count;
            } else if (_cigar[i].op == 'S') {
                op_counts_.s += _cigar[i].count;
            } else if (_cigar[i].op == 'H') {
                op_counts_.h += _cigar[i].count;
            }
        }
        double max_len = (double)std::max(op_counts_.aligned_qlen, op_counts_.aligned_tlen);
        op_counts_.identity_q =
            100.0 * ((op_counts_.aligned_qlen == 0) ? 1.0 : ((double)(op_counts_.eq)) / op_counts_.aligned_qlen);
        op_counts_.identity_t =
            100.0 * ((op_counts_.aligned_qlen == 0) ? 1.0 : ((double)(op_counts_.eq)) / op_counts_.aligned_tlen);
        op_counts_.identity_min = std::min(op_counts_.identity_q, op_counts_.identity_t);
        op_counts_.error_rate = (op_counts_.num_ops == 0)
                                    ? 0.0
                                    : ((double)(op_counts_.x + op_counts_.i + op_counts_.d)) /
                                          ((double)op_counts_.num_ops);
        op_counts_.error_rate_u =
            ((double)(op_counts_.x + op_counts_.iu + op_counts_.du)) /
            ((double)(op_counts_.eq + op_counts_.x + op_counts_.iu + op_counts_.du));
    }

   private:
    AlignmentResult();
    AlignmentResult(const AlignmentResult&) = delete;
    AlignmentResult& operator=(const AlignmentResult&) = delete;

    int64_t score_;
    int64_t edit_dist_;
    raptor::AlignmentPosition position_;  // There can be multiple alignments with the same score.
                                          // Only the first position and the corresponding alignment
    std::vector<raptor::CigarOp> cigar_;  // are reported
    int64_t max_score_;
    int64_t max_q_pos_;
    int64_t max_t_pos_;  // Maximum score in the alignment, and the coordinates on query and target.
    int64_t final_band_;           // Value of band k used in the final alignment.
    AlignmentReturnValue status_;  // Return value of the aligner.
    CigarOpCounts op_counts_;
};

}  // namespace raptor

#endif
