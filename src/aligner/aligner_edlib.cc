#include <aligner/aligner_edlib.h>
// #include "libs/edlibcigar.h"

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <aligner/aligner_util.hpp>

namespace raptor {

std::shared_ptr<AlignerBase> createAlignerEdlib(const raptor::AlignmentOptions& opt) {
    return std::shared_ptr<AlignerBase>(new AlignerEdlib(opt));
}

AlignerEdlib::AlignerEdlib(const raptor::AlignmentOptions& opt) : opt_(opt) {}

AlignerEdlib::~AlignerEdlib() {}

std::vector<raptor::CigarOp> EdlibAlignmentToCigar(unsigned char* aln, int aln_len) {
    std::vector<raptor::CigarOp> ret;

    if (aln_len <= 0) {
        return ret;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char op_to_cigar[] = {'=', 'I', 'D', 'X'};

    char prev_op = 0;  // Char of last move. 0 if there was no previous move.
    int count = 0;
    for (int i = 0; i <= aln_len; i++) {
        if (i == aln_len || (op_to_cigar[aln[i]] != prev_op && prev_op != 0)) {
            ret.emplace_back(raptor::CigarOp(prev_op, count));
            count = 0;
        }
        if (i < aln_len) {
            prev_op = op_to_cigar[aln[i]];
            count += 1;
        }
    }
    return ret;
}

std::shared_ptr<raptor::AlignmentResult> AlignerEdlib::Global(const char* qseq, int64_t qlen,
                                                              const char* tseq, int64_t tlen) {
    return RunEdlib(qseq, qlen, tseq, tlen, opt_, EDLIB_MODE_NW, EDLIB_TASK_PATH);
}

std::shared_ptr<raptor::AlignmentResult> AlignerEdlib::Extend(const char* qseq, int64_t qlen,
                                                              const char* tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result = RunEdlib(qseq, qlen, tseq, tlen, opt_, EDLIB_MODE_SHW, EDLIB_TASK_PATH);

    auto aln_vec = CigarToAlignmentArray(result->cigar());

    int64_t score = (opt_.end_bonus >= 0) ? opt_.end_bonus : 0;
    // bool stop_ext_on_zero_score = (opt_.zdrop == -99912345) ? false : true;         // Emulates local alignment if true. This might prevent end-to-end extensions, so allow a hidden feature to enable end-to-end extensions. These can severely reduce alignment score / identity.

    int64_t max_score = -1000000000;  // Large negative value;
    int64_t max_q = 0;
    int64_t max_t = 0;
    int64_t drop_q = 0, drop_t = 0, drop_score = 0;
    int64_t qpos = 0, tpos = 0;
    size_t final_i = 0;
    int64_t edit_dist = 0;
    bool stop_loop = false;
    char prev_op = 'Z';
    int32_t streak = 1;
    for (size_t i = 0; i < aln_vec.size(); ++i) {
        char op = aln_vec[i];
        streak = (op == prev_op) ? (streak + 1) : 1;
        // int32_t count = 1; // streak;
        final_i = i;
        // std::cerr << "op = " << op << ", prev_op = " << prev_op << ", count = " << 1 << ", i = " << i << ", score = " << score << "\n";
        prev_op = op;

        switch (op) {
            case '=':
                score += opt_.p.match;
                qpos += 1;
                tpos += 1;
                break;
            case 'X':
                score += opt_.p.mismatch;
                qpos += 1;
                tpos += 1;
                edit_dist += 1;
                break;
            case 'I':
                if (streak == 1) {
                    score += opt_.p.w[0].open;         // Gap open only if a new streak is found.
                }
                score += opt_.p.w[0].ext;  // This is the gap extend penalty.
                qpos += 1;
                edit_dist += 1;
                break;
            case 'D':
                if (streak == 1) {
                    score += opt_.p.w[0].open;         // Gap open only if a new streak is found.
                }
                score += opt_.p.w[0].ext;  // This is the gap extend penalty.
                tpos += 1;
                edit_dist += 1;
                break;
            default:
                // If an unknown op is found, don't do zdrop.
                max_score = -1000000000;
                qpos = qlen;
                tpos = tlen;
                score = result->score();
                final_i = result->cigar().size();
                stop_loop = true;
                break;
        }
        if (score > max_score) {
            max_score = score;
            max_q = qpos;
            max_t = tpos;
        }
        if (qpos >= qlen || tpos >= tlen) {
            // std::cerr << "Break 1!\n";
            break;
        }
        if (opt_.stop_ext_on_zero_score && score < 0) {
            // std::cerr << "Break 2!\n";
            break;
        }
        if (stop_loop == true) {
            // std::cerr << "Break 3!\n";
            break;
        }
        if (opt_.zdrop >= 0 && score < (max_score - opt_.zdrop)) {
            // std::cerr << "Break 4!\n";
            break;
        }
    }

    if (final_i == aln_vec.size()) {
        return result;
    }

    // std::cerr << "max_score = " << max_score << ", score = " << score << "\n";
    // std::cerr << "cigar = " << raptor::CigarToString(result->cigar(), false) << "\n";

    std::vector<raptor::CigarOp> new_cigar =
        AlignmentArrayToCigar((const unsigned char*)&aln_vec[0], final_i + 1);
    result->cigar(new_cigar);
    result->score(score);
    result->position(raptor::AlignmentPosition(0, qpos, 0, tpos));
    result->max_score(max_score);
    result->max_q_pos(max_q);
    result->max_t_pos(max_t);
    result->edit_dist(edit_dist);
    result->final_band(-1);
    result->status(raptor::AlignmentReturnValue::OK);

    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerEdlib::Local(const char* qseq, int64_t qlen,
                                                             const char* tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::NotImplementedYet);
    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerEdlib::Semiglobal(const char* qseq, int64_t qlen,
                                                                  const char* tseq, int64_t tlen) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::NotImplementedYet);
    return result;
}

std::shared_ptr<raptor::AlignmentResult> AlignerEdlib::RunEdlib(const char* qseq, int64_t qlen,
                                                                const char* tseq, int64_t tlen,
                                                                const raptor::AlignmentOptions& opt,
                                                                EdlibAlignMode edlib_mode,
                                                                EdlibAlignTask edlib_task) {
    auto result = raptor::createAlignmentResult();

    if (qseq == NULL || tseq == NULL || qlen <= 0 || tlen <= 0) {
        result->status(raptor::AlignmentReturnValue::InvalidOptions);
        return result;
    }

    result->status(raptor::AlignmentReturnValue::AlignerFailure);

    EdlibAlignResult edlib_result =
        edlibAlign((const char*)qseq, qlen, (const char*)tseq, tlen,
                   edlibNewAlignConfig(opt.bandwidth, edlib_mode, edlib_task));

    if (edlib_result.numLocations == 0) {
        edlibFreeAlignResult(edlib_result);
        result->status(raptor::AlignmentReturnValue::AlignerFailure);
        return result;
    }

    // result->score(-edlib_result.editDistance);
    result->edit_dist(edlib_result.editDistance);
    result->position(raptor::AlignmentPosition(0, qlen, 0, tlen));
    result->final_band(-1);
    result->status(raptor::AlignmentReturnValue::OK);
    result->max_q_pos(qlen);
    result->max_t_pos(tlen);
    result->cigar(EdlibAlignmentToCigar(edlib_result.alignment, edlib_result.alignmentLength));

    int64_t score = ScoreCigarAlignment(result->cigar(), opt.p.match, opt.p.mismatch, opt.p.w[0].open, opt.p.w[0].ext);

    result->score(score);
    result->max_score(result->score());

    edlibFreeAlignResult(edlib_result);

    return result;
}

}  // namespace raptor
