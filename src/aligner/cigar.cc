/*
 * cigar.cc
 *
 *  Created on: Jan 10, 2018
 *      Author: isovic
 */

#include <aligner/cigar.h>
#include <sstream>

namespace raptor {

std::vector<raptor::CigarOp> CigarStringToVector(const std::string& cigar_str) {
    std::vector<CigarOp> ret;

    CigarOp op;
    int64_t pos_ref = 0, pos_query = 0;

    const char* first_digit = NULL;
    for (size_t i = 0; i < cigar_str.size(); i++) {
        if (isalpha(cigar_str[i]) || cigar_str[i] == '=') {
            op.op = cigar_str[i];
            sscanf(first_digit, "%d", &op.count);
            ret.emplace_back(op);
            first_digit = NULL;
            if (IsCigarRef(op.op)) pos_ref += op.count;
            if (IsCigarQuery(op.op)) pos_query += op.count;
        } else if (first_digit == NULL) {
            first_digit = &(cigar_str[i]);
        }
    }

    return ret;
}

int64_t QueryLengthFromCigar(const std::vector<raptor::CigarOp>& split_cigar,
                             bool count_clipping_ops) {
    int64_t len = 0;
    if (count_clipping_ops) {
        for (size_t i = 0; i < split_cigar.size(); i++) {
            if (IsCigarQuery(split_cigar[i].op)) {
                len += split_cigar[i].count;
            }
        }
    } else {
        for (size_t i = 0; i < split_cigar.size(); i++) {
            if (IsCigarQuery(split_cigar[i].op) && split_cigar[i].op != 'S') {
                len += split_cigar[i].count;
            }
        }
    }
    return len;
}

int64_t ReferenceLengthFromCigar(const std::vector<raptor::CigarOp>& split_cigar) {
    int64_t len = 0;
    for (size_t i = 0; i < split_cigar.size(); i++) {
        if (IsCigarRef(split_cigar[i].op)) {
            len += split_cigar[i].count;
        }
    }
    return len;
}

std::vector<raptor::CigarOp> ConvertBasicToExtCIGAR(
    const char* qseq, int64_t qlen, const char* tseq, int64_t tlen,
    const std::vector<raptor::CigarOp>& basic_cigar) {
    std::vector<raptor::CigarOp> ret;

    int64_t qpos = 0, tpos = 0;
    for (size_t i = 0; i < basic_cigar.size(); i++) {
        char op = basic_cigar[i].op;
        int64_t count = basic_cigar[i].count;

        if (op != 'M') {
            ret.push_back(basic_cigar[i]);

            if (op == 'I' || op == 'S') {
                qpos += count;
            }
            if (op == 'D' || op == 'N') {
                tpos += count;
            }
        } else {
            char prev_m = 0;
            int64_t curr_count = 0;
            for (int64_t j = 0; j < count; j++) {
                char curr_m = (qseq[qpos] == tseq[tpos]) ? '=' : 'X';
                if (j == 0) {
                    prev_m = curr_m;
                }
                if (curr_m == prev_m) {
                    curr_count += 1;
                } else {
                    ret.push_back(raptor::CigarOp(prev_m, curr_count));
                    prev_m = curr_m;
                    curr_count = 1;
                }
                qpos += 1;
                tpos += 1;
            }
            if (curr_count > 0) {
                ret.push_back(raptor::CigarOp(prev_m, curr_count));
            }
        }
    }

    return ret;
}

int64_t EditDistFromExtCIGAR(const std::vector<raptor::CigarOp>& extended_cigar) {
    int64_t edit_dist = 0;
    for (size_t i = 0; i < extended_cigar.size(); i++) {
        char op = extended_cigar[i].op;
        if (op == 'M') {
            fprintf(stderr, "CIGAR is not in the extended format.");
            return -1;
        }
        if (op == 'X' || op == 'I' || op == 'D') {
            edit_dist += extended_cigar[i].count;
        }
    }
    return edit_dist;
}

int64_t MatchesFromExtCIGAR(const std::vector<raptor::CigarOp>& extended_cigar) {
    int64_t ret = 0;
    for (size_t i = 0; i < extended_cigar.size(); i++) {
        char op = extended_cigar[i].op;
        if (op == 'M') {
            fprintf(stderr, "CIGAR is not in the extended format.");
            return -1;
        }
        if (op == '=') {
            ret += extended_cigar[i].count;
        }
    }
    return ret;
}

std::vector<raptor::CigarOp> ExtractCigarBetweenQueryCoords(
    const std::vector<raptor::CigarOp>& cigar, int64_t qstart, int64_t qend) {
    std::vector<raptor::CigarOp> ret;

    // printf ("qstart = %ld, qend = %ld\n", qstart, qend);

    int64_t qpos = 0;

    for (auto& c : cigar) {
        int64_t qpos_next =
            (c.op == 'M' || c.op == '=' || c.op == 'X' || c.op == 'I' || c.op == 'S')
                ? (qpos + c.count)
                : qpos;

        // printf ("\n");
        // printf ("(1) Entered: c = %ld%c\n", c.count, c.op);
        // printf ("(2) qpos = %ld, qpos_next = %ld\n", qpos, qpos_next);

        if (qpos > qend) {
            break;
        }

        if (qpos_next < qstart) {
            qpos = qpos_next;
            continue;
        }

        int64_t b = 0, e = c.count;

        if (qstart >= qpos && qstart < qpos_next) {
            b = qstart - qpos;
        }
        if (qend >= qpos && qend < qpos_next) {
            e = qend - qpos;
        }

        if ((e - b) > 0) {
            ret.emplace_back(raptor::CigarOp(c.op, (e - b)));
        }

        qpos = qpos_next;
    }

    return ret;
}

std::string CigarToString(const std::vector<raptor::CigarOp>& cigar, bool skip_clipping_ops) {
    std::ostringstream ss;
    for (size_t i = 0; i < cigar.size(); i++) {
        if (skip_clipping_ops && (cigar[i].op == 'S' || cigar[i].op == 'H')) {
            continue;
        }

        ss << cigar[i].count << cigar[i].op;
    }
    return ss.str();
}

std::vector<raptor::CigarOp> AlignmentArrayToCigar(const unsigned char* aln, int aln_len) {
    std::vector<raptor::CigarOp> ret;

    if (aln_len <= 0) {
        return ret;
    }

    char prev_op = 0;  // Char of last move. 0 if there was no previous move.
    int count = 0;
    for (int i = 0; i <= aln_len; i++) {
        if (i == aln_len || (aln[i] != prev_op && prev_op != 0)) {
            ret.emplace_back(raptor::CigarOp(prev_op, count));
            count = 0;
        }
        if (i < aln_len) {
            prev_op = aln[i];
            count += 1;
        }
    }
    return ret;
}

std::vector<int8_t> CigarToAlignmentArray(const std::vector<raptor::CigarOp>& cigar) {
    int64_t len = 0;
    for (size_t i = 0; i < cigar.size(); i++) {
        len += cigar[i].count;
    }

    std::vector<int8_t> ret(len, 0);
    int64_t pos = 0;
    for (size_t i = 0; i < cigar.size(); i++) {
        for (size_t j = 0; j < cigar[i].count; j++) {
            ret[pos++] = cigar[i].op;
        }
    }
    return ret;
}

}  // namespace raptor
