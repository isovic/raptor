#include <aligner/sam_parser.h>

#include <sstream>

namespace raptor {

SamLine::SamLine() {}

SamLine::SamLine(const std::string& line) { ParseLine(line); }

SamLine::~SamLine() {}

int SamLine::ParseLine(const std::string& line) {
    if (line.size() == 0) {
        return 1;
    }

    std::stringstream ss(line);
    std::string cigar_string;
    ss >> qname >> flag >> rname >> pos >> mapq >> cigar_string >> rnext >> pnext >> tlen >> seq >>
        qual;

    cigar = CigarStringToVector(cigar_string);

    std::string all_optional;
    std::getline(ss, all_optional);
    Tokenize_(all_optional, '\t', optional);
    return 0;
}

bool SamLine::IsMapped() { return (!(flag & 4)); }

bool SamLine::IsReverse() { return ((flag & 16)); }

int SamLine::FindAlignmentPosition(int64_t& q_start, int64_t& q_end, int64_t& r_start,
                                   int64_t& r_end) {
    q_start = 0;
    q_end = seq.size();

    // Find query alignment start (skip the soft clipped bases).
    for (auto& c : cigar) {
        if (c.op == 'H') {
            continue;
        } else if (c.op == 'S') {
            q_start += c.count;
        } else {
            break;
        }
    }

    // Find query alignment end (skip the soft clipped bases).
    for (int64_t i = (cigar.size() - 1); i >= 0; i--) {
        auto& c = cigar[i];
        if (c.op == 'H') {
            continue;
        } else if (c.op == 'S') {
            q_end -= c.count;
        } else {
            break;
        }
    }

    // Find reference alignment start. (Convert from 1-based to 0-based).
    r_start = pos - 1;

    // Find reference alignment end.
    r_end = r_start + ReferenceLengthFromCigar(cigar);

    // Do not performe reverse complementing here, we do not know
    // the length of the reference.

    return 0;
}

void SamLine::Tokenize_(const std::string& str, const char delimiter,
                        std::vector<std::string>& words) {
    words.clear();
    std::stringstream ss(str);
    std::string line;
    while (std::getline(ss, line, delimiter)) {
        if (line.size() == 0) {
            continue;
        }
        words.push_back(line);
    }
}

}  // namespace raptor
