#ifndef SRC_SAM_PARSER_H_
#define SRC_SAM_PARSER_H_

#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <aligner/aligner_containers.h>
#include <aligner/cigar.h>

namespace raptor {

class SamLine {
   public:
    SamLine();
    SamLine(const std::string& line);
    ~SamLine();

    int ParseLine(const std::string& line);
    std::string YieldString();
    bool IsMapped();
    bool IsReverse();
    int FindAlignmentPosition(int64_t& q_start, int64_t& q_end, int64_t& r_start, int64_t& r_end);

    std::string qname;  // Field #1.
    uint32_t flag;      // Field #2.
    std::string rname;  // Field #3.
    int64_t pos;        // Field #4.
    int32_t mapq;       // Field #5.
    //  std::string cigar;  // Field #6.
    std::vector<raptor::CigarOp> cigar;
    std::string rnext;  // Field #7.
    int64_t pnext;      // Field #8.
    int64_t tlen;       // Field #9.
    std::string seq;    // Field #10.
    std::string qual;   // Field #11.

    // Optional fields in the SAM format:
    int64_t as;     // Alignment score.
    double evalue;  // E-value. There is no dedicated field in the SAM format, but GraphMap uses ZE
                    // to specify the E-value.
    std::vector<std::string> optional;  // Raw values (strings) of optional fields, not explicitly
                                        // converted to expected values;

   private:
    void Tokenize_(const std::string& str, const char delimiter, std::vector<std::string>& words);
};

}  // namespace raptor

#endif
