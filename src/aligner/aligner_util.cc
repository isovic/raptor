#include <aligner/aligner_util.hpp>

#include <sstream>
#include <string>
#include <log/log_tools.h>
#include <lib/edlib.h>
#include <aligner/cigar.h>

namespace raptor {

std::vector<int8_t> ConvertSeqAlphabet(const int8_t* seq, size_t seqlen,
                                       const uint8_t* conv_table) {
    std::vector<int8_t> ret(seqlen + 33);  // 32 for gaba
    for (size_t i = 0; i < seqlen; i++) {
        ret[i] = (int8_t)conv_table[(uint8_t)seq[i]];
    }
    return ret;
}

std::string EdlibAlignmentToPrettyString(const unsigned char* query, const int queryLength,
                                         const unsigned char* target, const int targetLength,
                                         const std::shared_ptr<raptor::AlignmentResult>& aln,
                                         int32_t row_width) {
    std::ostringstream ss;

    std::vector<int8_t> aln_array = CigarToAlignmentArray(aln->cigar());

    return ss.str();
}

// std::string EdlibAlignmentToPrettyString(const unsigned char* query, const int queryLength,
//                                          const unsigned char* target, const int targetLength,
//                                          const unsigned char* alignment, const int alignmentLength,
//                                          const int position, const int modeCode, int row_width) {
//     std::stringstream ss;

//     int tIdx = -1;
//     int qIdx = -1;
//     if (modeCode == EDLIB_MODE_HW) {
//         tIdx = position;
//         for (int i = 0; i < alignmentLength; i++) {
//             if (alignment[i] != EDLIB_I && alignment[i] != EDLIB_S) tIdx--;
//         }
//     }

//     for (int start = 0; start < alignmentLength; start += row_width) {
//         // target
//         ss << "T: ";
//         int startTIdx = 0;
//         for (int j = start; j < start + row_width && j < alignmentLength; j++) {
//             if (alignment[j] == EDLIB_I || alignment[j] == EDLIB_S) {
//                 ss << "_";

//             } else {
//                 ss << (char)target[++tIdx];
//             }
//             if (j == start) startTIdx = tIdx;
//         }
//         ss << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")\n";

//         // Mismatch
//         ss << "   ";
//         for (int j = start; j < start + row_width && j < alignmentLength; j++) {
//             if (alignment[j] == 0)
//                 ss << "|";
//             else if (alignment[j] == EDLIB_X)
//                 ss << "X";
//             else if (alignment[j] == EDLIB_S)
//                 ss << "-";
//             else
//                 ss << " ";
//         }
//         ss << "\n";

//         // query
//         ss << "Q: ";
//         int startQIdx = qIdx;
//         for (int j = start; j < start + row_width && j < alignmentLength; j++) {
//             if (alignment[j] == EDLIB_D) {
//                 ss << "_";
//             } else {
//                 ss << (char)query[++qIdx];
//             }
//             if (j == start) startQIdx = qIdx;
//         }
//         ss << " (" << std::max(startQIdx, 0) << " - " << qIdx << ")\n\n";
//     }

//     return ss.str();
// }

}  // namespace raptor
