/*
 * aligner_util.hpp
 *
 *  Created on: Jul 1, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_ALIGNER_ALIGNER_UTIL_H_
#define SRC_ALIGNER_ALIGNER_UTIL_H_

#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <aligner/sam_parser.h>

namespace raptor {

template <typename T>
std::vector<T> GenerateSimpleMatchMatrix(T match, T mismatch, size_t alphabet_size) {
    std::vector<T> matrix(alphabet_size * alphabet_size, mismatch);  // Set the mismatch score.
    // Goes to "-1" to allow for 'N' bases which should not match to themselves.
    for (size_t i = 0; i < (alphabet_size - 1); i++) {
        matrix[i * alphabet_size + i] = match;                // Set the match score.
        matrix[i * alphabet_size + alphabet_size - 1] = 0;    // Reset the last column to 0.
        matrix[(alphabet_size - 1) * alphabet_size + i] = 0;  // Reset the last row to 0.
    }
    return matrix;
}

std::vector<int8_t> ConvertSeqAlphabet(const int8_t* seq, size_t seqlen, const uint8_t* conv_table);

// std::string EdlibAlignmentToPrettyString(const unsigned char* query, const int queryLength,
//                                          const unsigned char* target, const int targetLength,
//                                          const unsigned char* alignment, const int alignmentLength,
//                                          const int position, const int modeCode, int row_width);

}  // namespace raptor

#endif
