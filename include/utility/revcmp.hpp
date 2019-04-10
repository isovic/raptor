/*
 * revcmp.hpp
 *
 *  Created on: Dec 9, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTILITY_REVCMP_HPP_
#define SRC_UTILITY_REVCMP_HPP_

#include <stdint.h>
#include <string>

#include <index/lookups.h>

namespace raptor {

inline std::string ReverseComplement(const std::string& in_seq);
inline std::string ReverseComplement(const char* in_seq, size_t seq_len);

inline std::string ReverseComplement(const std::string& in_seq) {
  return ReverseComplement(in_seq.c_str(), in_seq.size());
}

inline std::string ReverseComplement(const char* in_seq, size_t seq_len) {
  std::string ret(seq_len, ' ');
  for (int64_t i=0; i<(int64_t) seq_len; i++) {
    ret[i] = (char) mindex::nuc_to_complement[static_cast<int32_t>(in_seq[seq_len-i-1])];
  }
  return ret;
}

}

#endif
