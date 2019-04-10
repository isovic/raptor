/*
 * lis.hpp
 *
 *      Author: Ivan Sovic
 *      GitHub: @isovic
 *      Copyright: Ivan Sovic, 2017
 *      Licence: MIT
 *
 * A generic implementation of the Longest Increasing Subsequence algorithm
 * which allows for custom data types, provided that a suitable comparison
 * function is given.
 */

#ifndef SRC_LIS_H_
#define SRC_LIS_H_

#include <cstdint>
#include <algorithm>
#include <vector>
#include <functional>

namespace raptor {

template<class T>
std::vector<T> LIS(const std::vector<T> &v, size_t begin, size_t end,
                   std::function<bool(const T& a,
                                      const T& b)> comp =
                                      [](const T& a, const T& b)
                                      { return a < b; } ) {
  std::vector<T> lis;

  if (v.size() == 0) {
    return lis;
  }

  if (end < begin) {
    return lis;
  }

  size_t n = end - begin;

  std::vector<int64_t> tail(n + 1, 0);
  std::vector<int64_t> pred(n + 1, 0);
  int64_t len = 0;

  for (size_t i = 1; i < (n + 1); i++) {
    int64_t pred_id = 0;
    if (len == 0 || !comp(v[tail[1] + begin - 1], v[i + begin - 1])) {
        tail[0] = i;
        pred_id = 0;
    }
    else {
        // Find the longest LIS which ends in a value
        // smaller than the current one.
        int32_t left = 1;
        int32_t right = len + 1;
        while ((right - left) > 1) {
          int32_t mid = (right + left) / 2;
          if (comp(v[tail[mid] + begin - 1], v[i + begin - 1])) {
            left = mid;
          } else {
            right = mid;
          }
        }
        pred_id = left;
    }
    // Update the DP tables.
    pred[i] = tail[pred_id];
    if (pred_id == len || comp(v[i + begin - 1], v[tail[pred_id + 1] + begin - 1])) {
      tail[pred_id + 1] = i;
      len = std::max(len, pred_id + 1);
    }
  }

  // Backtrack.
  lis.reserve(len);
  int32_t pos = tail[len];
  int32_t curr_len = len;
  while (curr_len > 0) {
    lis.push_back(v[pos + begin - 1]);
    pos = pred[pos];
    curr_len -= 1;
  }
  std::reverse(lis.begin(), lis.end());

  return lis;
}

}

#endif
