/*
 * range.h
 *
 *  Created on: Jul 16, 2014
 *      Author: Ivan Sovic
 */

#ifndef SRC_RANGE_H_
#define SRC_RANGE_H_

#include <cstdint>

namespace raptor {

class Range {
 public:
  Range() : start(0), end(0) { }
  Range(int64_t _start, int64_t _end) : start(_start), end(_end) { }

  int64_t dist() const {
    return (end - start);
  }

  int64_t start = 0;
  int64_t end = 0;
};

}

#endif /* RANGE_H_ */
