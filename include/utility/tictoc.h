/*
 * TicToc.h
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_TicToc2_H_
#define SRC_CONSENSUS_TicToc2_H_

#include <time.h>

class TicToc {
 public:
  TicToc();
  ~TicToc();

  void start();
  void stop();
  double get_secs();
  double get_msecs();
  double get_secs_current();

 private:
  clock_t start_;
  clock_t end_;
};

#endif /* SRC_CONSENSUS_TicToc_H_ */
