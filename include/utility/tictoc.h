/*
 * TicToc.h
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_TicToc2_H_
#define SRC_CONSENSUS_TicToc2_H_

#include <time.h>
#include <chrono>

class TicToc {
  public:
    TicToc();
    ~TicToc();

    void start();
    void stop();
    double get_secs(bool current = false);
    double get_millisecs(bool current = false);
    double get_microsecs(bool current = false);
    double get_nanosecs(bool current = false);

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_;
};

#endif /* SRC_CONSENSUS_TicToc_H_ */
