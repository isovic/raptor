/*
 * TicToc.h
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSENSUS_TicToc2_H_
#define SRC_CONSENSUS_TicToc2_H_

#include <ctime>
#include <chrono>

class TicToc {
  public:
    TicToc();
    ~TicToc();

    void start();
    void stop();
    double get_secs(bool current = false) const;
    double get_millisecs(bool current = false) const;
    double get_microsecs(bool current = false) const;
    double get_nanosecs(bool current = false) const;

    double get_cpu_secs(bool current = false) const;
    double get_cpu_millisecs(bool current = false) const;
    double get_cpu_microsecs(bool current = false) const;
    double get_cpu_nanosecs(bool current = false) const;

  private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    std::chrono::time_point<std::chrono::high_resolution_clock> end_;
    std::clock_t start_cpu_;
    std::clock_t end_cpu_;

    inline double get_cpu_duration_(bool current, double factor) const;
};

#endif /* SRC_CONSENSUS_TicToc_H_ */
