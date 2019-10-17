/*
 * tictoc.cc
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */


#include <utility/tictoc.h>

TicToc::TicToc()
    : start_(), end_() {
    start();
    end_ = start_;
}

TicToc::~TicToc() {
}

void TicToc::start() {
    start_ = std::chrono::high_resolution_clock::now();
}

void TicToc::stop() {
    end_ = std::chrono::high_resolution_clock::now();
}

double TicToc::get_secs(bool current) {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_millisecs(bool current) {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_microsecs(bool current) {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_nanosecs(bool current) {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count();
    return elapsed;
}
