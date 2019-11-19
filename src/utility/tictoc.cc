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
    start_cpu_ = std::clock();
}

void TicToc::stop() {
    end_ = std::chrono::high_resolution_clock::now();
    end_cpu_ = std::clock();
}

double TicToc::get_secs(bool current) const {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_millisecs(bool current) const {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_microsecs(bool current) const {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_nanosecs(bool current) const {
    auto end = (current) ? std::chrono::high_resolution_clock::now() : end_;
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_).count();
    return elapsed;
}

double TicToc::get_cpu_duration_(bool current, double factor) const {
    auto end_cpu = (current) ? std::clock() : end_cpu_;
    double elapsed = factor * static_cast<double>(end_cpu - start_cpu_) / CLOCKS_PER_SEC;
    return elapsed;
}

double TicToc::get_cpu_secs(bool current) const {
    return get_cpu_duration_(current, 1.0);
}

double TicToc::get_cpu_millisecs(bool current) const {
    return get_cpu_duration_(current, 1000.0);
}

double TicToc::get_cpu_microsecs(bool current) const {
    return get_cpu_duration_(current, 1000000.0);
}

double TicToc::get_cpu_nanosecs(bool current) const {
    return get_cpu_duration_(current, 1000000000.0);
}
