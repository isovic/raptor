/*
 * tictoc.cc
 *
 *  Created on: Nov 29, 2016
 *      Author: isovic
 */


#include <utility/tictoc.h>

TicToc::TicToc() {
}

TicToc::~TicToc() {
}

void TicToc::start() {
  start_ = clock();
}

void TicToc::stop() {
  end_ = clock();
}

double TicToc::get_secs() {
  double elapsed_secs = double(end_ - start_) / CLOCKS_PER_SEC;
  return elapsed_secs;
}

double TicToc::get_msecs() {
  double elapsed_secs = double(end_ - start_) / CLOCKS_PER_SEC;
  return elapsed_secs * 1000.0;
}

double TicToc::get_secs_current() {
  double elapsed_secs = double(clock() - start_) / CLOCKS_PER_SEC;
  return elapsed_secs;
}
