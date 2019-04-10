/*
 * alignment_options.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNMENT_OPTIONS_H_
#define SRC_ALIGNER_ALIGNMENT_OPTIONS_H_

#include <stdint.h>
#include <aligner/global_margins.h>
#include <aligner/pairwise_penalties.h>

namespace raptor {

class AlignmentOptions {
 public:
  AlignmentOptions(const raptor::PiecewisePenalties& _p, int32_t _bandwidth,
                   bool _do_traceback, int32_t _zbandwidth, int32_t _zdrop, int32_t _end_bonus)
                    : p(_p), bandwidth(_bandwidth),
                      do_traceback(_do_traceback),
                      zbandwidth(_zbandwidth),
                      zdrop(_zdrop),
                      end_bonus(_end_bonus),
                      gm(),
                      stop_ext_on_zero_score(true) {
  }

  AlignmentOptions() :  p(), bandwidth(-1),
                        do_traceback(true),
                        zbandwidth(500),
                        zdrop(-1),
                        end_bonus(0),
                        gm(),
                        stop_ext_on_zero_score(true) {
  }

  raptor::PiecewisePenalties p; // Alignment scores/penalties.
  int32_t bandwidth;                // Band for banded alignment. If < 0, banded alignment is turned off.
  bool do_traceback;        // If traceback is not needed, then there is no need to alocate a large
                            // matrix to store directions.

  int32_t zbandwidth;       // Bandwidth for extension alignment.
  int32_t zdrop;            // Z-drop heuristic, applied on extension alignment.

  int32_t end_bonus;

  GlobalMargins gm;

  bool stop_ext_on_zero_score;  // Emulates local alignment if true. Enabling this might prevent end-to-end extensions. Disabling it can severely reduce alignment score / identity.
};

}

#endif
