/*
 * global_margin.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_GLOBAL_MARGINS_H_
#define SRC_ALIGNER_GLOBAL_MARGINS_H_

#include <stdint.h>

namespace raptor {

// If any global margin is true, then the corresponding will be penalized.
// Concretely, if top/left are true, then the first row/column will be initialized
// to the multiple of the gap extend penalty in global alignment.
// If bottom is false, the maximum of last row will be found instead of taking
// the bottom right corner for global alignment.
// If right is false, the maximum of last column will be found instead of taking
// the bottom right corner for global alignment.
class GlobalMargins {
   public:
    GlobalMargins() : top(true), left(true), bottom(true), right(true) {}
    GlobalMargins(bool _top, bool _left, bool _bottom, bool _right)
        : top(_top), left(_left), bottom(_bottom), right(_right) {}
    bool top, left, bottom, right;
};

}  // namespace raptor

#endif
