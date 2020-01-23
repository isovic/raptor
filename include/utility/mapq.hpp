/*
 * revcmp.hpp
 *
 *  Created on: Dec 9, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTILITY_MAPQ_HPP_
#define SRC_UTILITY_MAPQ_HPP_

#include <cmath>
#include <string>

namespace raptor {

/*
* @brief Calculates the mapping quality for the alignments, based on the number
* possible mapping locations.
*/
static inline int32_t CalcMapq(int64_t num_locations) {
    int32_t ret = (num_locations == 0) ? 255 :
                (num_locations == 1) ? 60:
                (int32_t) round(-10 * log10(1.0 - 1.0 / ((double) num_locations)));
    return ret;
}

}

#endif
