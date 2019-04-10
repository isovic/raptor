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

static inline int32_t CalcMapqFromScoreCounts(const std::vector<int32_t>& score_counts) {
    int32_t mapq = 0;

    if (score_counts.size() == 0) {
        return mapq;
    }

    // fprintf (stderr, "Score distribution:\n");

    // for (size_t i = 0; i < score_counts.size(); i++) {
    //     fprintf (stderr, "  [%02d] count: %d\n", i, score_counts[i]);
    // }

    int32_t first_non_unique = 0;
    for (int32_t i = 0; i < (int32_t) score_counts.size(); i++) {
        first_non_unique = i;
        if (score_counts[i] > 1) {
            break;
        }
    }

    int32_t bucket_span = (int32_t) (60.0 / (((float) (score_counts.size() - 1) / 2) + 1));
    int32_t mapq_bucket = bucket_span * (first_non_unique / 2); // (first_non_unique - 1) / 2;
    int32_t count = score_counts[first_non_unique];

    // All buckets have count 1.
    // if (count == 1) {
    //     mapq = 60;
    //     return mapq;
    // }

    // printf ("first_non_unique = %d, mapq_bucket = %d, count = %d, raptor::CalcMapq(count) = %d\n", first_non_unique, mapq_bucket, count, raptor::CalcMapq(count));

    mapq = mapq_bucket + raptor::CalcMapq(count);

    mapq = std::min(mapq, 60);
    mapq = std::max(0, mapq);

    // if ((first_non_unique + 1) == score_counts.size()) {
    //     mapq = 60;
    // } else {
    // }

    return mapq;
}

}

#endif
