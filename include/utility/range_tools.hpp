/*
 *  range_tools.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Ivan Sovic
 */

#include <algorithm>
#include <functional>
#include <tuple>
#include <vector>

#ifndef SRC_ISTL_RANGE_TOOLS_HPP_
#define SRC_ISTL_RANGE_TOOLS_HPP_

namespace istl {

/* Returns a vector of pair numbers: [start, end> indexes within the data vector,
 * which demarcate the range satisfied by the comparison function.
 * The comparison function should return true if the two elements are supposed to be
 * groupped together, otherwise false.
 * The first template class Q is the atomic data type.
 * The second template class T is an iterable container of atomic types, std::vector by default.
*/
template<class Q, class T = std::vector<Q>>
std::vector<std::pair<size_t, size_t>> FindRanges(const T& data, std::function<bool(const Q& a,
                                      const Q& b)> comp_eq =
                                      [](const Q& a, const Q& b)
                                      { return a == b; } ) {

    std::vector<std::pair<size_t, size_t>> ret;

    if (data.size() == 0) {
        return ret;
    }

    size_t iter = 0;
    size_t prev_end = 0;

    for (auto& data_point: data) {
        ++iter;
        size_t id = iter - 1;
        // Skip the first one.
        if (id == 0) {
            continue;
        }
        // If the same as before, continue crunching.
        if (comp_eq(data[id], data[id - 1])) {
            continue;
        }
        ret.emplace_back(std::make_pair(prev_end, id));
        prev_end = id;
    }

    if (prev_end != data.size()) {
        ret.emplace_back(std::make_pair(prev_end, data.size()));
    }

    return ret;
}

}

#endif

