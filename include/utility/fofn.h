/*
 * fofn.h
 *
 *  Created on: Mar 12, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTILITY_FOFN_HPP_
#define SRC_UTILITY_FOFN_HPP_

#include <string>
#include <vector>

namespace raptor {

std::vector<std::string> ParseFOFN(const std::string& in_path);

}  // namespace raptor

#endif
