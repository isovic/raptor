/*
 * revcmp.hpp
 *
 *  Created on: Dec 9, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTILITY_TRIM_HPP_
#define SRC_UTILITY_TRIM_HPP_

#include <string>

namespace raptor {

static std::string TrimToFirstWhitespace(const std::string &original_string) {
  std::string::size_type loc1 = original_string.find(' ', 0);
  std::string::size_type loc2 = original_string.find('\t', 0);

  std::string::size_type loc = loc1;
  if (loc == std::string::npos) {
      loc = loc2;
  }

  if (loc1 != std::string::npos && loc2 != std::string::npos) {
      loc = std::min(loc1, loc2);
  }

  if (loc != std::string::npos) {
    return original_string.substr(0, loc);

  } else {

    // There are no whitespaces in the string, do nothing and just report it as is.
    return original_string;
  }

  return std::string("");
}

}

#endif
