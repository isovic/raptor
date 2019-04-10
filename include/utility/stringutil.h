/*
 * stringutil.h
 *
 *  Created on: Jan 09, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_UTILITY_STRINGUTIL_HPP_
#define SRC_UTILITY_STRINGUTIL_HPP_

#include <string>
#include <sstream>
#include <vector>

namespace raptor {

std::vector<std::string> Tokenize(const std::string& str, const char delimiter);
std::vector<std::string> TokenizeToWhitespaces(const std::string& str);
std::string TrimToFirstDelimiter(const std::string& original_string, char delimiter);
std::string TrimToFirstSpace(const std::string& original_string);
std::string TrimToFirstWhiteSpace(const std::string& original_string);

}  // namespace raptor

#endif
