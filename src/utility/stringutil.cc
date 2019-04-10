/*
 *  stringutil.cc
 *
 *  Created on: Oct 11, 2014
 *      Author: ivan
 */

#include <utility/stringutil.h>
#include <algorithm>

namespace raptor {

std::vector<std::string> TokenizeToWhitespaces(const std::string& str) {
    std::vector<std::string> words;
    std::istringstream ss(str);
    std::string word;
    while (ss >> word) {
        words.emplace_back(word);
    }
    return words;
}

std::vector<std::string> Tokenize(const std::string& str, const char delimiter) {
    std::vector<std::string> words;
    std::istringstream ss(str);
    std::string line;
    while (std::getline(ss, line, delimiter)) {
        if (line.size() == 0) {
            continue;
        }
        words.push_back(line);
    }
    return words;
}

std::string TrimToFirstDelimiter(const std::string& original_string, char delimiter) {
    std::string::size_type loc = original_string.find(delimiter, 0);
    if (loc != std::string::npos) {
        return original_string.substr(0, loc);

    } else {
        // There is no spaces in the string, do nothing and just report it as is.
        return original_string;
    }
    return std::string("");
}

std::string TrimToFirstSpace(const std::string& original_string) {
    return TrimToFirstDelimiter(original_string, ' ');
}

std::string TrimToFirstWhiteSpace(const std::string& original_string) {
    auto tokens = TokenizeToWhitespaces(original_string);
    if (tokens.empty()) {
        return std::string("");
    }
    return tokens[0];
}

}
