/*
 *  stringutil.cc
 *
 *  Created on: Mar 12, 2019
 *      Author: Ivan Sovic
 */

#include <utility/fofn.h>
#include <algorithm>
#include <fstream>

namespace raptor {

std::vector<std::string> ParseFOFN(const std::string& in_path) {
    std::vector<std::string> ret;
    std::ifstream ifs(in_path);
    if (ifs.is_open() == false) {
        return ret;
    }
    std::string line;
    while (getline(ifs, line)) {
        ret.emplace_back(line);
    }
    return ret;
}

}
