/*
 * sequence_file_parser_utils.h
 *
 *  Created on: May 29, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_UTILS_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_UTILS_H_

#include <string>
#include <vector>
#include <unordered_map>

namespace mindex {

class HeaderTag {
   public:
   std::string name;
   std::string val;
};

using HeaderGroupType = std::unordered_map<std::string, std::unordered_map<std::string, std::vector<HeaderTag>>>;

HeaderGroupType ParseReadGroupAndProgramGroupFromSAMHeader(const std::string& header);

}

#endif
