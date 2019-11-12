/*
 * sequence_file_enums.h
 *
 *  Created on: Nov 12, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_PARAMS_SEQUENCE_FILE_UTILS_H_
#define SRC_PARAMS_SEQUENCE_FILE_UTILS_H_

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>
#include <sequences/sequence_file_enums.h>
#include <utility/files.hpp>
#include <utility/fofn.h>

namespace mindex {

std::vector<std::string> ExpandPathList(
                const mindex::SequenceFormat& apriori_in_fmt,
                const std::string& apriori_in_fmt_str,
                const std::vector<std::string>& in_paths);

bool ValidateInputFiles(
                    const mindex::SequenceFormat& apriori_in_fmt,
                    const std::vector<std::string>& paths);

bool IsInputFormatRaptorDB(const mindex::SequenceFormat& apriori_in_fmt,
                            const std::vector<std::string>& paths);

}

#endif
