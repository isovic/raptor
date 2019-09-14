/*
 * sequence_file_parser_utils.cc
 *
 *  Created on: May 29, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_utils.h>
#include <utility/stringutil.h>
#include <log/log_tools.h>

namespace mindex {

HeaderGroupType ParseReadGroupAndProgramGroupFromSAMHeader(const std::string& header) {
    HeaderGroupType header_groups;

    std::istringstream iss(header);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.size() < 3) {
            continue;
        }
        if (line[0] != '@') {
            break;
        }

        std::string field_name = line.substr(1, 2);

        if (field_name == "RG" || field_name == "PG") {
            std::vector<HeaderTag> tags;

            // Skip the first token, that's "@RG" or "@PG".
            auto tokens = raptor::Tokenize(line, '\t');
            bool is_ok = false;
            std::string field_id;
            for (size_t tid = 1; tid < tokens.size(); ++tid) {
                std::string tag_name = tokens[tid].substr(0, 2);
                std::string tag_val = tokens[tid].substr(3);
                tags.emplace_back(HeaderTag{tag_name, tag_val});
                if (tag_name == "ID") {
                    field_id = tag_val;
                    is_ok = true;
                }
            }

            if (is_ok == false) {
                WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Skipping faulty header line: '%s'.", line.c_str());
                continue;
            }

            header_groups[field_name][field_id] = std::move(tags);
        }
    }
    return header_groups;
}

}
