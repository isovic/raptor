/*
 * sequence_file_enums.cc
 *
 *  Created on: Nov 12, 2019
 *      Author: Ivan Sovic
 */

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
                const std::vector<std::string>& in_paths) {

    // In case the input contained FOFNs, collect all paths.
    std::vector<std::string> expanded_in_paths;

    for (size_t in_id = 0; in_id < in_paths.size(); ++in_id) {
        const auto& in_path = in_paths[in_id];
        mindex::SequenceFormat fmt = GetSequenceFormatFromPath(in_path, apriori_in_fmt);

        if (fmt == mindex::SequenceFormat::FOFN) {
            // This section parses the FOFN, and appends files to the list.
            std::vector<std::string> fofn_list = raptor::ParseFOFN(in_path);
            expanded_in_paths.insert(expanded_in_paths.end(), fofn_list.begin(), fofn_list.end());

        } else {
            // If the file wasn't a FOFN, just add it to the list.
            expanded_in_paths.emplace_back(in_path);
        }
    }

    return expanded_in_paths;
}

bool ValidateInputFiles(
                    const mindex::SequenceFormat& apriori_in_fmt,
                    const std::vector<std::string>& paths) {

    if (apriori_in_fmt == mindex::SequenceFormat::Unknown) {
        return false;
    }

    int32_t num_rdb_paths = 0;

    for (size_t in_id = 0; in_id < paths.size(); ++in_id) {
        const auto& in_path = paths[in_id];
        std::string in_path_ext = raptor::GetFileExtWithoutGZ(in_path);
        mindex::SequenceFormat curr_path_fmt = mindex::SequenceFormatFromString(in_path_ext);
        mindex::SequenceFormat fmt = (apriori_in_fmt == mindex::SequenceFormat::Auto) ? curr_path_fmt : apriori_in_fmt;

        if (fmt == mindex::SequenceFormat::Unknown) {
            // Exit if an unknown format is found.
            fprintf(stderr, "Unsupported input format '%s' for file '%s'!\n\n", in_path_ext.c_str(), in_path.c_str());
            return false;

        } else if (fmt == mindex::SequenceFormat::RaptorDB) {
            ++num_rdb_paths;
        }

        if (!raptor::FileExists(in_path.c_str())) {
            fprintf(stderr, "File does not exist: '%s'\n\n", in_path.c_str());
            return false;
        }
    }

    // Allow only one RaptorDB file to be loaded.
    if (num_rdb_paths > 1) {
        fprintf(stderr, "Only one RaptorDB sequence file can be loaded.\n\n");
        return false;
    } else if (num_rdb_paths == 1 && paths.size() > 1) {
        fprintf(stderr, "A RaptorDB input cannot be specified along side to another input file.\n\n");
        return false;
    }

    return true;
}

bool IsInputFormatRaptorDB(const mindex::SequenceFormat& apriori_in_fmt,
                            const std::vector<std::string>& paths) {
    if (paths.size() != 1) {
        return false;
    }
    mindex::SequenceFormat fmt = GetSequenceFormatFromPath(paths[0], apriori_in_fmt);
    return fmt == mindex::SequenceFormat::RaptorDB;
}


}
