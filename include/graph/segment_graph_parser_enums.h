/*
 * segment_graph_parser_enums_enums.h
 *
 *  Created on: Dec 23, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_GRAPH_SEGMENT_GRAPH_PARSER_ENUMS_H_
#define SRC_GRAPH_SEGMENT_GRAPH_PARSER_ENUMS_H_

namespace raptor {

enum class GraphFormat { GFA, BED, Unknown };

inline GraphFormat GraphFormatFromString(const std::string& format_str) {
    GraphFormat ret;
    if (format_str == "gfa") {
        ret = GraphFormat::GFA;
    } else if (format_str == "bed") {
        ret = GraphFormat::BED;
    } else {
        ret = GraphFormat::Unknown;
    }
    return ret;
}

}

#endif
