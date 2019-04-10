/*
 * output_formatter_enums.h
 *
 *  Created on: Dec 23, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_WRITER_OUTPUT_FORMATTTER_ENUMS_H_
#define SRC_WRITER_OUTPUT_FORMATTTER_ENUMS_H_

#include <string>

namespace raptor {

enum class OutputFormat { SAM, MHAP, PAF, GFA2, M4, CAGE_OVL, Unknown };

inline OutputFormat OutputFormatFromString(const std::string& format_str) {
    OutputFormat ret;
    if (format_str == "sam") {
        ret = OutputFormat::SAM;
    } else if (format_str == "paf") {
        ret = OutputFormat::PAF;
    } else if (format_str == "mhap") {
        ret = OutputFormat::MHAP;
    } else if (format_str == "m4") {
        ret = OutputFormat::M4;
    } else if (format_str == "gfa") {
        ret = OutputFormat::GFA2;
    } else if (format_str == "gfa2") {
        ret = OutputFormat::GFA2;
    } else if (format_str == "cage-ovl") {
        ret = OutputFormat::CAGE_OVL;
    } else {
        ret = OutputFormat::Unknown;
    }
    return ret;
}

}

#endif
