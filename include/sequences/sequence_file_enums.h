/*
 * sequence_file_enums.h
 *
 *  Created on: Dec 23, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_PARAMS_SEQUENCE_FILE_ENUMS_H_
#define SRC_PARAMS_SEQUENCE_FILE_ENUMS_H_

#include <algorithm>
#include <cstdlib>
#include <utility/files.hpp>

namespace mindex {

enum class BatchLoadType { MB, Coverage };
enum class SequenceFormat {
            Auto, Fasta, Fastq, SAM,
            #ifdef RAPTOR_COMPILED_WITH_PBBAM
                        BAM, XML,
            #endif
            GFA, GFA1, GFA2, RaptorDB, FOFN, Unknown };

inline SequenceFormat SequenceFormatFromString(const std::string& format_str) {
    SequenceFormat ret;
    if (format_str == "auto") {
        ret = SequenceFormat::Auto;
    } else if (format_str == "fasta" || format_str == "fa") {
        ret = SequenceFormat::Fasta;
    } else if (format_str == "fastq" || format_str == "fq") {
        ret = SequenceFormat::Fastq;
    } else if (format_str == "sam") {
        ret = SequenceFormat::SAM;
#ifdef RAPTOR_COMPILED_WITH_PBBAM
    } else if (format_str == "bam") {
        ret = SequenceFormat::BAM;
    } else if (format_str == "xml") {
        ret = SequenceFormat::XML;
#endif
    } else if (format_str == "gfa") {
        ret = SequenceFormat::GFA;
    } else if (format_str == "gfa1") {
        ret = SequenceFormat::GFA1;
    } else if (format_str == "gfa2") {
        ret = SequenceFormat::GFA2;
    } else if (format_str == "rdb") {
        ret = SequenceFormat::RaptorDB;
    } else if (format_str == "fofn") {
        ret = SequenceFormat::FOFN;
    } else {
        ret = SequenceFormat::Unknown;
    }
    return ret;
}

inline std::string SequenceFormatToString(const SequenceFormat& fmt) {
    std::string ret("unknown");

    switch(fmt) {
        case SequenceFormat::Auto:
            ret = "auto";
            break;
        case SequenceFormat::Fasta:
            ret = "fasta";
            break;
        case SequenceFormat::Fastq:
            ret = "fastq";
            break;
        case SequenceFormat::SAM:
            ret = "sam";
            break;
#ifdef RAPTOR_COMPILED_WITH_PBBAM
        case SequenceFormat::BAM:
            ret = "bam";
            break;
        case SequenceFormat::XML:
            ret = "xml";
            break;
#endif
        case SequenceFormat::GFA:
            ret = "gfa";
            break;
        case SequenceFormat::GFA1:
            ret = "gfa1";
            break;
        case SequenceFormat::GFA2:
            ret = "gfa2";
            break;
        case SequenceFormat::RaptorDB:
            ret = "rdb";
            break;
        case SequenceFormat::FOFN:
            ret = "fofn";
            break;
        default:
            ret = "unknown";
    }
    return ret;
}

/*
 * Deduce the format from the file extension.
*/
inline SequenceFormat GetSequenceFormatFromPath(const std::string& path) {
    std::string ext = raptor::GetFileExtWithoutGZ(path);
    return SequenceFormatFromString(ext);
}

/*
 * In case the input format has been specified explicitly before (i.e. it is different from Auto),
 * then return that one.
 * Otherwise, deduce the format from the path.
*/
inline mindex::SequenceFormat GetSequenceFormatFromPath(const std::string& path, const mindex::SequenceFormat& apriori_in_fmt) {
    if (apriori_in_fmt != mindex::SequenceFormat::Auto) {
        return apriori_in_fmt;
    }
    mindex::SequenceFormat fmt = mindex::GetSequenceFormatFromPath(path);
    return fmt;
}

}

#endif
