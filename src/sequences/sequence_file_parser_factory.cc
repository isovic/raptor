/*
 * sequence_file_parser_factory.cc
 *
 *  Created on: Apr 28, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_factory.h>
#include <index/sequence_file_enums.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParser(const std::string& in_path, mindex::SequenceFormat in_fmt) {

    if (in_fmt == mindex::SequenceFormat::Auto) {
        auto ext = raptor::GetFileExt(in_path);
        in_fmt = SequenceFormatFromString(ext);
    }

    mindex::SequenceFileParserBasePtr ret = nullptr;

    if (in_fmt == mindex::SequenceFormat::Fasta) {
        ret = mindex::createSequenceFileParserFastx(in_path);
    } else if (in_fmt == mindex::SequenceFormat::Fastq) {
        ret = mindex::createSequenceFileParserFastx(in_path);
    } else if (in_fmt == mindex::SequenceFormat::SAM) {
        ret = mindex::createSequenceFileParserSam(in_path);
    } else if (in_fmt == mindex::SequenceFormat::GFA1) {
        ret = mindex::createSequenceFileParserGfa1(in_path);
    } else if (in_fmt == mindex::SequenceFormat::GFA2) {
        ret = mindex::createSequenceFileParserGfa2(in_path);
    } else {
        ret = nullptr;
    }

    return ret;
}

}
