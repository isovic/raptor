/*
 * sequence_file_parser_factory.
 *
 *  Created on: Apr 28, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_FACTORY_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_FACTORY_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <sequences/sequence_file_parser_fastx.h>
#include <sequences/sequence_file_parser_gfa1.h>
#include <sequences/sequence_file_parser_gfa2.h>
#include <sequences/sequence_file_parser_sam.h>
#include <sequences/sequence_file_parser_bam.h>
#include <sequences/sequence_file_parser_pbxml.h>
#include <index/sequence_file_enums.h>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParser(const std::string& in_path, mindex::SequenceFormat in_fmt);

}

#endif
