/*
 * sequence_file_parser_base.h
 *
 *  Created on: Apr 27, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_BASE_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_BASE_H_

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <index/sequence.h>

namespace mindex {

class SequenceFileParserBase;

using SequenceFileParserBasePtr = std::unique_ptr<SequenceFileParserBase>;

using HeaderGroupType = std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::string>>>;

class SequenceFileParserBase {
public:
   // SequenceFileParserBase(const std::string& path);
   virtual ~SequenceFileParserBase() {}

   virtual bool Open(const std::string& in_path) = 0;
   virtual SequencePtr YieldSequence() = 0;
   virtual std::string GetFileHeaderAsString() const = 0;
   virtual int64_t GetFileOffset() const = 0;
   virtual std::string GetFilePath() const = 0;
   virtual bool FileSeek(int64_t abs_pos) = 0;
   virtual bool IsOpen() const = 0;
   virtual const HeaderGroupType& GetHeaderGroups() const = 0;
};

}

#endif
