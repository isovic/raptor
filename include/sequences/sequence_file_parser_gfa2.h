/*
 * sequence_file_parser_gfa2.h
 *
 *  Created on: Apr 28, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_GFA2_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_GFA2_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <index/sequence.h>

#include <zlib.h>

namespace mindex {

class SequenceFileParserGfa2;

// typedef std::unique_ptr<mindex::SequenceFileParserGfa2> SequenceFileParserGfa2Ptr;

mindex::SequenceFileParserBasePtr createSequenceFileParserGfa2(const std::string& in_path);

class SequenceFileParserGfa2 : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserBasePtr createSequenceFileParserGfa2(const std::string& in_path);
   ~SequenceFileParserGfa2();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;

private:
   SequenceFileParserGfa2();

   std::string path_;
   std::string header_;

   gzFile fp_gzip_;

};

}

#endif
