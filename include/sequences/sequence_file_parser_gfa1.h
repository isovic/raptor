/*
 * sequence_file_parser_gfa1.h
 *
 *  Created on: Apr 28, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_GFA1_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_GFA1_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <index/sequence.h>

#include <zlib.h>

namespace mindex {

class SequenceFileParserGfa1;

typedef std::unique_ptr<mindex::SequenceFileParserGfa1> SequenceFileParserGfa1Ptr;

mindex::SequenceFileParserGfa1Ptr createSequenceFileParserGfa1(const std::string& in_path);

class SequenceFileParserGfa1 : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserGfa1Ptr createSequenceFileParserGfa1(const std::string& in_path);
   ~SequenceFileParserGfa1();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;

private:
   SequenceFileParserGfa1();

   std::string path_;
   std::string header_;

   gzFile fp_gzip_;

};

}

#endif
