/*
 * sequence_file_parser_fastx.h
 *
 *  Created on: Apr 27, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_FASTX_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_FASTX_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <index/sequence.h>

#include <zlib.h>
#include <sequences/kseq.h>

KSEQ_INIT(gzFile, gzread)

namespace mindex {

class SequenceFileParserFastx;

typedef std::unique_ptr<mindex::SequenceFileParserFastx> SequenceFileParserFastxPtr;

mindex::SequenceFileParserFastxPtr createSequenceFileParserFastx(const std::string& in_path);

class SequenceFileParserFastx : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserFastxPtr createSequenceFileParserFastx(const std::string& in_path);
   ~SequenceFileParserFastx();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;

private:
   SequenceFileParserFastx();

   std::string path_;

   kseq_t *fp_kseq_;
   gzFile fp_gzip_;

};

}

#endif
