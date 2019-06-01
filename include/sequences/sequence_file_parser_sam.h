/*
 * sequence_file_parser_sam.h
 *
 *  Created on: Apr 28, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_SAM_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_SAM_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <sequences/sequence.h>

#include <zlib.h>

namespace mindex {

class SequenceFileParserSam;

// typedef std::unique_ptr<mindex::SequenceFileParserSam> SequenceFileParserSamPtr;

mindex::SequenceFileParserBasePtr createSequenceFileParserSam(const std::string& in_path);

class SequenceFileParserSam : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserBasePtr createSequenceFileParserSam(const std::string& in_path);
   ~SequenceFileParserSam();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;
   const HeaderGroupType& GetHeaderGroups() const;

private:
   SequenceFileParserSam();

   void ParseHeader_();
   void ParseReadGroupAndProgramGroupFromHeader_(const std::string& header);

   std::string path_;
   std::string header_;

   gzFile fp_gzip_;

   // The following is complicated, but it translates to:
   //   header_group_[group_name][ID] = map_of_tags;
   // For example:
   //   header_group_["RG"]["995ef077"] = {"ID": "995ef077", "PL": "PACBIO"};
   HeaderGroupType header_groups_;
};

}

#endif
