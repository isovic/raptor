/*
 * sequence_file_parser_pb_xml.h
 *
 *  Created on: May 29, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_PB_XML_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_PB_XML_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <sequences/sequence.h>

#include <pbbam/DataSet.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/QualityValues.h>
#include <pbbam/PbiFilterQuery.h>

namespace mindex {

class SequenceFileParserPbXml;

mindex::SequenceFileParserBasePtr createSequenceFileParserPbXml(const std::string& in_path);

class SequenceFileParserPbXml : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserBasePtr createSequenceFileParserPbXml(const std::string& in_path);
   ~SequenceFileParserPbXml();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;
   const HeaderGroupType& GetHeaderGroups() const;

private:
   SequenceFileParserPbXml();

   void ParseHeader_();
   void ParseReadGroupAndProgramGroupFromHeader_(const std::string& header);

   std::string path_;
   std::unique_ptr<PacBio::BAM::BamFile> bam_file_;
   std::unique_ptr<PacBio::BAM::DataSet> dataset_;
   std::string sam_header_;

   HeaderGroupType header_groups_;

   std::vector<std::string> files_;                            // List of files loaded from the Dataset.
   std::vector<std::vector<int64_t>> file_offsets_;            // For each file, list of all sequence offsets.
   int64_t open_file_id_;
   int64_t open_file_offset_id_;

   mindex::SequenceFileParserBasePtr bam_parser_;
};

}

#endif
