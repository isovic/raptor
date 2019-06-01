/*
 * sequence_file_parser_bam.h
 *
 *  Created on: May 03, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_BAM_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_BAM_H_

#ifdef RAPTOR_COMPILED_WITH_PBBAM

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <sequences/sequence_file_parser_base.h>
#include <sequences/sequence.h>

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/QualityValues.h>

namespace mindex {

class SequenceFileParserBam;

// typedef std::unique_ptr<mindex::SequenceFileParserBam> SequenceFileParserBamPtr;

mindex::SequenceFileParserBasePtr createSequenceFileParserBam(const std::string& in_path);

class SequenceFileParserBam : SequenceFileParserBase {
public:
   friend mindex::SequenceFileParserBasePtr createSequenceFileParserBam(const std::string& in_path);
   ~SequenceFileParserBam();

   bool Open(const std::string& path);

   SequencePtr YieldSequence();
   std::string GetFileHeaderAsString() const;
   int64_t GetFileOffset() const;
   std::string GetFilePath() const;
   bool FileSeek(int64_t pos);
   bool IsOpen() const;
   const HeaderGroupType& GetHeaderGroups() const;

private:
   SequenceFileParserBam();

   void ParseHeader_();
   void ParseReadGroupAndProgramGroupFromHeader_(const std::string& header);

   std::string path_;
   std::unique_ptr<PacBio::BAM::BamReader> bam_reader_;
   std::unique_ptr<PacBio::BAM::BamFile> bam_file_;
   std::unique_ptr<PacBio::BAM::DataSet> dataset_;
   std::string sam_header_;
   // The following is complicated, but it translates to:
   //   header_group_[group_name][ID] = map_of_tags;
   // For example:
   //   header_group_["RG"]["995ef077"] = {"ID": "995ef077", "PL": "PACBIO"};
   HeaderGroupType header_groups_;
};

}

#endif

#endif
