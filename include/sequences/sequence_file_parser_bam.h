/*
 * sequence_file_parser_bam.h
 *
 *  Created on: May 03, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_PARSER_BAM_H_
#define SRC_INDEX_SEQUENCE_FILE_PARSER_BAM_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_parser_base.h>
#include <index/sequence.h>

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

private:
   SequenceFileParserBam();

   void ParseHeader_();

   std::string path_;
   std::unique_ptr<PacBio::BAM::BamReader> bam_reader_;
   std::unique_ptr<PacBio::BAM::BamFile> bam_file_;
   std::unique_ptr<PacBio::BAM::DataSet> dataset_;
   std::string sam_header_;
};

}

#endif
