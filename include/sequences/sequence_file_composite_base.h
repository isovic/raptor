/*
 * sequence_file_composite_base.h
 *
 *  Created on: June 1, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_COMPOSITE_BASE_H_
#define SRC_INDEX_SEQUENCE_FILE_COMPOSITE_BASE_H_

#include <memory>
#include <string>
#include <vector>

#include <sequences/sequence.h>
#include <sequences/sequence_file.h>

namespace mindex {

class SequenceFileCompositeBase;

using SequenceFileCompositeBasePtr = std::unique_ptr<SequenceFileCompositeBase>;

class SequenceFileCompositeBase {
public:

   // // Allow opening one or more files, with manual format specification for all.
   // virtual bool Open(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts) = 0;

   // SequenceFileParserBase(const std::string& path);
   virtual ~SequenceFileCompositeBase() {}

   // Returns a vector of all file headers in the order of input files.
   virtual std::vector<std::string> GetFileHeaders() const = 0;

   // Currently open file path.
   virtual std::string GetCurrentFilePath() const = 0;

   // Internal ordinal ID of the currently open file.
   virtual int64_t GetCurrentFileId() const = 0;

   // Current offset in the currently open file.
   virtual int64_t GetCurrentFileOffset() const = 0;

   // The previous offset in the same current open file.
   virtual int64_t GetCurrentFilePreviousOffset() const = 0;

   // The file_id is one of the input files.
   virtual bool FileSeek(int64_t file_id, int64_t abs_pos) = 0;

   // The file_name must be one of the input files.
   virtual bool FileSeek(const std::string& file_name, int64_t abs_pos) = 0;

   // Return a single sequence from the file, the next one in line.
   virtual mindex::SequencePtr YieldSequence() = 0;

   // Also loads a single sequence, but it returns it as a SequenceFile.
   virtual mindex::SequenceFilePtr YieldBatchOfOne() = 0;

   // Loads all sequences from the input files.
   virtual mindex::SequenceFilePtr YieldBatchMB(int64_t batch_size_mb) = 0;

   // Loads all sequences from the input files.
   virtual mindex::SequenceFilePtr YieldAll() = 0;

   // Merges all header groups.
   virtual const HeaderGroupType& GetHeaderGroups() const = 0;

private:

   // // Allow opening one or more files, with auto format detection.
   // virtual bool Open_(const std::vector<std::string>& in_paths) = 0;

   // // Allow opening one or more files, where all should be in this specified format.
   // virtual bool Open_(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt) = 0;

   // // Allow opening one or more files, with manual format specification for all.
   // virtual bool Open_(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts) = 0;



   //  bool LoadAll(const std::vector<std::string>& headers, const std::vector<std::string>& seqs, bool convert_to_uppercase=true);
   //  bool LoadAll(bool convert_to_uppercase=true);
   //  bool LoadBatchMB(int64_t batch_size_mb, bool convert_to_uppercase=true);
   //  bool LoadBatchOfOne(bool convert_to_uppercase=true);

   //  mindex::SequencePtr FetchSequence(const std::string& qname);
   //  mindex::SequencePtr FetchSequence(int64_t qid);

   //  mindex::SequenceFilePtr FetchBlock(const int64_t block_id);
   //  mindex::SequenceFilePtr FetchAll();

};

}

#endif
