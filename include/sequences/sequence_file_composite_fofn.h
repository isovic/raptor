/*
 * sequence_file_composite_fofn.h
 *
 *  Created on: June 1, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_COMPOSITE_FOFN_H_
#define SRC_INDEX_SEQUENCE_FILE_COMPOSITE_FOFN_H_

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include <sequences/sequence_file_composite_base.h>
#include <sequences/sequence_file_parser_utils.h>

namespace mindex {

class SequenceFileCompositeFofn;

using SequenceFileCompositeFofnPtr = std::unique_ptr<SequenceFileCompositeFofn>;

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::string& fofn_path);
mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt);
mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

class SequenceFileCompositeFofn : SequenceFileCompositeBase {
public:
   friend mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::string& fofn_path);
   friend mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt);
   friend mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

   // // Allow opening one or more files, with manual format specification for all.
   // bool Open(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

   ~SequenceFileCompositeFofn();

   // Returns a vector of all file headers in the order of input files.
   std::vector<std::string> GetFileHeaders() const;

   // Currently open file path.
   std::string GetCurrentFilePath() const;

   // Internal ordinal ID of the currently open file.
   int64_t GetCurrentFileId() const;

   // Current offset in the currently open file.
   int64_t GetCurrentFileOffset() const;

   // The previous offset in the same current open file.
   int64_t GetCurrentFilePreviousOffset() const;

   // The file_id is one of the input files.
   bool FileSeek(int64_t file_id, int64_t abs_pos);

   // The file_name must be one of the input files.
   bool FileSeek(const std::string& file_name, int64_t abs_pos);

   // Return a single sequence from the file, the next one in line.
   mindex::SequencePtr YieldSequence();

   // Also loads a single sequence, but it returns it as a SequenceFile.
   mindex::SequenceFilePtr YieldBatchOfOne();

   // Loads all sequences from the input files.
   mindex::SequenceFilePtr YieldBatchMB(int64_t batch_size_mb);

   // Loads all sequences from the input files.
   mindex::SequenceFilePtr YieldAll();

   // Merges all header groups.
   const HeaderGroupType& GetHeaderGroups() const;

 private:
   SequenceFileCompositeFofn();
   SequenceFileCompositeFofn(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);
   SequenceFileCompositeFofn(const SequenceFileCompositeFofn&) = delete;
   SequenceFileCompositeFofn& operator=(const SequenceFileCompositeFofn&) = delete;

   bool Open_(const std::string& in_path, mindex::SequenceFormat in_fmt);

   mindex::SequenceFormat GetOpenFileFormat_() const;
   int64_t GetOpenFileTell_() const;
   bool ParseHeaders_();

   std::vector<std::string> files_;
   std::vector<mindex::SequenceFormat> formats_;
   int64_t curr_open_file_;
   mindex::SequenceFormat curr_input_fmt_;               // If the input format was auto, this is the resolved version.
   mindex::SequenceFileParserBasePtr parser_;
   bool convert_to_uppercase_;
   int64_t open_file_seq_offset_prev_;
   int64_t open_file_seq_offset_curr_;
   std::vector<std::string> headers_;
   HeaderGroupType header_groups_;
   int64_t batch_id_;
   int64_t num_loaded_seqs_;
};

}

#endif
