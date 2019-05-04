/*
 * sequence.h
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_FILE2_H_
#define SRC_SEQUENCES_SEQUENCE_FILE2_H_

#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <index/sequence.h>
#include <params/params_raptor.h>
#include <index/sequence_file_enums.h>
#include <sequences/sequence_file_parser_factory.h>

namespace mindex {

class SequenceFile;

typedef std::shared_ptr<mindex::SequenceFile> SequenceFilePtr;

mindex::SequenceFilePtr createSequenceFile();
mindex::SequenceFilePtr createSequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt);
mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

class SequenceFile {
public:
    /*
     * Creation.
    */
    friend mindex::SequenceFilePtr createSequenceFile();
    friend mindex::SequenceFilePtr createSequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt);
    friend mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

    bool LoadAll(const std::vector<std::string>& headers, const std::vector<std::string>& seqs, bool convert_to_uppercase=true);
    bool LoadAll(bool convert_to_uppercase=true);
    bool LoadBatchMB(int64_t batch_size_mb, bool convert_to_uppercase=true);
    bool LoadBatchOfOne(bool convert_to_uppercase=true);

    const mindex::SequencePtr& GetSeqByID(int64_t id) const;
    const mindex::SequencePtr& GetSeqByAbsID(int64_t abs_id) const;
    void Add(std::unique_ptr<mindex::Sequence> seq);

    bool Serialize(std::ostream& ofs, mindex::SequenceFormat to_fmt);

    std::string GetOpenFilePath() const;
    mindex::SequenceFormat GetOpenFileFormat() const;
    int64_t GetOpenFileTell() const;

    /*
     * Getters.
    */
    const std::vector<mindex::SequencePtr>& seqs() const {
        return seqs_;
    }
    int64_t batch_id() const {
        return batch_id_;
    }
    int64_t batch_start_seq_id() const {
        return batch_start_seq_id_;
    }
    int64_t total_size() const {
        return total_size_;
    }
    int64_t size() const {
        return seqs_.size();
    }
    const mindex::SequencePtr& operator[](int64_t id) const {
        return GetSeqByID(id);
    }
    int64_t GetOpenFileId() const {
        return curr_open_file_;
    }

    /*
     * Setters.
    */
    void seqs(std::vector<mindex::SequencePtr>& _seqs) {
        seqs_ = std::move(_seqs);
    }
    void batch_start_seq_id(int64_t val) {
        batch_start_seq_id_ = val;
    }

    /*
     * Modifiers.
    */
    std::vector<mindex::SequencePtr>& seqs() {
        return seqs_;
    }

private:
    SequenceFile();
    SequenceFile(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);
    SequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt);
    SequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

    SequenceFile(const SequenceFile&) = delete;
    SequenceFile& operator=(const SequenceFile&) = delete;

    bool Open_(const std::string& in_path, mindex::SequenceFormat in_fmt);
    mindex::SequencePtr YieldSequence_(bool convert_to_uppercase);



    std::vector<std::pair<std::string, mindex::SequenceFormat>> in_files_;
    int64_t curr_open_file_;
    mindex::SequenceFormat curr_input_fmt_;
    mindex::SequenceFileParserBasePtr parser_;

    std::vector<mindex::SequencePtr> seqs_;     // All sequences.
    int64_t batch_id_;              // The ID of the currently loaded batch.
    int64_t batch_start_seq_id_;    // Absolute ID of the first sequence in the loaded batch.
    int64_t total_size_;            // Sum of all sequence lengths.
    mindex::SequencePtr dummy_nullptr_seq_;
};

}

#endif
