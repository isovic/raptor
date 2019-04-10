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
#include <index/sequence_file_handlers.h>
#include <index/sequence_file_enums.h>
#include <index/sequence_deserializer.h>

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
//    friend mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);
    friend mindex::SequenceFilePtr createSequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt);
    friend mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

    bool Open(const std::string& in_path, mindex::SequenceFormat in_fmt);
    bool LoadAll(const std::vector<std::string>& headers, const std::vector<std::string>& seqs, bool convert_to_uppercase=true);
    bool LoadAll(bool convert_to_uppercase=true);
    // bool LoadBatchN();
    bool LoadBatchMB(int64_t batch_size_mb, bool convert_to_uppercase=true);
    bool LoadBatchOfOne(bool convert_to_uppercase=true);

    const mindex::SequencePtr& GetSeqByID(int64_t id) const {
        if (id < 0 || id >= seqs_.size()) {
            std::cerr << "[GetSeqByID] Warnning: Requested id is out of scope of current batch. id = " <<
                            id << ", batch_start_seq_id_ = " << batch_start_seq_id_ << ", batch_id_ = " << batch_id_ <<
                            ", seqs_.size() = " << seqs_.size() << " . Returning nullptr." << std::endl;
            return dummy_nullptr_seq_;
        }
        return seqs_[id];
    }

    const mindex::SequencePtr& GetSeqByAbsID(int64_t abs_id) const {
        int64_t id = abs_id - batch_start_seq_id_;
        return GetSeqByID(id);
    }

    void Add(std::unique_ptr<mindex::Sequence> seq) {
        if (seq == nullptr) {
            return;
        }
        int64_t new_id = static_cast<int64_t>(seqs_.size());
        seq->id(new_id);
        seq->abs_id(batch_start_seq_id_ + new_id);
        total_size_ += seq->len();
        seqs_.emplace_back(std::move(seq));
    }

    bool Serialize(std::ostream& ofs, mindex::SequenceFormat to_fmt);
    // bool Deserialize();

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

    std::string GetOpenFilePath() const {
        if (fp_handlers_ == nullptr) {
            return std::string();
        }
        return fp_handlers_->in_path();
    }

    mindex::SequenceFormat GetOpenFileFormat() const {
        if (curr_open_file_ < 0 || curr_open_file_ >= in_files_.size()) {
            return mindex::SequenceFormat::Unknown;
        }
        return curr_input_fmt_;
    }

    int64_t GetOpenFileTell() const {
        if (fp_handlers_ == nullptr) {
            return -2;
        }
        if (curr_open_file_ < 0 || curr_open_file_ >= in_files_.size()) {
            return -3;
        }
        int64_t ret_val = fp_handlers_->Tell();
        if (ret_val < 0) {
            return ret_val;
        }

        auto curr_fmt = GetOpenFileFormat();
        if (curr_fmt == mindex::SequenceFormat::Unknown) {
            return -4;
        }
        if (curr_fmt == mindex::SequenceFormat::Fasta || curr_fmt == mindex::SequenceFormat::Fastq) {
            if (fp_handlers_->fp_kseq == NULL) {
                return -5;
            }
            int64_t offset = static_cast<int64_t>(fp_handlers_->fp_kseq->f->end) - static_cast<int64_t>(fp_handlers_->fp_kseq->f->begin);
            if (offset < 0) {
                return -6;
            }
            // The "-1" because kseq parser loads the next '>' or '@' character too.
            int64_t offset_separator_char = (fp_handlers_->fp_kseq->last_char == '>' || fp_handlers_->fp_kseq->last_char == '@') ? (1) : 0;
            ret_val -= (offset + offset_separator_char);
        }

        // std::cerr << "gzip_tell = " << ret_val << ", fp_kseq = "
        //     << "[" << fp_handlers_->fp_kseq->f->begin
        //     << ", " << fp_handlers_->fp_kseq->f->end
        //     << ", " << fp_handlers_->fp_kseq->f->is_eof
        //     << "]\n";
        return ret_val;
    }

    int64_t GetOpenFileId() const {
        return curr_open_file_;
    }

private:
    SequenceFile();
    SequenceFile(const std::vector<std::string>& headers, const std::vector<std::string>& seqs);
    SequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt);
    SequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts);

    SequenceFile(const SequenceFile&) = delete;
    SequenceFile& operator=(const SequenceFile&) = delete;

    mindex::SequencePtr YieldSequence_(bool convert_to_uppercase) {
        if (in_files_.size() == 0) {
            WARNING_REPORT(ERR_OPENING_FILE, "No input files were specified, but YieldSequence_ was called. Returning nullptr.");
            return nullptr;
        }

        if (curr_open_file_ < 0 || fp_handlers_ == nullptr) {
            curr_open_file_ = 0;
            const std::string& next_in_path = std::get<0>(in_files_[curr_open_file_]);
            const mindex::SequenceFormat& next_in_fmt = std::get<1>(in_files_[curr_open_file_]);
            Open(next_in_path, next_in_fmt);
        }

        auto ret = mindex::SequenceDeserializer::DeserializeSequence(fp_handlers_, curr_input_fmt_, convert_to_uppercase);

        // If reading the input file did not yield a sequence,
        // it's either empty or something went wrong.
        // Open the next file in line.
        if (ret == nullptr) {
            ++curr_open_file_;
            if (curr_open_file_ < in_files_.size()) {
                const std::string& next_in_path = std::get<0>(in_files_[curr_open_file_]);
                const mindex::SequenceFormat& next_in_fmt = std::get<1>(in_files_[curr_open_file_]);
                Open(next_in_path, next_in_fmt);
            }
        }

        return ret;
    }

    std::vector<std::pair<std::string, mindex::SequenceFormat>> in_files_;
    int64_t curr_open_file_;
    mindex::SequenceFormat curr_input_fmt_;
    mindex::SequenceFileHandlersPtr fp_handlers_;
    std::vector<mindex::SequencePtr> seqs_;     // All sequences.
    int64_t batch_id_;              // The ID of the currently loaded batch.
    int64_t batch_start_seq_id_;    // Absolute ID of the first sequence in the loaded batch.
    int64_t total_size_;            // Sum of all sequence lengths.
    mindex::SequencePtr dummy_nullptr_seq_;
};

}

#endif
