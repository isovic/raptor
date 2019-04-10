/*
 * sequence_file.cc
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#include <index/sequence_file.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>

namespace mindex {

mindex::SequenceFilePtr createSequenceFile() {
    return mindex::SequenceFilePtr(new mindex::SequenceFile());
}

mindex::SequenceFilePtr createSequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt) {
    return mindex::SequenceFilePtr(new mindex::SequenceFile(in_path, in_fmt));
}

mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts) {
    return mindex::SequenceFilePtr(new mindex::SequenceFile(in_paths, in_fmts));
}

// mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& headers, const std::vector<std::string>& seqs) {
//     return mindex::SequenceFilePtr(new mindex::SequenceFile(headers, seqs));
// }

bool SequenceFile::Open(const std::string& in_path, mindex::SequenceFormat in_fmt) {
    // All file handlers and mem allocations should be released automatically upon destruction.
    fp_handlers_ = mindex::createSequenceFileHandlers(in_path);
    curr_input_fmt_ = in_fmt;

    if (in_fmt == mindex::SequenceFormat::Auto) {
        auto ext = raptor::GetFileExt(in_path);

        // If the file is Gzipped, the .gz will be in the ext.
        // E.g. the output from GetFileExt can be "fasta.gz".
        if (ext.size() >= 3 && ext.substr(ext.size() - 3) == ".gz") {
            ext = ext.substr(0, ext.size() - 3);
        }

        curr_input_fmt_ =
                    (ext == "fasta") ? mindex::SequenceFormat::Fasta :
                    (ext == "fa") ? mindex::SequenceFormat::Fasta :
                    (ext == "fastq") ? mindex::SequenceFormat::Fastq :
                    (ext == "fq") ? mindex::SequenceFormat::Fastq :
                    (ext == "gfa") ? mindex::SequenceFormat::GFA :
                    (ext == "gfa1") ? mindex::SequenceFormat::GFA1 :
                    (ext == "gfa2") ? mindex::SequenceFormat::GFA2 :
                    (ext == "sam") ? mindex::SequenceFormat::SAM :
                    mindex::SequenceFormat::Unknown;
    }

    return (fp_handlers_ != nullptr);
}

SequenceFile::SequenceFile()
                :
                    in_files_(),
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    fp_handlers_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr) {

}

SequenceFile::SequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt)
                :
                    in_files_{std::make_pair(in_path, in_fmt)},
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    fp_handlers_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr) {
}

SequenceFile::SequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts)
                :
                    in_files_{},
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    fp_handlers_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr) {
    if (in_paths.size() != in_fmts.size()) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "in_paths.size() = %ld and in_fmts.size() = %ld. Not setting the inputs.", in_paths.size(), in_fmts.size());
        return;
    }
    for (size_t i = 0; i < in_paths.size(); ++i) {
        in_files_.emplace_back(std::make_pair(in_paths[i], in_fmts[i]));
    }
}

bool SequenceFile::LoadAll(const std::vector<std::string>& headers, const std::vector<std::string>& seqs, bool convert_to_uppercase) {
    batch_start_seq_id_ += seqs_.size();
    seqs_.clear();

    if (headers.size() != seqs.size()) {
        return false;
    }

    int64_t total_size = 0;
    int64_t id = 0;
    int64_t abs_id = batch_start_seq_id_;

    for (size_t i = 0; i < seqs.size(); ++i) {
        auto seq = createSequence(headers[i], seqs[i]);
        total_size += seq->len();
        this->Add(std::move(seq));
    }

    total_size_ = total_size;

    return true;
}

bool SequenceFile::LoadAll(bool convert_to_uppercase) {
    return LoadBatchMB(-1, convert_to_uppercase);
}

bool SequenceFile::LoadBatchMB(int64_t batch_size_mb, bool convert_to_uppercase) {
    batch_start_seq_id_ += seqs_.size();
    seqs_.clear();

    mindex::SequencePtr seq(nullptr);
    int64_t total_size = 0;
    int64_t batch_size_b = batch_size_mb * 1024 * 1024;

    total_size_ = 0;

    while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr || curr_open_file_ < in_files_.size()) {
        if (seq == nullptr) {
            continue;
        }

        total_size += seq->len();
        this->Add(std::move(seq));

        if (batch_size_b > 0 && total_size > batch_size_b) {
            break;
        }
    }

    total_size_ = total_size;

    ++batch_id_;

    if (seqs_.size() == 0) {
        return false;
    }

    return true;
}

bool SequenceFile::LoadBatchOfOne(bool convert_to_uppercase) {
    total_size_ = 0;
    batch_start_seq_id_ += seqs_.size();
    seqs_.clear();

    mindex::SequencePtr seq(nullptr);

    // Find the first non-nullptr sequence.
    while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr || curr_open_file_ < in_files_.size()) {
        if (seq != nullptr) {
            break;
        }
    }

    if (seq == nullptr) {
        return false;
    }

    this->Add(std::move(seq));

    ++batch_id_;

    return true;
}

}
