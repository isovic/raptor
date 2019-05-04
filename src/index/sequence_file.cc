/*
 * sequence_file.cc
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#include <index/sequence_file.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>

namespace mindex {

mindex::SequenceFilePtr createSequenceFile() {
    return mindex::SequenceFilePtr(new mindex::SequenceFile());
}

mindex::SequenceFilePtr createSequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt) {
    auto ret = mindex::SequenceFilePtr(new mindex::SequenceFile(in_path, in_fmt));
    ret->ParseAndMergeHeaderGroups_();
    return ret;
}

mindex::SequenceFilePtr createSequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts) {
    auto ret = mindex::SequenceFilePtr(new mindex::SequenceFile(in_paths, in_fmts));
    ret->ParseAndMergeHeaderGroups_();
    return ret;
}

SequenceFile::SequenceFile()
                :
                    in_files_(),
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    parser_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr),
                    header_groups_(),
                    open_file_offset_prev_(0) {

}

SequenceFile::SequenceFile(const std::string& in_path, mindex::SequenceFormat in_fmt)
                :
                    in_files_{std::make_pair(in_path, in_fmt)},
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    parser_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr),
                    header_groups_(),
                    open_file_offset_prev_(0) {
}

SequenceFile::SequenceFile(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts)
                :
                    in_files_{},
                    curr_open_file_(-1),
                    curr_input_fmt_(mindex::SequenceFormat::Unknown),
                    parser_(nullptr),
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr),
                    header_groups_(),
                    open_file_offset_prev_(0) {
    if (in_paths.size() != in_fmts.size()) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "in_paths.size() = %ld and in_fmts.size() = %ld. Not setting the inputs.", in_paths.size(), in_fmts.size());
        return;
    }
    for (size_t i = 0; i < in_paths.size(); ++i) {
        in_files_.emplace_back(std::make_pair(in_paths[i], in_fmts[i]));
    }
}

const mindex::SequencePtr& SequenceFile::GetSeqByID(int64_t id) const {
    if (id < 0 || id >= seqs_.size()) {
        std::cerr << "[GetSeqByID] Warnning: Requested id is out of scope of current batch. id = " <<
                        id << ", batch_start_seq_id_ = " << batch_start_seq_id_ << ", batch_id_ = " << batch_id_ <<
                        ", seqs_.size() = " << seqs_.size() << " . Returning nullptr." << std::endl;
        return dummy_nullptr_seq_;
    }
    return seqs_[id];
}

const mindex::SequencePtr& SequenceFile::GetSeqByAbsID(int64_t abs_id) const {
    int64_t id = abs_id - batch_start_seq_id_;
    return GetSeqByID(id);
}

void SequenceFile::Add(std::unique_ptr<mindex::Sequence> seq) {
    if (seq == nullptr) {
        return;
    }
    int64_t new_id = static_cast<int64_t>(seqs_.size());
    seq->id(new_id);
    seq->abs_id(batch_start_seq_id_ + new_id);
    total_size_ += seq->len();
    seqs_.emplace_back(std::move(seq));
}

std::string SequenceFile::GetOpenFilePath() const {
    if (parser_ == nullptr) {
        return std::string();
    }
    return parser_->GetFilePath();
}

mindex::SequenceFormat SequenceFile::GetOpenFileFormat() const {
    if (curr_open_file_ < 0 || curr_open_file_ >= in_files_.size()) {
        return mindex::SequenceFormat::Unknown;
    }
    return curr_input_fmt_;
}

int64_t SequenceFile::GetOpenFileTell() const {
    auto curr_fmt = GetOpenFileFormat();
    if (curr_fmt == mindex::SequenceFormat::Unknown) {
        return -4;
    }
    if (parser_ == nullptr) {
        return -2;
    }
    if (parser_->IsOpen() == false) {
        return -3;
    }
    return parser_->GetFileOffset();
}

int64_t SequenceFile::GetOpenFilePrevTell() const {
    auto curr_fmt = GetOpenFileFormat();
    if (curr_fmt == mindex::SequenceFormat::Unknown) {
        return -4;
    }
    if (parser_ == nullptr) {
        return -2;
    }
    if (parser_->IsOpen() == false) {
        return -3;
    }
    return open_file_offset_prev_;
}

mindex::SequencePtr SequenceFile::YieldSequence_(bool convert_to_uppercase) {
    if (in_files_.size() == 0) {
        WARNING_REPORT(ERR_OPENING_FILE, "No input files were specified, but YieldSequence_ was called. Returning nullptr.");
        return nullptr;
    }

    if (curr_open_file_ < 0 || parser_ == nullptr) {
        curr_open_file_ = 0;
        const std::string& next_in_path = std::get<0>(in_files_[curr_open_file_]);
        const mindex::SequenceFormat& next_in_fmt = std::get<1>(in_files_[curr_open_file_]);

        Open_(next_in_path, next_in_fmt);
    }

    open_file_offset_prev_ = GetOpenFileTell();

    mindex::SequencePtr seq = nullptr;
    int64_t num_in_files = in_files_.size();

    while ((curr_open_file_ + 1) <= num_in_files && (seq = parser_->YieldSequence()) == nullptr) {
        ++curr_open_file_;
        if (curr_open_file_ >= in_files_.size()) {
            break;
        }
        const std::string& next_in_path = std::get<0>(in_files_[curr_open_file_]);
        const mindex::SequenceFormat& next_in_fmt = std::get<1>(in_files_[curr_open_file_]);
        Open_(next_in_path, next_in_fmt);
        open_file_offset_prev_ = GetOpenFileTell();
    }

    if (seq != nullptr && convert_to_uppercase) {
        seq->ToUppercase();
    }

    return seq;
}

bool SequenceFile::ParseAndMergeHeaderGroups_() {
    mindex::HeaderGroupType merged_groups;

    for (size_t fid = 0; fid < in_files_.size(); ++fid) {
        const std::string& in_path = std::get<0>(in_files_[fid]);
        mindex::SequenceFormat& in_fmt = std::get<1>(in_files_[fid]);

        if (in_fmt == mindex::SequenceFormat::Auto) {
            auto ext = raptor::GetFileExt(in_path);
            // If the file is Gzipped, the .gz will be in the ext.
            // E.g. the output from GetFileExt can be "fasta.gz".
            if (ext.size() >= 3 && ext.substr(ext.size() - 3) == ".gz") {
                ext = ext.substr(0, ext.size() - 3);
            }
            in_fmt = SequenceFormatFromString(ext);
        }

        auto parser = mindex::createSequenceFileParser(in_path, in_fmt);
        const auto& curr_groups = parser->GetHeaderGroups();

        for (const auto& it_field: curr_groups) {
            for (const auto& it_ids: it_field.second) {
                merged_groups[it_field.first][it_ids.first] = it_ids.second;
            }
        }
    }

    std::swap(header_groups_, merged_groups);

    return true;
}

bool SequenceFile::Open_(const std::string& in_path, mindex::SequenceFormat in_fmt) {
    if (in_fmt == mindex::SequenceFormat::Auto) {
        auto ext = raptor::GetFileExt(in_path);
        // If the file is Gzipped, the .gz will be in the ext.
        // E.g. the output from GetFileExt can be "fasta.gz".
        if (ext.size() >= 3 && ext.substr(ext.size() - 3) == ".gz") {
            ext = ext.substr(0, ext.size() - 3);
        }
        in_fmt = SequenceFormatFromString(ext);
    }

    parser_ = mindex::createSequenceFileParser(in_path, in_fmt);
    curr_input_fmt_ = in_fmt;

    return (parser_ != nullptr);
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

    // while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr || curr_open_file_ < in_files_.size()) {
    while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr) {
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
    // while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr || curr_open_file_ < in_files_.size()) {
    while ((seq = YieldSequence_(convert_to_uppercase)) != nullptr) {
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
