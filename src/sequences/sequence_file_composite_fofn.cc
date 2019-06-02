/*
 * sequence_file_composite_fofn.cc
 *
 *  Created on: Jun 02, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_composite_fofn.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>
#include <fstream>

namespace mindex {

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::string& fofn_path) {
    std::ifstream ifs(fofn_path);
    if (ifs.is_open() == false) {
        return nullptr;
    }
    std::string line;
    std::vector<std::string> in_paths;
    std::vector<mindex::SequenceFormat> in_fmts;
    while (std::getline(ifs, line)) {
        in_paths.emplace_back(line);
        in_fmts.emplace_back(mindex::SequenceFormat::Auto);
    }

    auto ret = mindex::SequenceFileCompositeBasePtr(new mindex::SequenceFileCompositeFofn(in_paths, in_fmts));

    if (ret == nullptr) {
        return nullptr;
    }

    return std::move(ret);
}

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt) {

    std::vector<mindex::SequenceFormat> in_fmts(in_paths.size(), in_fmt);

    auto ret = mindex::SequenceFileCompositeBasePtr(new mindex::SequenceFileCompositeFofn(in_paths, in_fmts));

    if (ret == nullptr) {
        return nullptr;
    }

    return std::move(ret);
}

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFofn(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts) {
    auto ret = mindex::SequenceFileCompositeBasePtr(new mindex::SequenceFileCompositeFofn(in_paths, in_fmts));

    if (ret == nullptr) {
        return nullptr;
    }

    return std::move(ret);
}


SequenceFileCompositeFofn::SequenceFileCompositeFofn()
    :
        files_(),
        formats_(),
        curr_open_file_(-1),
        curr_input_fmt_(mindex::SequenceFormat::Unknown),
        parser_(nullptr),
        convert_to_uppercase_(true),
        open_file_seq_offset_prev_(-1),
        open_file_seq_offset_curr_(-1),
        headers_(),
        header_groups_(),
        batch_id_(0),
        num_loaded_seqs_(0)
{

}

SequenceFileCompositeFofn::SequenceFileCompositeFofn(const std::vector<std::string>& in_paths, const std::vector<mindex::SequenceFormat>& in_fmts)
    :
        files_(in_paths),
        formats_(in_fmts),
        curr_open_file_(-1),
        curr_input_fmt_(mindex::SequenceFormat::Unknown),
        parser_(nullptr),
        convert_to_uppercase_(true),
        open_file_seq_offset_prev_(-1),
        open_file_seq_offset_curr_(-1),
        headers_(),
        header_groups_(),
        batch_id_(0),
        num_loaded_seqs_(0)
{
        ParseHeaders_();
}

SequenceFileCompositeFofn::~SequenceFileCompositeFofn() {

}

std::vector<std::string> SequenceFileCompositeFofn::GetFileHeaders() const {
    return headers_;
}

std::string SequenceFileCompositeFofn::GetCurrentFilePath() const {
    if (curr_open_file_ < 0 || parser_ == nullptr) {
        return {};
    }
    return files_[curr_open_file_];
}

int64_t SequenceFileCompositeFofn::GetCurrentFileId() const {
    return curr_open_file_;
}

int64_t SequenceFileCompositeFofn::GetCurrentFileOffset() const {
    return open_file_seq_offset_curr_;
}

int64_t SequenceFileCompositeFofn::GetCurrentFilePreviousOffset() const {
    return open_file_seq_offset_prev_;
}

bool SequenceFileCompositeFofn::FileSeek(int64_t file_id, int64_t abs_pos) {
    return false;
}

bool SequenceFileCompositeFofn::FileSeek(const std::string& file_name, int64_t abs_pos) {
    return false;
}

mindex::SequencePtr SequenceFileCompositeFofn::YieldSequence() {
    if (files_.size() == 0) {
        WARNING_REPORT(ERR_OPENING_FILE, "No input files were specified, but YieldSequence() was called. Returning nullptr.");
        return nullptr;
    }

    if (curr_open_file_ < 0 || parser_ == nullptr) {
        curr_open_file_ = 0;
        const std::string& next_in_path = files_[curr_open_file_];
        const mindex::SequenceFormat& next_in_fmt = formats_[curr_open_file_];

        Open_(next_in_path, next_in_fmt);
    }

    open_file_seq_offset_prev_ = GetOpenFileTell_();

    mindex::SequencePtr seq = nullptr;
    int64_t num_in_files = files_.size();

    while ((curr_open_file_ + 1) <= num_in_files && (seq = parser_->YieldSequence()) == nullptr) {
        ++curr_open_file_;
        if (curr_open_file_ >= files_.size()) {
            break;
        }
        const std::string& next_in_path = files_[curr_open_file_];
        const mindex::SequenceFormat& next_in_fmt = formats_[curr_open_file_];
        Open_(next_in_path, next_in_fmt);
        open_file_seq_offset_prev_ = GetOpenFileTell_();
    }

    open_file_seq_offset_curr_ = GetOpenFileTell_();

    if (seq != nullptr && convert_to_uppercase_) {
        seq->ToUppercase();
    }

    return seq;
}

mindex::SequenceFilePtr SequenceFileCompositeFofn::YieldBatchOfOne() {
    mindex::SequenceFilePtr seq_file = mindex::createSequenceFile();
    mindex::SequencePtr seq(nullptr);

    seq_file->batch_start_seq_id(num_loaded_seqs_);
    seq_file->batch_id(batch_id_);

    while ((seq = YieldSequence()) != nullptr) {
        seq_file->Add(std::move(seq));
        num_loaded_seqs_ += 1;
        break;
    }

    if (seq_file->seqs().size() == 0) {
        return nullptr;
    }

    ++batch_id_;

    return seq_file;
}

mindex::SequenceFilePtr SequenceFileCompositeFofn::YieldBatchMB(int64_t batch_size_mb) {
    mindex::SequenceFilePtr seq_file = mindex::createSequenceFile();

    int64_t total_size = 0;
    int64_t batch_size_b = batch_size_mb * 1024 * 1024;

    mindex::SequencePtr seq(nullptr);

    seq_file->batch_start_seq_id(num_loaded_seqs_);
    seq_file->batch_id(batch_id_);

    while ((seq = YieldSequence()) != nullptr) {
        total_size += seq->len();
        seq_file->Add(std::move(seq));
        num_loaded_seqs_ += 1;

        if (batch_size_b > 0 && total_size > batch_size_b) {
            break;
        }
    }

    if (seq_file->seqs().size() == 0) {
        return nullptr;
    }

    ++batch_id_;

    return seq_file;
}

mindex::SequenceFilePtr SequenceFileCompositeFofn::YieldAll() {
    mindex::SequenceFilePtr seq_file = mindex::createSequenceFile();

    mindex::SequencePtr seq(nullptr);

    seq_file->batch_start_seq_id(num_loaded_seqs_);
    seq_file->batch_id(batch_id_);

    while ((seq = YieldSequence()) != nullptr) {
        seq_file->Add(std::move(seq));
        num_loaded_seqs_ += 1;
    }

    if (seq_file->seqs().size() == 0) {
        return nullptr;
    }

    ++batch_id_;

    return seq_file;
}

const HeaderGroupType& SequenceFileCompositeFofn::GetHeaderGroups() const {
    return header_groups_;
}

mindex::SequenceFormat SequenceFileCompositeFofn::GetOpenFileFormat_() const {
    if (curr_open_file_ < 0 || curr_open_file_ >= files_.size()) {
        return mindex::SequenceFormat::Unknown;
    }
    return curr_input_fmt_;
}

bool SequenceFileCompositeFofn::Open_(const std::string& in_path, mindex::SequenceFormat in_fmt) {
    if (in_fmt == mindex::SequenceFormat::Auto) {
        in_fmt = GetSequenceFormatFromPath(in_path);
    }

    parser_ = mindex::createSequenceFileParser(in_path, in_fmt);
    curr_input_fmt_ = in_fmt;

    return (parser_ != nullptr);
}

int64_t SequenceFileCompositeFofn::GetOpenFileTell_() const {
    auto curr_fmt = GetOpenFileFormat_();
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

bool SequenceFileCompositeFofn::ParseHeaders_() {
    mindex::HeaderGroupType merged_groups;
    std::vector<std::string> headers;

    for (size_t fid = 0; fid < files_.size(); ++fid) {
        const std::string& in_path = files_[fid];
        mindex::SequenceFormat& in_fmt = formats_[fid];

        if (in_fmt == mindex::SequenceFormat::Auto) {
            in_fmt = GetSequenceFormatFromPath(in_path);
        }

        auto parser = mindex::createSequenceFileParser(in_path, in_fmt);
        headers.emplace_back(parser->GetFileHeaderAsString());
        const auto& curr_groups = parser->GetHeaderGroups();

        for (const auto& it_field: curr_groups) {
            for (const auto& it_ids: it_field.second) {
                merged_groups[it_field.first][it_ids.first] = it_ids.second;
            }
        }
    }

    std::swap(header_groups_, merged_groups);
    std::swap(headers, headers_);

    return true;
}

}
