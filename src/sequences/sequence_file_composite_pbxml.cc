/*
 * sequence_file_composite_pbxml.cc
 *
 *  Created on: Jun 02, 2019
 *      Author: Ivan Sovic
 */

#ifdef RAPTOR_COMPILED_WITH_PBBAM

#include <sequences/sequence_file_composite_pbxml.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>
#include <fstream>

#include <pbbam/PbiRawData.h>
#include <pbbam/PbiIndexedBamReader.h>

namespace mindex {

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositePbXml(const std::string& in_path) {
    auto ret = mindex::SequenceFileCompositeBasePtr(new mindex::SequenceFileCompositePbXml(in_path));

    if (ret == nullptr) {
        return nullptr;
    }

    return std::move(ret);
}

SequenceFileCompositePbXml::SequenceFileCompositePbXml()
    :
        in_xml_path_(),
        files_(),
        file_offsets_(),
        open_file_id_(-1),
        open_file_offset_id_(0),
        dataset_(nullptr),
        parser_(nullptr),
        convert_to_uppercase_(true),
        open_file_seq_offset_prev_(-1),
        open_file_seq_offset_curr_(-1),
        headers_(),
        header_groups_(),
        merged_header_(),
        batch_id_(0),
        num_loaded_seqs_(0)
{

}

SequenceFileCompositePbXml::SequenceFileCompositePbXml(const std::string& in_path)
    :
        in_xml_path_(in_path),
        files_(),
        file_offsets_(),
        open_file_id_(-1),
        open_file_offset_id_(0),
        dataset_(nullptr),
        parser_(nullptr),
        convert_to_uppercase_(true),
        open_file_seq_offset_prev_(-1),
        open_file_seq_offset_curr_(-1),
        headers_(),
        header_groups_(),
        merged_header_(),
        batch_id_(0),
        num_loaded_seqs_(0)
{
        Open_(in_path);
}

SequenceFileCompositePbXml::~SequenceFileCompositePbXml() {

}

std::vector<std::string> SequenceFileCompositePbXml::GetFiles() const {
    return files_;
}

std::vector<std::string> SequenceFileCompositePbXml::GetFileHeaders() const {
    return headers_;
}

std::string SequenceFileCompositePbXml::GetCurrentFilePath() const {
    if (open_file_id_ < 0 || parser_ == nullptr) {
        return {};
    }
    return files_[open_file_id_];
}

int64_t SequenceFileCompositePbXml::GetCurrentFileId() const {
    return open_file_id_;
}

int64_t SequenceFileCompositePbXml::GetCurrentFileOffset() const {
    return open_file_seq_offset_curr_;
}

int64_t SequenceFileCompositePbXml::GetCurrentFilePreviousOffset() const {
    return open_file_seq_offset_prev_;
}

bool SequenceFileCompositePbXml::FileSeek(int64_t file_id, int64_t abs_pos) {
    return false;
}

bool SequenceFileCompositePbXml::FileSeek(const std::string& file_name, int64_t abs_pos) {
    return false;
}

mindex::SequencePtr SequenceFileCompositePbXml::YieldSequence() {
    mindex::SequencePtr seq = nullptr;

    if (files_.size() == 0) {
        return seq;
    }
    if (open_file_id_ >= static_cast<int64_t>(files_.size())) {
        return seq;
    }

    bool do_open = false;
    if (open_file_id_ < 0) {
        open_file_id_ = 0;
        open_file_offset_id_ = 0;
        do_open = true;
    }
    if (open_file_offset_id_ >= static_cast<int64_t>(file_offsets_[open_file_id_].size())) {
        ++open_file_id_;
        open_file_offset_id_ = 0;
        if (open_file_id_ >= static_cast<int64_t>(files_.size())) {
            return seq;
        }
        do_open = true;
    }

    if (do_open) {
        parser_ = mindex::createSequenceFileParser(files_[open_file_id_], mindex::SequenceFormat::BAM);
    }

    parser_->FileSeek(file_offsets_[open_file_id_][open_file_offset_id_]);
    seq = parser_->YieldSequence();

    open_file_seq_offset_prev_ = file_offsets_[open_file_id_][open_file_offset_id_];
    open_file_seq_offset_curr_ = parser_->GetFileOffset();

    ++open_file_offset_id_;

    return seq;
}

mindex::SequenceFilePtr SequenceFileCompositePbXml::YieldBatchOfOne() {
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

mindex::SequenceFilePtr SequenceFileCompositePbXml::YieldBatchMB(int64_t batch_size_mb) {
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

mindex::SequenceFilePtr SequenceFileCompositePbXml::YieldAll() {
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

const HeaderGroupType& SequenceFileCompositePbXml::GetHeaderGroups() const {
    return header_groups_;
}

bool SequenceFileCompositePbXml::Open_(const std::string& in_path) {
    files_.clear();
    file_offsets_.clear();
    open_file_id_ = -1;
    open_file_offset_id_ = -1;
    parser_ = nullptr;

    dataset_ = std::make_unique<PacBio::BAM::DataSet>(in_path);

    const PacBio::BAM::PbiIndexCache pbiCache = PacBio::BAM::MakePbiIndexCache(*dataset_);
    const PacBio::BAM::PbiFilter filter = PacBio::BAM::PbiFilter::FromDataSet(*dataset_);

    size_t file_number = 0;
    for (const PacBio::BAM::BamFile& bam : dataset_->BamFiles()) {
        const auto& fn = bam.Filename();

        const std::shared_ptr<PacBio::BAM::PbiRawData>& file_index = pbiCache->at(file_number);

        PacBio::BAM::PbiIndexedBamReader reader{filter, bam, file_index};
        const PacBio::BAM::IndexResultBlocks& blocks = reader.IndexBlocks();

        std::vector<int64_t> offsets;
        const auto& pbi_offsets = file_index->BasicData().fileOffset_;
        int64_t num_added = 0;
        for (const auto& block : blocks) {
            // if (block.numReads_ <= 0) {
                // ERROR_REPORT(ERR_FILE_READ_DATA, "block.numReads_ = %ld", block.numReads_);
            // }

            for (size_t i = 0; i < block.numReads_; ++i) {
                offsets.emplace_back(pbi_offsets.at(block.firstIndex_ + i));
            }

            num_added += block.numReads_;
        }

        if (num_added) {
            files_.emplace_back(fn);
            file_offsets_.emplace_back(std::move(offsets));
        }

        ++file_number;
    }

    if (files_.size() == 0) {
        return false;
    }

    open_file_id_ = 0;
    open_file_offset_id_ = 0;
    parser_ = mindex::createSequenceFileParser(files_[0], mindex::SequenceFormat::BAM);

    ParseHeaders_();

    return true;
}

int64_t SequenceFileCompositePbXml::GetOpenFileTell_() const {
    if (parser_ == nullptr) {
        return -2;
    }
    if (parser_->IsOpen() == false) {
        return -3;
    }
    return parser_->GetFileOffset();
}

bool SequenceFileCompositePbXml::ParseHeaders_() {
    PacBio::BAM::BamHeader merged_header;
    std::vector<std::string> headers;
    mindex::HeaderGroupType merged_groups;
    bool header_initialized = false;

    // Merge all the read groups and additional PacBio info.
    const auto& bam_files = dataset_->BamFiles();
    for (auto& bam_file : bam_files) {
        auto header = bam_file.Header();

        // Collect string representations of the headers.
        std::string header_sam = header.ToSam();
        headers.emplace_back(header_sam);

        // Get all the groups from the header.
        auto curr_groups = mindex::ParseReadGroupAndProgramGroupFromSAMHeader(header_sam);
        for (const auto& it_field: curr_groups) {
            for (const auto& it_ids: it_field.second) {
                merged_groups[it_field.first][it_ids.first] = it_ids.second;
            }
        }

        // Use Pbbam's merge capability to merge the headers.
        if (!header_initialized) {
            merged_header = header.DeepCopy();
            header_initialized = true;
        } else {
            merged_header += header;
        }
    }

    // Update the members.
    merged_header_ = merged_header.ToSam();
    std::swap(header_groups_, merged_groups);
    std::swap(headers, headers_);

    return true;
}

}

#endif
