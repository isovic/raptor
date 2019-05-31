/*
 * sequence_file_parser_bam.cc
 *
 *  Created on: May 03, 2019
 *      Author: Ivan Sovic
 */

#ifdef RAPTOR_COMPILED_WITH_PBBAM

#include <sequences/sequence_file_parser_pbxml.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/PbiRawData.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <pbbam/PbiBasicTypes.h>

#include <sequences/sequence_file_parser_bam.h>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParserPbXml(const std::string& in_path) {
    auto ret = mindex::SequenceFileParserBasePtr(new mindex::SequenceFileParserPbXml());
    bool rv = ret->Open(in_path);
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

SequenceFileParserPbXml::SequenceFileParserPbXml()
    : path_(), bam_file_(nullptr), dataset_(nullptr),
        sam_header_(), header_groups_(), files_(), file_offsets_(),
        open_file_id_(-1), open_file_offset_id_(-1), bam_parser_(nullptr) {
}

SequenceFileParserPbXml::~SequenceFileParserPbXml() {

}

bool SequenceFileParserPbXml::Open(const std::string& path) {
    path_ = path;

    files_.clear();
    file_offsets_.clear();
    open_file_id_ = -1;
    open_file_offset_id_ = -1;
    bam_parser_ = nullptr;

    dataset_ = std::make_unique<PacBio::BAM::DataSet>(path);

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
        for (const auto& block : blocks) {
            if (block.numReads_ <= 0) {
                ERROR_REPORT(ERR_FILE_READ_DATA, "block.numReads_ = %ld", block.numReads_);
            }

            for (size_t i = 0; i < block.numReads_; ++i) {
                offsets.emplace_back(pbi_offsets.at(block.firstIndex_ + i));
            }
        }

        files_.emplace_back(fn);
        file_offsets_.emplace_back(std::move(offsets));

        ++file_number;
    }

    if (files_.size() == 0) {
        return false;
    }

    open_file_id_ = 0;
    open_file_offset_id_ = 0;
    bam_parser_ = mindex::createSequenceFileParserBam(files_[open_file_id_]);

    ParseHeader_();

    ParseReadGroupAndProgramGroupFromHeader_(sam_header_);

    return true;
}

void SequenceFileParserPbXml::ParseReadGroupAndProgramGroupFromHeader_(const std::string& header) {
    header_groups_ = mindex::ParseReadGroupAndProgramGroupFromSAMHeader(header);
}

void SequenceFileParserPbXml::ParseHeader_() {
    PacBio::BAM::BamHeader merged_header;
    bool header_initialized = false;

    // Merge all the read groups and additional PacBio info.
    const auto& bam_files = dataset_->BamFiles();
    for (auto& bam_file : bam_files) {
        auto header = bam_file.Header();
        if (!header_initialized) {
            merged_header = header.DeepCopy();
            header_initialized = true;
        } else {
            merged_header += header;
        }
    }

    sam_header_ = merged_header.ToSam();
}

SequencePtr SequenceFileParserPbXml::YieldSequence() {
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
        bam_parser_ = mindex::createSequenceFileParserBam(files_[open_file_id_]);
    }

    // std::cerr << "[open_file_id_ = " << open_file_id_ << ", open_file_offset_id_ = " << open_file_offset_id_ << "] File: '" << files_[open_file_id_] << ", offset: " << file_offsets_[open_file_id_][open_file_offset_id_] << "\n";

    bam_parser_->FileSeek(file_offsets_[open_file_id_][open_file_offset_id_]);
    seq = bam_parser_->YieldSequence();
    ++open_file_offset_id_;

    return seq;
}

const HeaderGroupType& SequenceFileParserPbXml::GetHeaderGroups() const {
    return header_groups_;
}

std::string SequenceFileParserPbXml::GetFileHeaderAsString() const {
    return sam_header_;
}

int64_t SequenceFileParserPbXml::GetFileOffset() const {
    if (bam_parser_) {
        return bam_parser_->GetFileOffset();
    }
    return -1;
}

std::string SequenceFileParserPbXml::GetFilePath() const {
    return path_;
}

bool SequenceFileParserPbXml::IsOpen() const {
    if (bam_parser_ == nullptr) {
        return false;
    }
    if (open_file_id_ < 0) {
        return false;
    }
    return true;
}

bool SequenceFileParserPbXml::FileSeek(int64_t abs_pos) {
    if (bam_parser_) {
        return bam_parser_->FileSeek(abs_pos);
    }
    return true;
}

}

#endif
