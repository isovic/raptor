/*
 * sequence_file_parser_bam.cc
 *
 *  Created on: May 03, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_bam.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/BamRecord.h>
#include <pbbam/SamTagCodec.h>
#include <utility/stringutil.h>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParserBam(const std::string& in_path) {
    auto ret = mindex::SequenceFileParserBasePtr(new mindex::SequenceFileParserBam());
    bool rv = ret->Open(in_path);
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

SequenceFileParserBam::SequenceFileParserBam()
    : path_(), bam_reader_(nullptr), bam_file_(nullptr), dataset_(nullptr), sam_header_() {
}

SequenceFileParserBam::~SequenceFileParserBam() {

}

bool SequenceFileParserBam::Open(const std::string& path) {
    path_ = path;

    bam_file_ = std::make_unique<PacBio::BAM::BamFile>(path);
    bam_reader_ = std::make_unique<PacBio::BAM::BamReader>(*bam_file_);
    dataset_ = std::make_unique<PacBio::BAM::DataSet>(*bam_file_);

    ParseHeader_();

    return true;
}

void SequenceFileParserBam::ParseHeader_() {
    auto input_header = bam_reader_->Header();
    sam_header_ = input_header.ToSam();
}

SequencePtr SequenceFileParserBam::YieldSequence() {
    mindex::SequencePtr seq = nullptr;

    int64_t id = -1;
    int64_t abs_id = -1;

    PacBio::BAM::BamRecord record;
    bool rv_get = bam_reader_->GetNext(record);

    if (rv_get == false) {
        return nullptr;
    }

    // The "PacBio::BAM::Orientation::NATIVE" will always return the FWD strand of the sequence.
    std::string seq_str = record.Sequence(PacBio::BAM::Orientation::NATIVE);
    PacBio::BAM::QualityValues qual_data = record.Qualities(PacBio::BAM::Orientation::NATIVE);
    std::string qual_str = qual_data.Fastq();

    seq = mindex::createSequence();
    seq->data().insert(seq->data().end(), (int8_t *) &seq_str[0], (int8_t *) (&seq_str[0] + seq_str.size()));

    if (qual_data.size() > 0) {
        seq->qual().insert(seq->qual().end(), (int8_t *) &qual_str[0], (int8_t *) (&qual_str[0] + qual_str.size()));
    }

    seq->header(record.FullName());
    seq->id(id);
    seq->abs_id(abs_id);

    std::string joined_tags = PacBio::BAM::SamTagCodec::Encode(record.Impl().Tags());
    std::istringstream iss(joined_tags);
    std::string tag_str;
    while ((iss >> tag_str)) {
        auto tokens = raptor::Tokenize(tag_str, ':');
        raptor::SamTag tag(tokens[0], tokens[1], tokens[2]);
        seq->AddTag(tag);
    }

    return seq;
}

std::string SequenceFileParserBam::GetFileHeaderAsString() const {
    return sam_header_;
}

int64_t SequenceFileParserBam::GetFileOffset() const {
    return bam_reader_->VirtualTell();
}

std::string SequenceFileParserBam::GetFilePath() const {
    return path_;
}

bool SequenceFileParserBam::IsOpen() const {
    return bam_reader_ != nullptr;
}

bool SequenceFileParserBam::FileSeek(int64_t abs_pos) {
    bam_reader_->VirtualSeek(abs_pos);
    return true;
}

}
