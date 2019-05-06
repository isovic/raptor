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

#include <boost/lexical_cast.hpp>

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
    ParseReadGroupAndProgramGroupFromHeader_(sam_header_);

    return true;
}

void SequenceFileParserBam::ParseReadGroupAndProgramGroupFromHeader_(const std::string& header) {
    std::istringstream iss(header);
    std::string line;
    while (std::getline(iss, line)) {
        if (line.size() < 3) {
            continue;
        }
        if (line[0] != '@') {
            break;
        }

        std::string field_name = line.substr(1, 2);

        if (field_name == "RG" || field_name == "PG") {
            std::vector<HeaderTag> tags;

            // Skip the first token, that's "@RG" or "@PG".
            auto tokens = raptor::Tokenize(line, '\t');
            bool is_ok = false;
            std::string field_id;
            for (size_t tid = 1; tid < tokens.size(); ++tid) {
                std::string tag_name = tokens[tid].substr(0, 2);
                std::string tag_val = tokens[tid].substr(3);
                tags.emplace_back(HeaderTag{tag_name, tag_val});
                if (tag_name == "ID") {
                    field_id = tag_name;
                    is_ok = true;
                }
            }

            if (is_ok == false) {
                WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Skipping faulty header line: '%s'.", line.c_str());
                continue;
            }

            header_groups_[field_name][field_id] = std::move(tags);
        }
    }
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

    std::string joined_tags = EncodeSamTags_(record.Impl().Tags());
    std::istringstream iss(joined_tags);
    std::string tag_str;
    while ((iss >> tag_str)) {
        auto tokens = raptor::Tokenize(tag_str, ':');
        if (tokens.size() != 3) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "SAM tag has unexpected number of tokens. tokens.size() = %ld, expected = 3. Skipping.", static_cast<int32_t>(tokens.size()));
            continue;
        }
        raptor::SamTag tag(tokens[0], tokens[1], tokens[2]);
        seq->AddTag(tag);
    }

    return seq;
}

const HeaderGroupType& SequenceFileParserBam::GetHeaderGroups() const {
    return header_groups_;
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




template <typename T>
inline void appendSamValue(const T& value, std::string& result)
{
    result.append(boost::lexical_cast<std::string>(value));
}

template <typename T>
inline void appendSamValue_8bit(const T& value, std::string& result)
{
    result.append(boost::lexical_cast<std::string>(static_cast<int>(value)));
}

template <typename T>
void appendSamMultiValue(const T& container, std::string& result)
{
    for (const auto x : container) {
        result.push_back(',');
        appendSamValue(x, result);
    }
}

template <typename T>
void appendSamMultiValue_8bit(const T& container, std::string& result)
{
    for (const auto x : container) {
        result.push_back(',');
        appendSamValue_8bit(x, result);
    }
}

std::string SequenceFileParserBam::EncodeSamTags_(const PacBio::BAM::TagCollection& tags) const {
    std::string result;
    result.reserve(1024);

    for (const auto& tagIter : tags) {
        const auto& name = tagIter.first;
        if (name.size() != 2)
            throw std::runtime_error{"SamTagCodec: malformatted tag name: " + name};

        const auto& tag = tagIter.second;
        if (tag.IsNull()) continue;

        // tab separator
        if (!result.empty()) result.push_back('\t');

        // "<TAG>:"
        result.append(name);
        result.push_back(':');

        // "<TYPE>:<DATA>" for printable, ASCII char
        if (tag.HasModifier(PacBio::BAM::TagModifier::ASCII_CHAR)) {
            const auto c = tag.ToAscii();
            if (c != '\0') {
                result.push_back('A');
                result.push_back(':');
                result.push_back(c);
                continue;
            }
        }

        // "<TYPE>:<DATA>" for all other data
        switch (tag.Type()) {
            case PacBio::BAM::TagDataType::INT8:
                result.push_back('i');
                result.push_back(':');
                appendSamValue_8bit(tag.ToInt8(), result);
                break;
            case PacBio::BAM::TagDataType::UINT8:
                result.push_back('i');
                result.push_back(':');
                appendSamValue_8bit(tag.ToUInt8(), result);
                break;
            case PacBio::BAM::TagDataType::INT16:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToInt16(), result);
                break;
            case PacBio::BAM::TagDataType::UINT16:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToUInt16(), result);
                break;
            case PacBio::BAM::TagDataType::INT32:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToInt32(), result);
                break;
            case PacBio::BAM::TagDataType::UINT32:
                result.push_back('i');
                result.push_back(':');
                appendSamValue(tag.ToUInt32(), result);
                break;
            case PacBio::BAM::TagDataType::FLOAT:
                result.push_back('f');
                result.push_back(':');
                appendSamValue(tag.ToFloat(), result);
                break;

            case PacBio::BAM::TagDataType::STRING: {
                result.push_back(tag.HasModifier(PacBio::BAM::TagModifier::HEX_STRING) ? 'H' : 'Z');
                result.push_back(':');
                result.append(tag.ToString());
                break;
            }

            case PacBio::BAM::TagDataType::INT8_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('c');
                appendSamMultiValue_8bit(tag.ToInt8Array(), result);
                break;
            case PacBio::BAM::TagDataType::UINT8_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('C');
                appendSamMultiValue_8bit(tag.ToUInt8Array(), result);
                break;
            case PacBio::BAM::TagDataType::INT16_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('s');
                appendSamMultiValue(tag.ToInt16Array(), result);
                break;
            case PacBio::BAM::TagDataType::UINT16_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('S');
                appendSamMultiValue(tag.ToUInt16Array(), result);
                break;
            case PacBio::BAM::TagDataType::INT32_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('i');
                appendSamMultiValue(tag.ToInt32Array(), result);
                break;
            case PacBio::BAM::TagDataType::UINT32_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('I');
                appendSamMultiValue(tag.ToUInt32Array(), result);
                break;
            case PacBio::BAM::TagDataType::FLOAT_ARRAY:
                result.push_back('B');
                result.push_back(':');
                result.push_back('f');
                appendSamMultiValue(tag.ToFloatArray(), result);
                break;

            default:
                throw std::runtime_error{"SamTagCodec: unsupported tag-type encountered: " +
                                         std::to_string(static_cast<uint16_t>(tag.Type()))};
        }
    }

    return result;
}

}
