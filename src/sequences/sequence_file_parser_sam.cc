/*
 * sequence_file_parser_sam.cc
 *
 *  Created on: Apr 28, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_sam.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParserSam(const std::string& in_path) {
    auto ret = mindex::SequenceFileParserBasePtr(new mindex::SequenceFileParserSam());
    bool rv = ret->Open(in_path);
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

SequenceFileParserSam::SequenceFileParserSam()
    : path_(), header_(), fp_gzip_(NULL), header_groups_() {
}

SequenceFileParserSam::~SequenceFileParserSam() {
    if (fp_gzip_ != NULL) {
        gzclose(fp_gzip_);
        fp_gzip_ = NULL;
    }
}

SequencePtr SequenceFileParserSam::YieldSequence() {
    int64_t id = -1;
    int64_t abs_id = -1;

    mindex::SequencePtr seq = nullptr;

    std::string line;

    while (raptor::ReadGZLine(fp_gzip_, line)) {
        if (line.size() == 0) {
            continue;
        }
        if (line[0] == '@') {
            continue;
        }

        std::istringstream ss(line);
        std::string qname, rname, cigar, rnext, seq_data, qual_data;
        int32_t flag = 0, mapq = 255;
        int64_t rpos = 0, pnext = 0, tlen = 0;

        ss >> qname >> flag >> rname >> rpos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq_data >> qual_data;

        if (seq_data.size() == 0 || seq_data == "*") {
            continue;
        }

        seq = mindex::createSequence();
        seq->data().insert(seq->data().end(), (int8_t *) &seq_data[0], (int8_t *) (&seq_data[0] + seq_data.size()));

        if (qual_data.size() > 0) {
            seq->qual().insert(seq->qual().end(), (int8_t *) &qual_data[0], (int8_t *) (&qual_data[0] + qual_data.size()));
        }

        seq->header(qname);
        seq->id(id);
        seq->abs_id(abs_id);

        // Parse additional tags of the alignment line.
        std::string tag_str;
        while (ss >> tag_str) {
            auto tokens = raptor::Tokenize(tag_str, ':');
            if (tokens.size() != 3) {
                WARNING_REPORT(ERR_UNEXPECTED_VALUE, "SAM tag has unexpected number of tokens. tokens.size() = %ld, expected = 3. Skipping.", static_cast<int32_t>(tokens.size()));
                continue;
            }
            raptor::SamTag tag(tokens[0], tokens[1], tokens[2]);
            seq->AddTag(tag);
        }

        // Always report the sequence in the FWD strand.
        if (flag & 0x10) {
            seq->ReverseComplement();
        }

        // We found a valid sequence. Yield.
        break;
    }

    return seq;
}

void SequenceFileParserSam::ParseHeader_() {
    std::string header;

    int64_t prev_pos = GetFileOffset();

    std::string line;
    while (raptor::ReadGZLine(fp_gzip_, line)) {
        if (line.size() == 0) {
            continue;
        }
        if (line[0] != '@') {
            break;
        }
        header += line + "\n";
        prev_pos = GetFileOffset();
    }
    FileSeek(prev_pos);
    header_ = header;
}

void SequenceFileParserSam::ParseReadGroupAndProgramGroupFromHeader_(const std::string& header) {
    header_groups_ = mindex::ParseReadGroupAndProgramGroupFromSAMHeader(header);
}

std::string SequenceFileParserSam::GetFileHeaderAsString() const {
    // FASTA/FASTQ doesn't have a header. The header_ is initialized
    // in the constructor.
    return header_;
}

int64_t SequenceFileParserSam::GetFileOffset() const {
    if (IsOpen() == false) {
        return -1;
    }
    return gztell(fp_gzip_);
}

const HeaderGroupType& SequenceFileParserSam::GetHeaderGroups() const {
    return header_groups_;
}

std::string SequenceFileParserSam::GetFilePath() const {
    return path_;
}

bool SequenceFileParserSam::IsOpen() const {
    if (path_.empty() || fp_gzip_ == nullptr) {
        return false;
    }
    return true;
}

bool SequenceFileParserSam::FileSeek(int64_t abs_pos) {
    if (IsOpen() == false) {
        return false;
    }
    gzseek(fp_gzip_, abs_pos, SEEK_SET);
    return true;
}

bool SequenceFileParserSam::Open(const std::string& path) {
    path_ = path;
    fp_gzip_ = gzopen(path_.c_str(), "r");

    if (fp_gzip_ == NULL) {
        FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", path_.c_str());
        return false;
    }

    ParseHeader_();
    ParseReadGroupAndProgramGroupFromHeader_(header_);

    return true;
}

}
