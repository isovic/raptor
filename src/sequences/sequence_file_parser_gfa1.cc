/*
 * sequence_file_parser_gfa2.cc
 *
 *  Created on: Apr 28, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_gfa1.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>
#include <iostream>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParserGfa1(const std::string& in_path) {
    auto ret = mindex::SequenceFileParserBasePtr(new mindex::SequenceFileParserGfa1());
    bool rv = ret->Open(in_path);
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

SequenceFileParserGfa1::SequenceFileParserGfa1()
    : path_(), fp_gzip_(NULL) {
}

SequenceFileParserGfa1::~SequenceFileParserGfa1() {
    if (fp_gzip_ != NULL) {
        gzclose(fp_gzip_);
        fp_gzip_ = NULL;
    }
}

SequencePtr SequenceFileParserGfa1::YieldSequence() {
    int64_t id = -1;
    int64_t abs_id = -1;

    mindex::SequencePtr seq = nullptr;

    std::string line;
    while (raptor::ReadGZLine(fp_gzip_, line)) {
        if (line.size() == 0) {
            continue;
        }
        if (line[0] == 'S') {
            std::istringstream ss(line);
            std::string keyword, header, seq_data;
            ss >> keyword >> header >> seq_data;

            if (seq_data.size() == 0 || seq_data == "*") {
                continue;
            }

            seq = mindex::createSequence();
            seq->data().insert(seq->data().end(), (int8_t *) &seq_data[0], (int8_t *) (&seq_data[0] + seq_data.size()));
            seq->header(header);
            seq->id(id);
            seq->abs_id(abs_id);

            // We found a valid sequence. Yield.
            break;
        }
    }

    return seq;
}

std::string SequenceFileParserGfa1::GetFileHeaderAsString() const {
    // FASTA/FASTQ doesn't have a header. The header_ is initialized
    // in the constructor.
    return header_;
}

int64_t SequenceFileParserGfa1::GetFileOffset() const {
    if (IsOpen() == false) {
        return -1;
    }
    return gztell(fp_gzip_);
}

std::string SequenceFileParserGfa1::GetFilePath() const {
    return path_;
}

bool SequenceFileParserGfa1::IsOpen() const {
    if (path_.empty() || fp_gzip_ == nullptr) {
        return false;
    }
    return true;
}

bool SequenceFileParserGfa1::FileSeek(int64_t abs_pos) {
    if (IsOpen() == false) {
        return false;
    }
    gzseek(fp_gzip_, abs_pos, SEEK_SET);
    return true;
}

bool SequenceFileParserGfa1::Open(const std::string& path) {
    path_ = path;
    fp_gzip_ = gzopen(path_.c_str(), "r");

    if (fp_gzip_ == NULL) {
        FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", path_.c_str());
        return false;
    }

    return true;
}

}
