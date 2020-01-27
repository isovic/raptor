/*
 * sequence_file_parser_fastx.cc
 *
 *  Created on: Apr 28, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file_parser_fastx.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>

namespace mindex {

mindex::SequenceFileParserBasePtr createSequenceFileParserFastx(const std::string& in_path) {
    auto ret = mindex::SequenceFileParserBasePtr(new mindex::SequenceFileParserFastx());
    bool rv = ret->Open(in_path);
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

SequenceFileParserFastx::SequenceFileParserFastx()
    : path_(), fp_kseq_(NULL), fp_gzip_(NULL), header_groups_() {
}

SequenceFileParserFastx::~SequenceFileParserFastx() {
    if (fp_kseq_ != NULL) {
        kseq_destroy(fp_kseq_);
        fp_kseq_ = NULL;
    }

    if (fp_gzip_ != NULL) {
        gzclose(fp_gzip_);
        fp_gzip_ = NULL;
    }
}

SequencePtr SequenceFileParserFastx::YieldSequence() {
    int64_t id = -1;
    int64_t abs_id = -1;

    int32_t l = kseq_read(fp_kseq_);

    if (l < 0) {
        return nullptr;
    }

    auto seq = mindex::createSequence();

    std::string header;
    if (fp_kseq_->name.l > 0) {
        header += std::string(fp_kseq_->name.s);
    }
    if (fp_kseq_->comment.l > 0) {
        header += std::string(" ") + std::string(fp_kseq_->comment.s);
    }

    seq->data().resize(fp_kseq_->seq.l);
    memcpy(&seq->data()[0], reinterpret_cast<int8_t *>(fp_kseq_->seq.s), fp_kseq_->seq.l);

    if (fp_kseq_->qual.l > 0) {
        seq->qual().resize(fp_kseq_->qual.l);
        memcpy(&seq->qual()[0], reinterpret_cast<int8_t *>(fp_kseq_->qual.s), fp_kseq_->qual.l);
    }

    seq->header(header);
    seq->id(id);
    seq->abs_id(abs_id);

    return seq;
}

std::string SequenceFileParserFastx::GetFileHeaderAsString() const {
    // FASTA/FASTQ doesn't have a header. The header_ is initialized
    // in the constructor.
    return {};
}

int64_t SequenceFileParserFastx::GetFileOffset() const {
    if (IsOpen() == false) {
        return -1;
    }

    int64_t ret_val = gztell(fp_gzip_);

    if (ret_val < 0) {
        return ret_val;
    }
    if (fp_kseq_ == NULL) {
        return -5;
    }

    int64_t offset = static_cast<int64_t>(fp_kseq_->f->end) - static_cast<int64_t>(fp_kseq_->f->begin);

    if (offset < 0) {
        return -6;
    }

    // The "-1" because kseq parser loads the next '>' or '@' character too.
    int64_t offset_separator_char = (fp_kseq_->last_char == '>' || fp_kseq_->last_char == '@') ? (1) : 0;
    ret_val -= (offset + offset_separator_char);

    return ret_val;
}

const HeaderGroupType& SequenceFileParserFastx::GetHeaderGroups() const {
    return header_groups_;
}

std::string SequenceFileParserFastx::GetFilePath() const {
    return path_;
}

bool SequenceFileParserFastx::IsOpen() const {
    if (path_.empty() || fp_kseq_ == nullptr || fp_gzip_ == nullptr) {
        return false;
    }
    return true;
}

bool SequenceFileParserFastx::FileSeek(int64_t abs_pos) {
    if (IsOpen() == false) {
        return false;
    }
    gzseek(fp_gzip_, abs_pos, SEEK_SET);
    kseq_destroy(fp_kseq_);
    fp_kseq_ = kseq_init(fp_gzip_);
    return true;
}

bool SequenceFileParserFastx::Open(const std::string& path) {
    path_ = path;
    fp_gzip_ = gzopen(path_.c_str(), "r");

    if (fp_gzip_ == NULL) {
        FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", path_.c_str());
        return false;
    }

    fp_kseq_ = kseq_init(fp_gzip_);
    return true;
}

}
