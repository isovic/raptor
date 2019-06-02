/*
 * sequence_file.cc
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#include <sequences/sequence_file.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>

namespace mindex {

mindex::SequenceFilePtr createSequenceFile() {
    return mindex::SequenceFilePtr(new mindex::SequenceFile());
}

SequenceFile::SequenceFile()
                :
                    seqs_(),
                    batch_id_(-1),
                    batch_start_seq_id_(0),
                    total_size_(0),
                    dummy_nullptr_seq_(nullptr),
                    header_groups_()
{

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

}
