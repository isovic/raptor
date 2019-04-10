/*
 * sequence.cc
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#include <index/sequence.h>
#include <index/lookups.h>
#include <sstream>
#include <string>

namespace mindex {

mindex::SequencePtr createSequence() {
    return mindex::SequencePtr(new mindex::Sequence());
}

mindex::SequencePtr createSequence(const std::string& header, const seq_t* seq, size_t seq_len) {
    return mindex::SequencePtr(new mindex::Sequence(header, seq, seq_len));
}

mindex::SequencePtr createSequence(const std::string& header, const std::string& seq) {
    return mindex::SequencePtr(new mindex::Sequence(header, seq));
}

mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data) {
    return mindex::SequencePtr(new mindex::Sequence(header, data));
}

mindex::SequencePtr createSequence(const std::string& header, const std::string& seq, const std::string& qual) {
    return mindex::SequencePtr(new mindex::Sequence(header, seq, qual));
}

mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data, const std::vector<mindex::seqqual_t>& qual) {
    return mindex::SequencePtr(new mindex::Sequence(header, data, qual));
}

Sequence::Sequence()
                        : header_(), data_(), qual_(), id_(0), abs_id_(0), header_hash_(0) {

}

Sequence::Sequence(const std::string& _header, const seq_t* seq, size_t seq_len)
                        : header_(), data_(seq, seq + seq_len), qual_(), id_(0), abs_id_(0), header_hash_(0) {
    // Call the setter, which calculates the hash.
    header(_header);
}

Sequence::Sequence(const std::string& _header, const std::string& seq)
                        : header_(), data_(seq.begin(), seq.end()), qual_(), id_(0), abs_id_(0), header_hash_(0) {
    // Call the setter, which calculates the hash.
    header(_header);
}

Sequence::Sequence(const std::string& _header, const std::vector<mindex::seq_t>& data)
                        : header_(), data_(data), qual_(), id_(0), abs_id_(0), header_hash_(0) {
    // Call the setter, which calculates the hash.
    header(_header);
}

Sequence::Sequence(const std::string& _header, const std::string& seq, const std::string& qual)
                        : header_(), data_(seq.begin(), seq.end()), qual_(qual.begin(), qual.end()), id_(0), abs_id_(0), header_hash_(0) {
    // Call the setter, which calculates the hash.
    header(_header);
}

Sequence::Sequence(const std::string& _header, const std::vector<mindex::seq_t>& data, const std::vector<mindex::seqqual_t>& qual)
                        : header_(), data_(data), qual_(qual), id_(0), abs_id_(0), header_hash_(0) {
    // Call the setter, which calculates the hash.
    header(_header);
}

void Sequence::ReverseComplement() {
    size_t max_i = data_.size() / 2;
    for (size_t i = 0, j = (data_.size() - 1); i < max_i; i++, j--) {
        data_[i] = mindex::nuc_to_complement[data_[i]];
        data_[j] = mindex::nuc_to_complement[data_[j]];
        std::swap(data_[i], data_[j]);
    }
    if ((data_.size() & 0x01) != 0) {
        data_[max_i] = mindex::nuc_to_complement[data_[max_i]];
    }
}

void Sequence::ToUppercase() {
    for (size_t i = 0; i < data_.size(); ++i) {
        data_[i] = std::toupper(data_[i]);
    }
}

std::string Sequence::GetTrimmedHeader() const {
    std::istringstream iss(header());
    std::string trimmed;
    iss >> trimmed;
    return trimmed;
}

}
