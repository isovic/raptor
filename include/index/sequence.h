/*
 * sequence.h
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_H_
#define SRC_SEQUENCES_SEQUENCE_H_

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <functional>
#include <containers/region/region_aligned.h>

#ifdef RAPTOR_COMPILED_WITH_PBBAM
#include <pbbam/BamRecord.h>
#endif

namespace mindex {

// Abstract the underlying data type for possible future updates.
typedef int8_t seq_t;
typedef int8_t seqqual_t;

class Sequence;

typedef std::unique_ptr<mindex::Sequence> SequencePtr;

mindex::SequencePtr createSequence();
mindex::SequencePtr createSequence(const std::string& header, const seq_t* seq, size_t seq_len);
mindex::SequencePtr createSequence(const std::string& header, const std::string& seq);
mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data);
mindex::SequencePtr createSequence(const std::string& header, const std::string& seq, const std::string& qual);
mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data, const std::vector<mindex::seqqual_t>& qual);

class Sequence {
public:
    /*
     * Creation.
    */
    friend mindex::SequencePtr createSequence();
    friend mindex::SequencePtr createSequence(const std::string& header, const seq_t* seq, size_t seq_len);
    friend mindex::SequencePtr createSequence(const std::string& header, const std::string& seq);
    friend mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data);
    friend mindex::SequencePtr createSequence(const std::string& header, const std::string& seq, const std::string& qual);
    friend mindex::SequencePtr createSequence(const std::string& header, const std::vector<mindex::seq_t>& data, const std::vector<mindex::seqqual_t>& qual);

    void ReverseComplement();
    void ToUppercase();

    std::string GetTrimmedHeader() const;

    /*
     * Const getters.
    */
    const std::string& header() const {
        return header_;
    }
    const std::vector<mindex::seq_t>& data() const {
        return data_;
    }
    const std::vector<mindex::seqqual_t>& qual() const {
        return qual_;
    }
    int64_t id() const {
        return id_;
    }
    int64_t abs_id() const {
        return abs_id_;
    }
    int64_t len() const {
        return data_.size();
    }
    size_t header_hash() const {
        return header_hash_;
    }

    /*
     * Modifier getters.
     * Useful to prevent large memory allocations.
    */
    std::string& header() {
        return header_;
    }
    std::vector<mindex::seq_t>& data() {
        return data_;
    }
    std::vector<mindex::seqqual_t>& qual() {
        return qual_;
    }
    const std::vector<raptor::SamTag>& tags() const {
        return tags_;
    }
#ifdef RAPTOR_COMPILED_WITH_PBBAM
    const std::unique_ptr<PacBio::BAM::BamRecord>& apriori_bam() const {
        return apriori_bam_;
    }
#endif

    /*
     * Setters.
    */
    void header(const std::string& _header) {
        header_ = _header;
        header_hash_ = std::hash<std::string>{}(_header);
    }
    void data(const std::vector<mindex::seq_t>& _data) {
        data_ = _data;
    }
    void qual(const std::vector<mindex::seqqual_t>& _qual) {
        qual_ = _qual;
    }
    void id(int64_t _id) {
        id_ = _id;
    }
    void abs_id(int64_t _abs_id) {
        abs_id_ = _abs_id;
    }
    void tags(const std::vector<raptor::SamTag>& val) {
        tags_ = val;
    }
    void AddTag(const raptor::SamTag& val) {
        tags_.emplace_back(val);
    }
#ifdef RAPTOR_COMPILED_WITH_PBBAM
    void apriori_bam(std::unique_ptr<PacBio::BAM::BamRecord> val) {
        apriori_bam_ = std::move(val);
    }
#endif

    std::string GetSequenceAsString() const {
        return std::string((const char *) data_.data(), data_.size());
    }

    std::string GetSubsequenceAsString(int32_t start, int32_t end) const {
        if (start < 0 || start > data_.size() || end < start || end < 0 || end > data_.size()) {
            return std::string();
        }
        return std::string((const char *) data_.data() + start, (end - start));
    }

    std::string GetQualityAsString() const {
        return std::string((const char *) qual_.data(), qual_.size());
    }

private:
    Sequence();
    Sequence(const std::string& header, const seq_t* seq, size_t seq_len);
    Sequence(const std::string& header, const std::string& seq);
    Sequence(const std::string& header, const std::vector<mindex::seq_t>& data);
    Sequence(const std::string& header, const std::string& seq, const std::string& qual);
    Sequence(const std::string& header, const std::vector<mindex::seq_t>& data, const std::vector<mindex::seqqual_t>& qual);

    Sequence(const Sequence&) = delete;
    Sequence& operator=(const Sequence&) = delete;

    std::string header_;
    std::vector<mindex::seq_t> data_;
    std::vector<mindex::seqqual_t> qual_;
    int64_t id_;
    int64_t abs_id_;
    size_t header_hash_;
    // std::unordered_map<std::string, int32_t> tag_id_;
    std::vector<raptor::SamTag> tags_;

#ifdef RAPTOR_COMPILED_WITH_PBBAM
    std::unique_ptr<PacBio::BAM::BamRecord> apriori_bam_;
#endif
};

}

#endif
