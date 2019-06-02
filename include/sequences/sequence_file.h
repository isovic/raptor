/*
 * sequence.h
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_FILE2_H_
#define SRC_SEQUENCES_SEQUENCE_FILE2_H_

#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <sequences/sequence.h>
#include <params/params_raptor.h>
#include <sequences/sequence_file_enums.h>
#include <sequences/sequence_file_parser_factory.h>

namespace mindex {

class SequenceFile;

typedef std::shared_ptr<mindex::SequenceFile> SequenceFilePtr;

mindex::SequenceFilePtr createSequenceFile();

class SequenceFile {
public:
    /*
     * Creation.
    */
    friend mindex::SequenceFilePtr createSequenceFile();

    bool LoadAll(const std::vector<std::string>& headers, const std::vector<std::string>& seqs, bool convert_to_uppercase=true);

    const mindex::SequencePtr& GetSeqByID(int64_t id) const;
    const mindex::SequencePtr& GetSeqByAbsID(int64_t abs_id) const;
    void Add(std::unique_ptr<mindex::Sequence> seq);

    bool Serialize(std::ostream& ofs, mindex::SequenceFormat to_fmt);

    /*
     * Getters.
    */
    const std::vector<mindex::SequencePtr>& seqs() const {
        return seqs_;
    }
    int64_t batch_id() const {
        return batch_id_;
    }
    int64_t batch_start_seq_id() const {
        return batch_start_seq_id_;
    }
    int64_t total_size() const {
        return total_size_;
    }
    int64_t size() const {
        return seqs_.size();
    }
    const mindex::SequencePtr& operator[](int64_t id) const {
        return GetSeqByID(id);
    }
    const HeaderGroupType& header_groups() const {
        return header_groups_;
    }
    const HeaderGroupType& GetHeaderGroups() const {
        return header_groups_;
    }

    /*
     * Setters.
    */
    void seqs(std::vector<mindex::SequencePtr>& _seqs) {
        seqs_ = std::move(_seqs);
    }
    void batch_start_seq_id(int64_t val) {
        batch_start_seq_id_ = val;
    }
    void batch_id(int64_t val) {
        batch_id_ = val;
    }

    /*
     * Modifiers.
    */
    std::vector<mindex::SequencePtr>& seqs() {
        return seqs_;
    }

private:
    SequenceFile();
    SequenceFile(const SequenceFile&) = delete;
    SequenceFile& operator=(const SequenceFile&) = delete;

    std::vector<mindex::SequencePtr> seqs_;     // All sequences.
    int64_t batch_id_;              // The ID of the currently loaded batch.
    int64_t batch_start_seq_id_;    // Absolute ID of the first sequence in the loaded batch.
    int64_t total_size_;            // Sum of all sequence lengths.
    mindex::SequencePtr dummy_nullptr_seq_;
    HeaderGroupType header_groups_;
};

}

#endif
