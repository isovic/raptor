/*
 * sequence.h
 *
 *  Created on: Jan 15, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_RANDOM_ACCESS_SEQUENCE_FILE_H_
#define SRC_SEQUENCES_RANDOM_ACCESS_SEQUENCE_FILE_H_

#include <cstdint>
#include <deque>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <index/sequence.h>
#include <params/params_raptor.h>
#include <index/sequence_file_handlers.h>
#include <index/sequence_file_enums.h>
#include <index/sequence_file.h>
#include <index/sequence_deserializer.h>

namespace mindex {

class RandomAccessSequenceFile;

typedef std::shared_ptr<mindex::RandomAccessSequenceFile> RandomAccessSequenceFilePtr;

mindex::RandomAccessSequenceFilePtr createRandomAccessSequenceFile(const std::string& in_path, size_t max_num_streams);

class RaptorDBRecordSequence {
  public:
    RaptorDBRecordSequence()
        : id(0), name(), seq_len(0), file_id(0), data_start(0), data_len(0) {}
    int64_t id;
    std::string name;
    int64_t seq_len;
    int64_t file_id;
    int64_t data_start;
    int64_t data_len;
};

class RaptorDBRecordFile {
  public:
    RaptorDBRecordFile()
        : id(0), path(), format_str(), format(mindex::SequenceFormat::Unknown) {}
    int64_t id;
    std::string path;
    std::string format_str;
    mindex::SequenceFormat format;
};

class RaptorDBRecordBlock {
  public:
    RaptorDBRecordBlock()
        : id(0), seq_id_start(0), seq_id_end(0), bases_in_block(0) {}
    int64_t id;
    int64_t seq_id_start;
    int64_t seq_id_end;
    int64_t bases_in_block;
};

class RandomAccessSequenceFile {
public:
    /*
     * Creation.
    */
    friend mindex::RandomAccessSequenceFilePtr createRandomAccessSequenceFile(const std::string& in_path, size_t max_num_streams);
    bool LoadDB(const std::string& in_path);

    mindex::SequencePtr FetchSequence(const std::string& qname);
    mindex::SequencePtr FetchSequence(int64_t qid);

    mindex::SequenceFilePtr FetchBlock(const int64_t block_id);
    mindex::SequenceFilePtr FetchAll();

    const std::vector<RaptorDBRecordFile>& db_files() const {
        return db_files_;
    }
    const std::vector<RaptorDBRecordSequence>& db_seqs() const {
        return db_seqs_;
    }
    const std::vector<RaptorDBRecordBlock>& db_blocks() const {
        return db_blocks_;
    }

private:
    RandomAccessSequenceFile(size_t max_num_streams);

    RandomAccessSequenceFile(const RandomAccessSequenceFile&) = delete;
    RandomAccessSequenceFile& operator=(const RandomAccessSequenceFile&) = delete;

    mindex::SequencePtr FetchSequence_(int64_t db_seq_id);

    size_t max_num_streams_;
    float db_version_;
    std::vector<RaptorDBRecordFile> db_files_;
    std::vector<RaptorDBRecordSequence> db_seqs_;
    std::vector<RaptorDBRecordBlock> db_blocks_;
    std::unordered_map<int64_t, size_t> seq_id_to_vec_;         // The first is the seq_id from DB, and second is the ordinal number in the vector.
    std::unordered_map<int64_t, size_t> file_id_to_vec_;        // The first is the file_id from DB, and second is the ordinal number in the vector.
    std::unordered_map<int64_t, size_t> block_id_to_vec_;       // The first is the block_id from DB, and second is the ordinal number in the vector.
    std::unordered_map<std::string, size_t> qname_to_vec_;      // The first is the sequence qname from DB, and second is the ordinal number in the vector.
    std::vector<std::unique_ptr<mindex::SequenceFileHandlers>> streams_;
    std::deque<int64_t> fid_stream_priority_;
};

}

#endif
