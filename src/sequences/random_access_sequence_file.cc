/*
 * sequence_file.cc
 *
 *  Created on: Jan 15, 2019
 *      Author: Ivan Sovic
 */

#include <sequences/random_access_sequence_file.h>
#include <utility/stringutil.h>
#include <utility/files.hpp>
#include <log/log_tools.h>

namespace mindex {

mindex::RandomAccessSequenceFilePtr createRandomAccessSequenceFile(const std::string& in_path, size_t max_num_streams) {
    auto ret = mindex::RandomAccessSequenceFilePtr(new mindex::RandomAccessSequenceFile(max_num_streams));
    bool rv = ret->LoadDB(in_path);
    ret->ParseAndMergeHeaderGroups_();
    if (rv == false) {
        return nullptr;
    }
    return std::move(ret);
}

RandomAccessSequenceFile::RandomAccessSequenceFile(size_t max_num_streams)
                :
                    max_num_parsers_(max_num_streams),
                    db_version_(-1.0f), db_files_(), db_seqs_(), db_blocks_(),
                    seq_id_to_vec_(), file_id_to_vec_(), block_id_to_vec_(), qname_to_vec_(),
                    parsers_(), fid_stream_priority_(),
                    header_groups_(),
                    convert_to_uppercase_(true) {
}

bool RandomAccessSequenceFile::LoadDB(const std::string& in_path) {
    std::ifstream ifs(in_path);
    if (ifs.is_open() == false) {
        WARNING_REPORT(ERR_OPENING_FILE, "Could not open file: '%s', returning empty handed.", in_path.c_str());
        return false;
    }

    db_version_ = 0.0f;
    db_files_.clear();
    db_seqs_.clear();
    db_blocks_.clear();
    parsers_.clear();
    fid_stream_priority_.clear();

    std::string line;
    while(std::getline(ifs, line)) {
        if (line.size() == 0) {
            continue;
        }

        std::istringstream iss(line);
        std::string keyword;

        if (line[0] == 'V') {
            iss >> keyword >> db_version_;
        } else if (line[0] == 'F') {
            RaptorDBRecordFile record;
            iss >> keyword >> record.id >> record.path >> record.format_str;
            record.format = mindex::SequenceFormatFromString(record.format_str);
            file_id_to_vec_[record.id] = db_files_.size();
            db_files_.emplace_back(record);
        } else if (line[0] == 'S') {
            RaptorDBRecordSequence record;
            iss >> keyword >> record.id >> record.name >> record.seq_len >> record.file_id >> record.data_start >> record.data_len;
            seq_id_to_vec_[record.id] = db_seqs_.size();
            qname_to_vec_[record.name] = db_seqs_.size();
            db_seqs_.emplace_back(record);
        } else if (line[0] == 'B') {
            RaptorDBRecordBlock record;
            iss >> keyword >> record.id >> record.seq_id_start >> record.seq_id_end >> record.bases_in_block;
            block_id_to_vec_[record.id] = db_blocks_.size();
            db_blocks_.emplace_back(record);
        }
    }

    parsers_ = std::vector<mindex::SequenceFileParserBasePtr>(db_files_.size());

    // Sanity check that all formats are known.
    for (const auto& db_file: db_files_) {
        if (db_file.format == mindex::SequenceFormat::Unknown) {
            ERROR_REPORT(ERR_WRONG_FILE_TYPE, "Unknown format '%s' for file '%s' in RaptorDB: '%s'.",
                                db_file.format_str.c_str(), db_file.path.c_str(), in_path.c_str());
        }
    }

    // std::cerr << "db_version_ = " << db_version_ << "\n";
    // std::cerr << "db_files_.size() = " << db_files_.size() << "\n";
    // std::cerr << "db_seqs_.size() = " << db_seqs_.size() << "\n";
    // std::cerr << "db_blocks_.size() = " << db_blocks_.size() << "\n";

    return true;
}

mindex::SequencePtr RandomAccessSequenceFile::FetchSequence(const std::string& qname) {
    auto it_qname_to_vec = qname_to_vec_.find(qname);
    if (it_qname_to_vec == qname_to_vec_.end()) {
        return nullptr;
    }
    size_t db_seq_id = it_qname_to_vec->second;
    return FetchSequence_(db_seq_id);
}

mindex::SequencePtr RandomAccessSequenceFile::FetchSequence(int64_t qid) {
    auto it_seq_id_to_vec = seq_id_to_vec_.find(qid);
    if (it_seq_id_to_vec == seq_id_to_vec_.end()) {
        return nullptr;
    }
    size_t db_seq_id = it_seq_id_to_vec->second;
    return FetchSequence_(db_seq_id);
}

mindex::SequencePtr RandomAccessSequenceFile::FetchSequence_(int64_t db_seq_id) {
    const auto& db_seq = db_seqs_[db_seq_id];

    auto it_file_id_to_vec = file_id_to_vec_.find(db_seq.file_id);
    if (it_file_id_to_vec == file_id_to_vec_.end()) {
        return nullptr;
    }
    size_t db_file_id = it_file_id_to_vec->second;
    const auto& db_file = db_files_[db_file_id];

    if (parsers_[db_file_id] == nullptr) {
        parsers_[db_file_id] = std::move(mindex::createSequenceFileParser(db_file.path, db_file.format));
        if (parsers_[db_file_id] == nullptr) {
            // If this is still null, then something is wrong.
            ERROR_REPORT(ERR_FILE_NOT_FOUND, "Problem opening file '%s' from RaptorDB.", db_file.path.c_str());
        }
        fid_stream_priority_.push_back(db_file_id);
    }

    while (fid_stream_priority_.size() > max_num_parsers_) {
        parsers_[fid_stream_priority_.front()] = nullptr;
        fid_stream_priority_.pop_front();
    }

    if (parsers_[db_file_id] == nullptr) {
        return nullptr;
    }

    parsers_[db_file_id]->FileSeek(db_seq.data_start);

    mindex::SequencePtr ret = parsers_[db_file_id]->YieldSequence(convert_to_uppercase_);

    if (ret == nullptr) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Deserialized sequence is nullptr! db_seq_id = %ld, db_seq.name = %s, db_seq.data_start = %ld, db_seq.seq_len = %ld\n", db_seq_id, db_seq.name.c_str(), db_seq.data_start, db_seq.seq_len);
        return nullptr;
    }

    ret->abs_id(db_seq.id);
    ret->id(db_seq.id);

    // std::cerr << db_seq.id << "\t" << db_seq.name << "\t" << db_seq.seq_len << "\n";
    // std::cerr << db_file.id << "\t" << db_file.path << "\t" << db_file.format_str << "\n";
    // std::cerr << "streams_.size() = " << streams_.size() << "\n";
    // std::cerr << "fid_stream_priority_.size() = " << fid_stream_priority_.size() << "\n";
    // std::cerr << ">" << ret->header() << "\n" << ret->GetSequenceAsString() << "\n";

    return std::move(ret);
}

mindex::SequenceFilePtr RandomAccessSequenceFile::FetchBlock(const int64_t block_id) {
    mindex::SequenceFilePtr seq_file = mindex::createSequenceFile();

    if (block_id < 0 || block_id >= static_cast<int64_t>(db_blocks().size())) {
        WARNING_REPORT(
            ERR_UNEXPECTED_VALUE,
            "Invalid block_id value. "
            "block_id = %ld, db_blocks().size() = %lu. "
            "Returning an empty SequenceFile.\n",
            block_id, db_blocks().size());
        return seq_file;
    }

    int64_t start_id = db_blocks()[block_id].seq_id_start;
    int64_t end_id = db_blocks()[block_id].seq_id_end;

    seq_file->batch_start_seq_id(start_id);

    for (int64_t i = start_id; i < end_id; ++i) {
        // We have to use the public interface here, in case sequences are out
        // of order after loading (order of appearance can be different than
        // the naming order).
        // The public interface does the ID lookup first, and then fetches
        // the sequecne.
        int64_t seq_id = db_seqs()[i].id;
        mindex::SequencePtr seq = FetchSequence(seq_id);
        if (seq == nullptr) {
            WARNING_REPORT(
                ERR_UNEXPECTED_VALUE,
                "Could not load sequence from RaptorDB, seq_id = %ld. "
                "Returning an empty SequenceFile.\n", seq_id);
            return mindex::createSequenceFile();
        }
        seq_file->Add(std::move(seq));
    }

    return seq_file;
}

mindex::SequenceFilePtr RandomAccessSequenceFile::FetchAll() {
    mindex::SequenceFilePtr seq_file = mindex::createSequenceFile();

    int64_t start_id = 0;
    int64_t end_id = db_seqs().size();

    for (int64_t i = start_id; i < end_id; ++i) {
        // The FetchSequence_ doesn't perform the lookup first, so
        // it's faster this way than to use the public interface.
        mindex::SequencePtr seq = FetchSequence_(i);
        if (seq == nullptr) {
            WARNING_REPORT(
                ERR_UNEXPECTED_VALUE,
                "Could not load sequence from RaptorDB, seq_id = %ld. "
                "Returning an empty SequenceFile.\n", db_seqs()[i].id);
            return mindex::createSequenceFile();
        }
        seq_file->Add(std::move(seq));
    }

    return seq_file;
}

bool RandomAccessSequenceFile::ParseAndMergeHeaderGroups_() {
    mindex::HeaderGroupType merged_groups;

    for (size_t fid = 0; fid < db_files_.size(); ++fid) {

        const std::string& in_path = db_files_[fid].path;
        mindex::SequenceFormat& in_fmt = db_files_[fid].format;

        if (in_fmt == mindex::SequenceFormat::Auto) {
            auto ext = raptor::GetFileExt(in_path);
            // If the file is Gzipped, the .gz will be in the ext.
            // E.g. the output from GetFileExt can be "fasta.gz".
            if (ext.size() >= 3 && ext.substr(ext.size() - 3) == ".gz") {
                ext = ext.substr(0, ext.size() - 3);
            }
            in_fmt = SequenceFormatFromString(ext);
        }

        auto parser = mindex::createSequenceFileParser(in_path, in_fmt);
        const auto& curr_groups = parser->GetHeaderGroups();

        for (const auto& it_field: curr_groups) {
            for (const auto& it_ids: it_field.second) {
                merged_groups[it_field.first][it_ids.first] = it_ids.second;
            }
        }
    }

    std::swap(header_groups_, merged_groups);

    return true;
}

}
