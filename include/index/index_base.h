/*
 * index_base.h
 *
 *  Created on: Sep 26, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_INDEX_BASE_H_
#define SRC_RAPTOR_INDEX_BASE_H_

#include <cstdint>
#include <string>
#include <vector>
#include <sequences/sequence_file.h>
#include <index/seed_hit.hpp>
#include <index/index_types.h>

namespace mindex {

class IndexBase {
   public:
    virtual void AddSequence(const std::string& seq, const std::string& header) = 0;
    virtual void AddSequence(const std::vector<int8_t>& seq, const std::string& header) = 0;
    virtual void AddSequence(const int8_t* seq, size_t seq_len, const char* header, size_t header_len) = 0;
    virtual bool AddSequences(const std::vector<std::string>& seqs, const std::vector<std::string>& headers) = 0;
    virtual void SetSequenceFile(mindex::SequenceFilePtr seq_file) = 0;
    virtual std::string FetchSeqAsString(size_t seq_id, bool rev_cmp) const = 0;
    virtual std::string FetchSeqAsString(size_t seq_id, size_t start, size_t end, bool rev_cmp) const = 0;
    virtual const int8_t* FetchRawSeq(size_t seq_id) const = 0;
    virtual int BuildIndex() = 0;
    virtual std::vector<mindex::SeedHitPacked> CollectHits(const int8_t* seq, size_t seq_len, indid_t seq_id = 0) = 0;
    virtual std::vector<mindex::SeedHitPacked> CollectHits(const std::string& seq) = 0;
    virtual int Load(const std::string& index_path) = 0;
    virtual int Store(const std::string& index_path) = 0;

    /*
     * Getters.
    */
    virtual const std::shared_ptr<mindex::IndexParams> params() const = 0;
    virtual std::string header(size_t seq_id) const = 0;
    virtual ind_t len(size_t seq_id) const = 0;
    virtual const mindex::SequenceFilePtr& seqs() const = 0;
    virtual int64_t total_len() const = 0;
    virtual int64_t num_seqs() const = 0;
    virtual const std::vector<mindex128_t>& seeds() const = 0;

    // size_t occ_max() const { return occ_max_; }
    // size_t occ_max_after() const { return occ_max_after_; }
    // size_t occ_cutoff() const { return occ_cutoff_; }
    // double occ_singletons() const { return occ_singletons_; }
    // double occ_avg_before() const { return occ_avg_before_; }
    // double occ_avg_after() const { return occ_avg_after_; }
    // double spacing() const { return spacing_; }
};

}  // namespace mindex

#endif
