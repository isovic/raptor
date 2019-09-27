/*
 * dense_index.h
 *
 *  Created on: Sep 27, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_INDEX_DENSE_INDEX_H_
#define SRC_INDEX_DENSE_INDEX_H_

#include <cstdint>
#include <deque>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include <index/index_base.h>
#include <index/index_params.h>
#include <index/index_types.h>
#include <index/seed.hpp>
#include <index/seed_hit.hpp>
#include <sequences/sequence_file.h>

namespace mindex {

class DenseIndex;

std::unique_ptr<mindex::IndexBase> createDenseIndex(std::shared_ptr<mindex::IndexParams> params);

class DenseIndex : mindex::IndexBase {
   public:
    friend std::unique_ptr<mindex::IndexBase> createDenseIndex(
        std::shared_ptr<mindex::IndexParams> params);

    ~DenseIndex();

    void AddSequence(const std::string& seq, const std::string& header) override;

    void AddSequence(const std::vector<int8_t>& seq, const std::string& header) override;

    void AddSequence(const int8_t* seq, size_t seq_len, const char* header, size_t header_len) override;

    bool AddSequences(const std::vector<std::string>& seqs, const std::vector<std::string>& headers) override;

    void SetSequenceFile(mindex::SequenceFilePtr seq_file) override;

    std::string FetchSeqAsString(size_t seq_id, bool rev_cmp) const override;

    std::string FetchSeqAsString(size_t seq_id, size_t start, size_t end, bool rev_cmp) const override;

    const int8_t* FetchRawSeq(size_t seq_id) const override;

    int BuildIndex() override;

    std::vector<mindex::SeedHitPacked> CollectHits(const int8_t* seq, size_t seq_len,
                                                    indid_t seq_id = 0) override;

    std::vector<mindex::SeedHitPacked> CollectHits(const std::string& seq) override;

    int Load(const std::string& index_path) override;

    int Store(const std::string& index_path) override;

    /*
     * Getters.
    */
    const std::shared_ptr<mindex::IndexParams> params() const override { return params_; }

    // const std::vector<int8_t>& data() const { return data_; }

    std::string header(size_t seq_id) const override {
        const mindex::SequencePtr& seq = seqs()->GetSeqByID(seq_id);
        if (seq == nullptr) {
            std::cerr << "[DenseIndex::header] Warning: seq_id does not match a sequence in index! Returning empty.";
            return std::string();
        }
        return seq->header();
    }

    // const std::vector<int64_t>& starts() const { return starts_; }

    ind_t len(size_t seq_id) const override {
        const mindex::SequencePtr& seq = seqs()->GetSeqByID(seq_id);
        if (seq == nullptr) {
            std::cerr << "[DenseIndex::len] Warning: seq_id does not match a sequence in index! Returning empty.";
            return 0;
        }
        return seq->len();
    }

    void seqs(mindex::SequenceFilePtr val) override { seqs_ = std::move(val); }

    const mindex::SequenceFilePtr& seqs() const override { return (seqs_ == nullptr) ? dummy_nullptr_seqs_ : seqs_; }

    int64_t total_len() const override { return (seqs_ == nullptr) ? 0 : seqs_->total_size(); }

    int64_t num_seqs() const override { return (seqs_ == nullptr) ? 0 : seqs_->size(); }

    const std::vector<mindex128_t>& seeds() const override { return seeds_; }

    int32_t k() const {
        if (params_ != nullptr) {
            return params_->k;
        }
        return 0;
    }

    int32_t w() const {
        if (params_ != nullptr) {
            return params_->w;
        }
        return 0;
    }

    size_t occ_max() const { return occ_max_; }

    size_t occ_max_after() const { return occ_max_after_; }

    size_t occ_cutoff() const { return occ_cutoff_; }

    double occ_singletons() const { return occ_singletons_; }

    double occ_avg_before() const { return occ_avg_before_; }

    double occ_avg_after() const { return occ_avg_after_; }

    double spacing() const { return spacing_; }

   private:
    DenseIndex(std::shared_ptr<mindex::IndexParams> params);
    DenseIndex(const DenseIndex&) = delete;
    DenseIndex& operator=(const DenseIndex&) = delete;

    /*
     * For a given sequence, generates and returns a set of minimizers.
     * The `minimizers` vector is never cleared but only inserted to, so minimizers
     * can be accumulated through multiple calls to this method.
     * This method takes care to not generate the minimizers across the non-ACTG regions.
     * The seq_offset is the distance from the beginning of the sequence, used for the seed position.
     *
    */
    static int GenerateMinimizers_(std::vector<mindex128_t>& minimizers,
                                    const int8_t* seq, ind_t seq_len, ind_t seq_offset,
                                    indid_t seq_id, int32_t k, int32_t w, bool use_rc,
                                    bool homopolymer_suppression,
                                    int32_t max_homopolymer_run,
                                    ind_t seq_start, ind_t seq_end);

    /*
     * Takes the vector of minimizers, initializes bucket for hashes, and constructs the hashes.
     */
    void ConstructHash_(const std::vector<mindex128_t>& minimizers);

    /*
     * Used by ConstructHash_ to precalculate the number of minimizers which fall within
     * each of the buckets, so that the size of individual hashes can be initialized ahead
     * of filling the hashing.
     */
    static std::vector<int64_t> BucketPrecount_(const std::vector<mindex128_t>& minimizers,
                                                int32_t num_buckets, int32_t bucket_shift,
                                                int32_t bucket_mask, int64_t& total_key_count);

    /*
     * Counts occurrence of each key in seeds_, sorts them and returns the statistics
     * which include the maximum occurrence, the occurrence cutoff (based on freq_percentil),
     * count of singleton keys (keys with only one hit) and the average occurrence.
     */
    void CalcOccurrenceThreshold_(size_t& occ_max, size_t& occ_cutoff, size_t& occ_singletons,
                                  double& occ_avg) const;

    // /*
    //  * This method linearly processes the minimizers by skipping those with occurrence
    //  * count abouve the cutoff. It shifts the minimizers towards left as much as possible,
    //  * and resizes the minimizers vector to fit.
    //  */
    // void RemoveFrequentMinimizers_(std::vector<mindex128_t>& minimizers, size_t cutoff) const;

    /*
     * Linear pass over minimizers to count the number of distinct keys
     * present in the array.
     */
    size_t CountKeys_(const std::vector<mindex128_t>& minimizers) const;

    /*
     * A helper method to determine the number of bits to shift a key, so
     * that the remaining bits are used as the hash bucket ID.
     * TODO: This needs to be reimplemented differently.
     */
    int32_t CalcBucketShift_() const;

    int32_t CalcBucketMask_() const;

    /*
     * For a given position on the sequence, this method finds the first following
     * position which is not part of a homopolymer. If the homopolymer is above
     * the specified maximum length, it returns false, otherwise it returns true.
     * This method is basically a wrapper method used when building kmers for
     * indexing, in case homopolymer suppresion is specified.
    */
    bool SkipHP_(const int8_t* seq, size_t seq_len, int32_t pos, int32_t max_homopolymer_len, int32_t& ret_pos) const;

    /*
     * For every position of the given sequence, it performs the same procedure
     * as when generating the minimizers. Instead of generating minimizers, this
     * method builds a relation map between the kmer start and end position,
     * so that reverse complement positions (when collecting hits) can be compensated.
    */
    std::vector<int32_t> CalcHPSpan_(const int8_t* seq, size_t seq_len, int32_t k, bool homopolymer_suppression, int32_t max_homopolymer_len) const;

    /*
     * FindKey_ wraps lookup of the key in the hashes. Since there are
     * Multiple hashes used, this method locates the correct one and
     * attempts to return an iterator to the queried element.
     */
    inline SeedHashType2::iterator FindKey_(const minkey_t& key, bool& is_found) {
        // Decode the hash bucket.
        int32_t bucket = (key >> bucket_shift_) & bucket_mask_;

        // Sanity check.
        if (bucket >= hashes_.size()) {
            fprintf (stderr, "ERROR: Determined bucket exceeds the size of hashes_. bucket = %d, "
                         "hashes_.size() = %lu\n",
                         bucket, hashes_.size());
        }

        // Find and add all the key hits.
        auto it = hashes_[bucket].find(key);

        // Since there are multiple possible hashes that were looked-up,
        // we provide a special bool value so that a user can
        // verify if the key was found.
        is_found = (it != hashes_[bucket].end());

        return it;
    }

    static inline minkey_t InvertibleHash_(minkey_t key) {
        /*
        Credit: Heng Li, Minimap2.
        */
        key = (~key + (key << 21));
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8));
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4));
        key = key ^ key >> 28;
        key = (key + (key << 31));
        return key;
    }

    void VerboseSeeds_(std::ostream& os);

    std::shared_ptr<mindex::IndexParams> params_;

    mindex::SequenceFilePtr seqs_;

    // This is used as a hack for a const reference getter.
    mindex::SequenceFilePtr dummy_nullptr_seqs_;

    // Statistics on seed counts.
    size_t occ_max_;
    size_t occ_max_after_;
    size_t occ_cutoff_;
    double occ_singletons_;
    double occ_avg_before_;
    double occ_avg_after_;
    size_t num_keys_;
    double spacing_;

    int32_t bucket_shift_;
    int32_t bucket_mask_;

    std::vector<mindex128_t> seeds_;
    std::vector<SeedSpan> hash_ranges_;
    std::vector<SeedHashType2> hashes_;
};

}  // namespace mindex

#endif
