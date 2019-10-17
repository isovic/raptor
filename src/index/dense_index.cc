/*
 * dense_index.cc
 *
 *  Created on: Sep 27, 2019
 *      Author: Ivan Sovic
 */

#include <index/dense_index.h>

#include <ctime>
#include <cmath>
#include <iostream>

#include <index/kxsort.h>
#include <index/lookups.h>
#include <utility/trim.hpp>
#include <utility/memtime.h>
#include <utility/tictoc.h>

#include <log/log_system.h>

static const int32_t NUM_HASH_BUCKET_BITS = 10;
static const int32_t NUM_HASH_BUCKETS = pow(2, NUM_HASH_BUCKET_BITS);

namespace mindex {

std::unique_ptr<mindex::IndexBase> createDenseIndex(
    std::shared_ptr<mindex::IndexParams> params) {
    auto ret = std::unique_ptr<mindex::IndexBase>(new DenseIndex(params));
    ret->seqs(mindex::createSequenceFile());
    return std::move(ret);
}

DenseIndex::DenseIndex(std::shared_ptr<mindex::IndexParams> params)
    : params_(params),
      seqs_(nullptr),
      dummy_nullptr_seqs_(nullptr),
      occ_max_(0),
      occ_max_after_(0),
      occ_cutoff_(0),
      occ_singletons_(0.0),
      occ_avg_before_(0),
      occ_avg_after_(0),
      num_keys_(0),
      spacing_(0.0),
      bucket_shift_(0),
      bucket_mask_(0x3FF) {

    hashes_.resize(NUM_HASH_BUCKETS);

    for (auto& hash : hashes_) {
        #ifdef MINIMIZER_INDEX2_USING_DENSEHASH
                hash.set_empty_key(
                    MINIMIZER_INDEX_EMPTY_HASH_KEY);  // Densehash requires this to be defined on top.
        #endif
    }

    bucket_shift_ = CalcBucketShift_();
    bucket_mask_ = CalcBucketMask_();
}

DenseIndex::~DenseIndex() {}

void DenseIndex::AddSequence(const std::string& seq, const std::string& header) {
    if (seq.size() == 0) {
        return;
    }
    AddSequence((const int8_t*)&seq[0], seq.size(), (const char*)&header[0], header.size());
}

void DenseIndex::AddSequence(const std::vector<int8_t>& seq, const std::string& header) {
    if (seq.size() == 0) {
        return;
    }
    AddSequence((const int8_t*)&seq[0], seq.size(), (const char*)&header[0], header.size());
}

bool DenseIndex::AddSequences(const std::vector<std::string>& seqs, const std::vector<std::string>& headers) {
    if (seqs.size() != headers.size()) {
        std::cerr << "Warning: the seqs and the headers vectors are not of same length. Cannot add to index. Skipping." << std::endl;
        return false;
    }
    for (size_t i = 0; i < seqs.size(); i++) {
        AddSequence(seqs[i], headers[i]);
    }
    return true;
}

void DenseIndex::AddSequence(const int8_t* seq, size_t seq_len, const char* header,
                                 size_t header_len) {
    mindex::SequencePtr new_seq = mindex::createSequence(std::string(header, header_len), seq, seq_len);
    if (seqs_ == nullptr) {
        seqs_ = mindex::createSequenceFile();
    }
    seqs_->Add(std::move(new_seq));
}

void DenseIndex::SetSequenceFile(mindex::SequenceFilePtr seq_file) {
    seqs_ = seq_file;
}

int DenseIndex::BuildIndex() {
    seeds_.clear();
    if (seqs_ == nullptr) {
        return 1;
    }

    // An empiric estimate to try to avoid reallocation.
    seeds_.reserve(total_len() * 2.3 / params_->w);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
        LOG_ALL("Building the index.\n");
        LOG_ALL("  - num_seqs() = %llu\n", num_seqs());
        LOG_ALL("  - total_len() = %llu\n", total_len());
        // LOG_ALL("  - data_.capacity() = %llu\n", data_.size() );
        LOG_ALL("  - seeds_.size() = %llu\n", seeds_.size());
        LOG_ALL("  - seeds_.capacity() = %llu\n", seeds_.capacity());
    #endif

    int64_t seed_spacing_sum = 0;
    ind_t region_rstart = 0;
    ind_t region_rend = -1;
    if (params_->is_region_specified) {
        region_rstart = params_->region_rstart;
        region_rend = params_->region_rend;
    }

    // Generate minimizers for each separate seq.
    clock_t time_collect_minimizers_before = clock();
    for (size_t i = 0; i < seqs_->size(); i++) {
        const mindex::SequencePtr& seq = seqs_->GetSeqByID(i);
        if (seq == nullptr) {
            continue;
        }

        if (seq->len() < params_->min_tlen) {
            continue;
        }

        if (params_->is_region_specified) {
            if (raptor::TrimToFirstWhitespace(seq->header()) != params_->region_rname) {
                continue;
            }
        }

        #ifdef RAPTOR_TESTING_MODE
                if (seqs_->size() < 30 || (seqs_->size() >= 30 && (i % 1000) == 0)) {
                    LOG_ALL("[%llu] Generating minimizers for sequence '%s', len = %lld, abs_id = %lld\n", i, seq->header().c_str(), seq->len(), seq->abs_id());
                    LOG_ALL("    Indexing specific region: region_rname = '%s', region_rstart = %lld, region_rend = %lld\n", params_->region_rname.c_str(), region_rstart, region_rend);
                    LOG_ALL("    (before) seeds_.size() = %llu, seeds_.capacity() = %llu\n", seeds_.size(), seeds_.capacity());
                    LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
                }
        #endif

        size_t num_seeds_before = seeds_.size();

        GenerateSeeds_(seeds_, &seq->data()[0], seq->len(), 0, i, params_->k, params_->w,
                            !params_->index_only_fwd_strand,
                            params_->homopolymer_suppression, params_->max_homopolymer_len,
                            region_rstart, region_rend);

        // This part is used to calculate average seed spacing.
        // Useful info, not used for anything functionall atm.
        ind_t prev_pos =
            (num_seeds_before > 0) ? (mindex::Seed::DecodePos(seeds_[num_seeds_before])) : 0;

        for (size_t seed_id = (num_seeds_before + 1); seed_id < seeds_.size(); ++seed_id) {
            ind_t pos = mindex::Seed::DecodePos(seeds_[seed_id]);
            seed_spacing_sum += (pos - prev_pos);
            prev_pos = pos;
        }

        #ifdef RAPTOR_TESTING_MODE
                if (seqs_->size() < 30) {
                    LOG_ALL("    (after) seeds_.size() = %llu, seeds_.capacity() = %llu\n", seeds_.size(), seeds_.capacity());
                    LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
                }
        #endif
    }
    spacing_ = (seeds_.size() > 0) ? (((double)seed_spacing_sum) / ((double)seeds_.size())) : 0.0;
    clock_t time_collect_minimizers_after = clock();
    double time_collect_minimizers = ((double)(time_collect_minimizers_after - time_collect_minimizers_before)) / CLOCKS_PER_SEC;
    LOG_ALL("Collected %llu minimizers in %.2f seconds.\n", seeds_.size(), time_collect_minimizers);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Average minimizer spacing: %lf\n", spacing_);
        LOG_ALL("Sorting minimizers.\n");
    #endif

    clock_t before_sort = clock();
    // std::sort(seeds_.begin(), seeds_.end());
    kx::radix_sort(seeds_.begin(), seeds_.end());
    clock_t after_sort = clock();
    double sort_time = ((double)(after_sort - before_sort)) / CLOCKS_PER_SEC;
    LOG_ALL("Sorted %llu minimizers in %.2f seconds.\n", seeds_.size(), sort_time);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
        LOG_ALL("Counting keys.\n");
    #endif

    num_keys_ = CountKeys_(seeds_);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("There are %lld minimizer keys.\n", num_keys_);
    #endif

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Calculating occurrence thresholds.\n");
    #endif

    size_t occ_num_singletons = 0;
    CalcOccurrenceThreshold_(occ_max_, occ_cutoff_, occ_num_singletons, occ_avg_before_);
    occ_singletons_ = ((double)occ_num_singletons) / ((double)seeds_.size());

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("From CalcOccurrenceThreshold_: occ_max =  %lld, occ_cutoff = %lld, occ_num_singletons = %lld, occ_avg_before_ = %lld\n", occ_max_, occ_cutoff_, occ_num_singletons, occ_avg_before_);
        LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
    #endif

    // #ifdef RAPTOR_TESTING_MODE
    //     LOG_ALL("Removing frequent minimizers.\n");
    // #endif

    //     RemoveFrequentMinimizers_(seeds_, occ_cutoff_);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Constructing the hash.\n");
        fflush(stderr);
    #endif

    clock_t before_hash = clock();
    ConstructHash_(seeds_);
    clock_t after_hash = clock();
    double hash_time = ((double)(after_hash - before_hash)) / CLOCKS_PER_SEC;
    LOG_ALL("Constructed the hash in %.2f seconds.\n", hash_time);

    #ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Finally, there are %lld keys, and %llu minimizer seeds.\n", num_keys_, seeds_.size());
    #endif

    {
        size_t occ_cutoff = 0, occ_singletons = 0;

        CalcOccurrenceThreshold_(occ_max_after_, occ_cutoff, occ_singletons, occ_avg_after_);

        #ifdef RAPTOR_TESTING_MODE
                LOG_ALL("Finally, from CalcOccurrenceThreshold_: occ_max =  %lld, occ_cutoff = %lld, occ_num_singletons = %lld, occ_avg_before_ = %lld\n", occ_max_after_, occ_cutoff_, occ_num_singletons, occ_avg_after_);
                LOG_ALL("Memory usage: %.2f GB\n", ((double)raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
        #endif
    }

	LOG_ALL("Seed statistics: avg = %.2f, max_before = %lu, max_after = %lu, cutoff = %lu, singletons = %.2f%%, spacing = %.2f\n",
                occ_avg_before(), occ_max(), occ_max_after(), occ_cutoff(), 100.0 * occ_singletons(), spacing());

    return 0;
}

std::vector<int64_t> DenseIndex::BucketPrecount_(const std::vector<mindex128_t>& minimizers,
                                                     int32_t num_buckets, int32_t bucket_shift,
                                                     int32_t bucket_mask,
                                                     int64_t& total_key_count) {
    std::vector<int64_t> bucket_precount(num_buckets, 0);

    // If the input is empty, just return.
    if (minimizers.size() == 0) {
        return bucket_precount;
    }

    // Count all the keys.
    total_key_count = 0;

    // Preload the first key. This is OK because we return
    // if minimizers.size() == 0.
    minkey_t prev_key = (uint64_t) mindex::Seed::DecodeKey(minimizers[0]);

    for (size_t i = 0; i < minimizers.size(); i++) {
        minkey_t key = (uint64_t) mindex::Seed::DecodeKey(minimizers[i]);
        if (key < prev_key) {
            std::cerr << "prev: " << prev_key << std::endl;
            std::cerr << "key: " << key << std::endl;
            fprintf(stderr, "ERROR: Keys are not sorted in ascending order! Exiting.");
            exit(1);
        }

        if (key != prev_key) {
            int32_t bucket = (prev_key >> bucket_shift) & bucket_mask;
            bucket_precount[bucket] += 1;
            total_key_count += 1;
        }
        prev_key = key;
    }
    {
        // Store the last hash key in the list.
        int32_t bucket = (prev_key >> bucket_shift) & bucket_mask;
        bucket_precount[bucket] += 1;
        total_key_count += 1;
    }

    return bucket_precount;
}

void DenseIndex::ConstructHash_(const std::vector<mindex128_t>& minimizers) {
    // Cleanup.
    for (auto& hash : hashes_) {
        hash.clear();
    }
    hash_ranges_.clear();

    if (minimizers.size() == 0) {
        return;
    }

    // Count the size of each hash.
    int64_t total_key_count = 0;
    auto bucket_precount =
        BucketPrecount_(minimizers, NUM_HASH_BUCKETS, bucket_shift_, bucket_mask_, total_key_count);

    // Reserve the space for hashes.
    for (size_t i = 0; i < hashes_.size(); ++i) {
        #if defined MINIMIZER_INDEX2_USING_SPARSEHASH or defined MINIMIZER_INDEX2_USING_DENSEHASH
                hashes_[i].resize(bucket_precount[i]);
        #elif defined MINIMIZER_INDEX2_USING_UNORDERED_MAP
                hashes_[i].reserve(bucket_precount[i]);
        #endif
    }
    hash_ranges_.reserve(total_key_count);

    // Count all the keys.
    num_keys_ = 0;

    // This hold the beginning of a streak of the same
    // key and the length of the streak.
    SeedSpan new_hash_val(0, 0);

    // Preload the first key. This is OK because we return
    // if minimizers.size() == 0.
    Seed mm = mindex::Seed(minimizers[0]);
    minkey_t prev_key = mm.key;
    // minkey_t prev_key = mindex::Seed::DecodeKey(minimizers[0]);

    // Initialize the hash.
    for (size_t i = 0; i < minimizers.size(); i++) {
        mm = mindex::Seed(minimizers[i]);
        minkey_t key = mm.key;
        // minkey_t key = mindex::Seed::DecodeKey(minimizers[i]);

        if (key == prev_key) {
            new_hash_val.end = i + 1;

        } else {
            int32_t bucket = (prev_key >> bucket_shift_) & bucket_mask_;
            hashes_[bucket][prev_key] = hash_ranges_.size();  // new_hash_val.start;
            hash_ranges_.emplace_back(new_hash_val);
            new_hash_val.start = i;
            new_hash_val.end = i + 1;
            num_keys_ += 1;
        }

        prev_key = key;
    }
    {  // Store the last hash key in the list.
        int32_t bucket = (prev_key >> bucket_shift_) & bucket_mask_;
        hashes_[bucket][prev_key] = hash_ranges_.size();
        hash_ranges_.emplace_back(new_hash_val);
        num_keys_ += 1;
    }
}

size_t DenseIndex::CountKeys_(const std::vector<mindex128_t>& minimizers) const {
    if (minimizers.size() == 0) {
        return 0;
    }

    size_t key_count = 1;   // There is at least one key.

    // Preload the first key. This is OK because we return
    // if minimizers.size() == 0.
    minkey_t prev_key = mindex::Seed::DecodeKey(minimizers[0]);

    for (size_t i = 0; i < seeds_.size(); i++) {
        minkey_t key = mindex::Seed::DecodeKey(minimizers[i]);
        if (key != prev_key) {
            key_count += 1;
        }
        prev_key = key;
    }

    return key_count;
}

void DenseIndex::CalcOccurrenceThreshold_(size_t& occ_max, size_t& occ_cutoff,
                                              size_t& occ_singletons, double& occ_avg) const {
    occ_max = 0;
    occ_cutoff = 0;
    occ_singletons = 0;
    occ_avg = 0;

    if (seeds_.size() == 0) {
        return;
    }

    size_t sum_before = 0;

    // This will hold all key occurrence counts.
    std::vector<int32_t> hit_counts(num_keys_ + 1, 0);

    // Preload the first key. This is OK because we return
    // if minimizers.size() == 0.
    minkey_t prev_key = mindex::Seed::DecodeKey(seeds_[0]);
    size_t curr_key_id = 0;
    size_t prev_start = 0;
    // Initialize the hash.
    for (size_t i = 0; i < seeds_.size(); i++) {
        minkey_t key = mindex::Seed::DecodeKey(seeds_[i]);
        if (key != prev_key) {
            size_t dist = i - prev_start;
            hit_counts[curr_key_id] = dist;
            if (dist == 1) {
                occ_singletons += 1;
            }
            sum_before += dist;
            prev_start = i;
            ++curr_key_id;
        }
        prev_key = key;
    }
    {  // Add the last key.
        size_t dist = seeds_.size() - prev_start;
        hit_counts[curr_key_id] = dist;
        sum_before += dist;
        if (dist == 1) {
            occ_singletons += 1;
        }
        ++curr_key_id;
    }

    // Sort the counts to get the percentil values.
    std::sort(hit_counts.begin(), hit_counts.end());

    // Calculate the cutoff threshold.
    occ_max = hit_counts.back();
    int64_t cutoff_id = (size_t)floor(num_keys_ * (1.0 - params_->freq_percentil));
    occ_cutoff = (size_t)hit_counts[cutoff_id];
    occ_cutoff = std::max((size_t) 1, std::max(occ_cutoff, (size_t)params_->min_occ_cutoff));

    occ_avg = ((double)sum_before) / ((double)hit_counts.size());
}

// void DenseIndex::RemoveFrequentMinimizers_(std::vector<mindex128_t>& minimizers,
//                                                size_t cutoff) const {
//     if (minimizers.size() == 0) {
//         return;
//     }

//     // Preload the first key. This is OK because we return
//     // if minimizers.size() == 0.
//     minkey_t prev_key = mindex::Seed::DecodeKey(minimizers[0]);
//     size_t prev_start = 0;

//     size_t dest_id = 0;
//     for (size_t source_id = 0; source_id < minimizers.size(); source_id++) {
//         minkey_t key = mindex::Seed::DecodeKey(minimizers[source_id]);

//         // If above threshold, rewind.
//         if (key != prev_key) {
//             size_t dist = source_id - prev_start;
//             if (dist > cutoff) {
//                 dest_id -= dist;
//             }
//             prev_start = source_id;
//         }

//         minimizers[dest_id] = minimizers[source_id];
//         ++dest_id;

//         prev_key = key;
//     }

//     {
//         size_t dist = minimizers.size() - prev_start;
//         if (dist > cutoff) {
//             dest_id -= dist;
//         }
//     }

//     minimizers.resize(dest_id);
// }

int DenseIndex::Load(const std::string& index_path) { return 1; }

int DenseIndex::Store(const std::string& index_path) { return 1; }

std::vector<mindex::SeedHitPacked> DenseIndex::CollectHits(const std::string& seq) {
    return CollectHits((const int8_t*)&seq[0], seq.size(), 0);
}

int32_t DenseIndex::CalcBucketShift_() const {
    return std::max(0, (k() * 2 - NUM_HASH_BUCKET_BITS));
}

int32_t DenseIndex::CalcBucketMask_() const {
    return (1 << NUM_HASH_BUCKET_BITS) - 1;
}

bool DenseIndex::SkipHP_(const int8_t* seq, size_t seq_len, int32_t pos, int32_t max_homopolymer_len, int32_t& ret_pos) const {
    ret_pos = pos;
    int8_t b = seq[pos];
    int32_t hp_len = 0;
    if ((pos + 1) < seq_len && seq[pos + 1] == b) {
        for (hp_len = 1; hp_len < max_homopolymer_len && (pos + hp_len) < seq_len; ++hp_len) {
            if (seq[pos + hp_len] != b) {
                break;
            }
        }
        pos += hp_len - 1;
        if (hp_len == max_homopolymer_len) {
            return false;
        }
    }
    ret_pos = pos;
    return true;
}

std::vector<int32_t> DenseIndex::CalcHPSpan_(const int8_t* seq, size_t seq_len, int32_t k, bool homopolymer_suppression, int32_t max_homopolymer_len) const {
    std::vector<int32_t> ret(seq_len, 0);

    std::deque<ind_t> base_pos;

    ind_t last_i = 0;
    // Initialize the first k bases.
    for (ind_t i = 0; i < (seq_len); i++) {
        // int8_t b = seq[i];

        // If enabled, skip same bases.
        if (homopolymer_suppression) {
            int32_t new_i = 0;
            bool hp_len_good = SkipHP_(seq, seq_len, i, max_homopolymer_len, new_i);
            if (hp_len_good == false) {
                continue;
            }
            i = new_i;
        }

        // Add the base to the buffer.
        base_pos.push_back(i);
        last_i = i + 1;
        if (base_pos.size() >= k) {
            break;
        }
    }

    for (ind_t i = (last_i); i < seq_len; i++) {
        // int8_t b = seq[i];

        // If enabled, skip same bases.
        if (homopolymer_suppression) {
            int32_t new_i = 0;
            bool hp_len_good = SkipHP_(seq, seq_len, i, max_homopolymer_len, new_i);
            if (hp_len_good == false) {
                continue;
            }
            i = new_i;
        }

        int32_t first_pos = base_pos.front();
        int32_t last_pos = base_pos.back();

        // Add the base to the buffer.
        base_pos.push_back(i);
        base_pos.pop_front();

        ret[first_pos] = last_pos - first_pos + 1;
    }

    return ret;
}

void DenseIndex::VerboseSeeds_(std::ostream& os) {
    minkey_t prev_key = mindex::Seed::DecodeKey(seeds_[0]);
    size_t curr_key_id = 0;
    size_t prev_start = 0;

    os << "bucket_shift = " << bucket_shift_ << ", bucket_mask = " << bucket_mask_ << std::endl;

    // Initialize the hash.
    for (size_t i = 0; i < seeds_.size(); i++) {
        // minkey_t key = mindex::Seed::DecodeKey(seeds_[i]);

        Seed mm(seeds_[i]);
        int32_t bucket = (mm.key >> bucket_shift_) & bucket_mask_;
        os << "[" << i << "] seq_id = " << mm.seq_id << ", pos = " << mm.pos << ", flag = " << mm.flag << ", key = " << mm.key << ", bucket = " << bucket << ", hashes_[bucket].size() = " << hashes_[bucket].size() << std::endl;

        // if (key != prev_key) {
        //     size_t dist = i - prev_start;
        //     hit_counts[curr_key_id] = dist;
        //     if (dist == 1) {
        //         occ_singletons += 1;
        //     }
        //     sum_before += dist;
        //     prev_start = i;
        //     ++curr_key_id;
        // }
        // prev_key = key;
    }
    os << "hashes_.size() = " << hashes_.size() << std::endl;
    // for (size_t i = 0; i < hashes_.size(); i++) {
    //     os << "  [" << i << "] hashes_[i].size() = " << hashes_[i].size() << std::endl;
    // }

}



// #define EXPERIMENTAL_QUERY_MASK
std::vector<mindex::SeedHitPacked> DenseIndex::CollectHits(const int8_t* seq,
                                                                    size_t seq_len,
                                                                    indid_t seq_id) {
    std::vector<mindex::SeedHitPacked> hits;

    if (seq == NULL || seq_len == 0) {
        return hits;
    }

    // TicToc tt1;
    // tt1.start();

    // Queries are usually very small so we can
    // afford to preallocate a lot of space.
    hits.reserve(seq_len);

    // For lookup, don't index the reverse complement;
    std::vector<mindex128_t> seeds;

    // We must enforce w = 1 in the query sequence because these are not minimizers.
    int ret_val = GenerateSeeds_(
        seeds, seq, seq_len, 0, seq_id, params_->k, 1, !params_->index_only_fwd_strand,
        params_->homopolymer_suppression, params_->max_homopolymer_len, 0, -1);

    // tt1.stop();
    // std::cerr << "Time for GenerateMinimizers: " << tt1.get_millisecs() << std::endl;

#ifdef EXPERIMENTAL_QUERY_MASK
    kx::radix_sort(seeds.begin(), seeds.end());
#endif

    if (ret_val) {
        return hits;
    }

    int32_t k_len = k();

    // TicToc tt2;
    // tt2.start();

    // To properly handle the reverse complement hit coordinates, the reverse minimizer
    // position _end_ location needs to be known. One way would be to store that with
    // the minimizer when generating seeds.
    // Instead, we calculate the end position on the fly, so that the number of bits
    // in a minimizes's position is not reduced.
    std::vector<int32_t> hp_span = (params_->homopolymer_suppression == false)
                                    ? std::vector<int32_t>(seq_len, k_len)
                                    : CalcHPSpan_(seq, seq_len, k_len, params_->homopolymer_suppression, params_->max_homopolymer_len);

    // tt2.stop();
    // std::cerr << "Time for CalcHPSpan_: " << tt2.get_millisecs() << std::endl;

    // TicToc tt3;
    // tt3.start();

    // VerboseSeeds_(std::cerr);

    for (size_t min_id = 0; min_id < seeds.size(); ++min_id) {
        // Decode the minimizer.
        mindex::Seed minimizer(seeds[min_id]);
        bool minimizer_is_rev = minimizer.is_rev();

        // Query mask is an additional filter which will be encoded in the seed hit.
        // It contains fields to mark whether the hit is a tandem repeat of a kmer,
        // or other interesting info which should impact the sorting order.
        int32_t query_mask = 0x0;
        #ifdef EXPERIMENTAL_QUERY_MASK
                if (min_id > 0 && (minimizer.key >> 8) == (mindex::Seed::DecodeKey(seeds[min_id - 1]) >> 8)) {
                    query_mask |= MINIMIZER_HIT_TANDEM_FLAG;
                }
                if ((min_id + 1) < seeds.size() && (minimizer.key >> 8) == (mindex::Seed::DecodeKey(minimizers[min_id + 1]) >> 8)) {
                    query_mask |= MINIMIZER_HIT_TANDEM_FLAG;
                }
        #endif

        // Find and add all the key hits.
        bool is_found = false;
        auto it = FindKey_(minimizer.key, is_found);  // Determines the correct bucket and looks up the key.

        if (is_found) {
            auto& seed_range = hash_ranges_[it->second];
            // Skip very frequent seeds.
            if ((seed_range.end - seed_range.start) <= occ_cutoff_) {
                for (ind_t seed_id = seed_range.start; seed_id < seed_range.end; seed_id++) {
                    mindex::Seed hit(seeds_[seed_id]);
                    bool indexed_minimizer_is_rev = hit.is_rev();
                    ind_t pos = minimizer.pos;
                    bool hit_is_rev = false;

                    if (minimizer_is_rev != indexed_minimizer_is_rev) {
                        hit_is_rev = true;
                        const mindex::SequencePtr& seq = seqs_->GetSeqByID(hit.seq_id);
                        if (seq == nullptr) {
                            std::cerr << "Warning: hit.seq_id does not match a sequence in index! Skipping hit.\n";
                            continue;
                        }
                        hit.pos = seq->len() - (hit.pos + 1);
                        pos += hp_span[pos] - 1;
                    }

                    hits.emplace_back(mindex::SeedHitPacked(hit.seq_id, hit_is_rev, hit.pos, query_mask, pos));
                }
            }
        }
    }

    // tt3.stop();
    // std::cerr << "Time for processing hits: " << tt3.get_millisecs() << std::endl;

    return hits;
}

int DenseIndex::GenerateSeeds_(std::vector<mindex128_t>& minimizers,
                                        const int8_t* seq, ind_t seq_len, ind_t seq_offset, // The seq_offset is the distance from the beginning of the sequence, used for the seed position.
                                        indid_t seq_id, int32_t k, int32_t w, bool use_rc,
                                        bool homopolymer_suppression,
                                        int32_t max_homopolymer_len,
                                        ind_t seq_start, ind_t seq_end) {
    // Sanity check that the seq is not NULL;
    if (!seq) {
        return 1;
    }
    // Sanity check for the size of k. It can't
    // be too large, or it won't fit the uint64_t buffer.
    if (k <= 0 || k >= 31) {
        return 2;
    }
    // Not technically an error if the seq_len is 0,
    // but there's nothing to do, so return.
    if (seq_len == 0) {
        return 0;
    }

    seq_start = (seq_start < 0) ? 0 : seq_start;
    seq_end = (seq_end <= 0 || seq_end > seq_len) ? seq_len : seq_end;

    if (seq_start >= seq_len) {
        return 3;
    }
    if (seq_start >= seq_end) {
        return 4;
    }

    const minkey_t mask = (((uint64_t)1) << (2 * k)) - 1;   // Mask the number of required bits for the k-mer.
    minkey_t buffer = 0x0;     // Holds the current 2-bit seed representation.
    minkey_t buffer_rc = 0x0;  // Holds the reverse complement 2-bit seed at the same position.
    int32_t num_bases_in = 0;   // Number of bases added to the buffer.
    ind_t kmer_span = 0;
    std::deque<ind_t> hp_events;

    for (ind_t pos = seq_start; pos < seq_end; ++pos) {
        int8_t b = seq[pos];
        mindex::Seed new_seed; // (UINT64_T_MAX, INT32_T_MAX, INT32_T_MAX, INT8_T_MAX);
        bool new_seed_set = false;

        if (is_nuc[b]) {
            // If enabled, skip same bases.
            if (homopolymer_suppression) {
                // In this case, find the stretch of homopolymers (at least 1 base), and
                // add that stretch to "seed_span" (the seed_span is iteratively increased every
                // loop iteration).
                int32_t hp_len = 1;
                for (hp_len = 1; hp_len < max_homopolymer_len && (pos + hp_len) < seq_len; ++hp_len) {
                    if (seq[pos + hp_len] != b) {
                        break;
                    }
                }
                hp_events.push_back(hp_len);
                pos += hp_len - 1;
                kmer_span += hp_len;
                // We need to keep track of the length of each HP event. They contribute as a single
                // base, but in reality, the span is larger. This deque keeps track of all hp events
                // in the current k-mer.
                while (hp_events.size() > k) {
                    kmer_span -= hp_events.front();
                    hp_events.pop_front();
                }
                if (hp_len >= max_homopolymer_len) {
                    continue;
                }
            } else {
                kmer_span = std::min(num_bases_in + 1, k);
            }

            // Add the base to the buffer.
            buffer = ((buffer << 2) | ((((uint64_t)nuc_to_twobit[b]))));
            buffer_rc = (buffer_rc >> 2) | ((((uint64_t)nuc_to_twobit_complement[b])) << (k * 2 - 2));
            // Calculate the seed key.
            minkey_t key = buffer & mask;
            minkey_t rev_key = buffer_rc & mask;
            // Skip symmetric k-mers.
            if (buffer == buffer_rc) {
                continue;
            }

            // Determine the orientation of the key.
            int8_t flag = MINIMIZER_FLAG_DEFAULT_FWD;
            if (use_rc && rev_key < key) {
                std::swap(key, rev_key);
                flag = MINIMIZER_FLAG_IS_REV;
            }

            ++num_bases_in;

            if (num_bases_in >= k) {    // Minimap2 has another condition here: "kmer_span < 256". That's because it encodes the kmer span into the seed definition as 8 bits.
                int32_t kmer_start = (pos + 1) - kmer_span; // The 'pos' is the current position which is inclusive. We need to add a +1 to make it non-inclusive, so that the start position is calculated properly.
                new_seed = mindex::Seed(key, seq_id, kmer_start + seq_offset, flag);
                new_seed_set = true;
            }
        } else {
            // If we encountered a non-nucleotide base, this is enough to reset the minimizer coding.
            // First, write out the current minimizer so it doesn't get lost in the void.
            // We only need to write one and not loop through all of them, because when we encounter a seed
            // with the same key then we write out the previous minimizer out. This means that we only need
            // to write out the current minimizer.
            num_bases_in = 0;
            buffer = 0;
            buffer_rc = 0;
            hp_events.clear();
        }

        // Select kmers with a step 'w'.
        if (new_seed_set && (new_seed.pos % w) == 0) {
            minimizers.emplace_back(new_seed.to_128t());
        }
    }

    return 0;
}

const int8_t* DenseIndex::FetchRawSeq(size_t seq_id) const {
    const mindex::SequencePtr& seq = seqs_->GetSeqByID(seq_id);
    if (seq == nullptr) {
        return NULL;
    }

    return &seq->data()[0];
}

std::string DenseIndex::FetchSeqAsString(size_t seq_id, size_t start, size_t end,
                                             bool rev_cmp) const {
    const mindex::SequencePtr& seq = seqs_->GetSeqByID(seq_id);
    if (seq == nullptr) {
        return std::string("");
    }

    if (end <= start) {
        return std::string("");
    }

    if (start >= seq->len()) {
        return std::string("");
    }

    if (end > seq->len()) {
        return std::string("");
    }

    std::string ret;
    if (rev_cmp == false) {
        ret = std::string((const char*)(&seq->data()[start]), (end - start));
    } else {
        size_t seq_len = end - start;
        ret = std::string(seq_len, ' ');
        int64_t end_pos = end - 1;
        for (int64_t i = 0; i < seq_len; i++) {
            int64_t rev_pos = end_pos - i;
            ret[i] = (char)nuc_to_complement[seq->data()[rev_pos]];
        }
    }

    return ret;
}

std::string DenseIndex::FetchSeqAsString(size_t seq_id, bool rev_cmp) const {
    const mindex::SequencePtr& seq = seqs_->GetSeqByID(seq_id);
    if (seq == nullptr) {
        return std::string("");
    }

    return FetchSeqAsString(seq_id, 0, seq->len(), rev_cmp);
}

}  // namespace mindex
