# Introduction
There were several different approaches to calculating the minimizers that I've tried throughout this project.
It would be a shame to get them lost through time, so this document will serve as a reminder.

These are:
1. GenerateMinimizersWithQueue_ (deprecated on Sep 05, 2019)
2. GenerateMinimizersWithRMQ_ (work started on Dec 16, 2017, commit 290de8e421f62a84fc1c3c12245c2d4378be9a4f)
3. Minheap-Minqueue (src/minimizer_index2/minimizer_generator2_minqueue.cc) (somewhere around this time: Dec 09, 2017, commit d16776b6092c3dd59e1b0e47ebd0fbb29bd35490)

It looks like the last commit with all of these methods is this one:
commit 82f073b7ac721d776a2e5e49db2d647f00599420 (HEAD)
Author: Ivan Sovic <ivan.sovic@gmail.com>
Date:   Sat Dec 16 01:17:24 2017 +0100

    Small cleanup in the output formatter.
    
## 1. Queue-based method
The best one until recently was the queue-based approach which wraps the minimizer generation in a generator class. The generator takes in a new seed, and outputs the minimizers for the current window right away.
The method was deprecated on Sep 04, 2019 in favor of a slightly faster approach.
The original implementation is given here.

- Header for the minimizer generator:
```
#include <index/minimizer_index_util.h>

#include <algorithm>
#include <sstream>

namespace mindex {

std::string MinimizerKeyToString(minkey_t seed, int32_t k) {
    std::stringstream ss;
    for (int32_t i = 0; i < k; i++) {
        ss << (char)twobit_to_nuc[seed & 0x03];
        seed = seed >> 2;
    }
    std::string seed_str;
    seed_str = ss.str();
    std::reverse(seed_str.begin(), seed_str.end());
    return seed_str;
}

}  // namespace mindex

```

- Code for the minimizer generator:
```
#include <index/minimizer_generator.h>

namespace mindex {

std::unique_ptr<MinimizerGenerator> createMinimizerGenerator(int32_t window_len) {
    return std::unique_ptr<MinimizerGenerator>(new MinimizerGenerator(window_len));
}

MinimizerGenerator::MinimizerGenerator(int32_t window_len)
    : window_len_(window_len), num_added_keys_(0) {}

MinimizerGenerator::~MinimizerGenerator() {}

int MinimizerGenerator::Yield(const mindex::Minimizer& seed_in,
                              std::vector<mindex::Minimizer>& seeds_out) {
    if (num_added_keys_ < window_len_) {
        // Remove smaller elements if any.
        while ((!q_.empty()) && seed_in.key <= std::get<0>(q_.back()).key) {
            q_.pop_back();
        }

    } else {
        // Return one or more elements of the window
        // with the minimal key value.
        seeds_out.clear();
        if (q_.size() > 0) {
            minkey_t key = std::get<0>(q_.front()).key;
            auto it_next = q_.begin();
            while (it_next != q_.end() && std::get<0>(*(it_next)).key == key) {
                seeds_out.emplace_back(std::get<0>(*(it_next)));
                ++it_next;
            }
        }

        // Remove the elements which are out of this window.
        while ((!q_.empty()) && std::get<1>(q_.front()) <= (num_added_keys_ - window_len_)) {
            q_.pop_front();
        }
        // Remove smaller elements if any.
        while ((!q_.empty()) && seed_in.key <= std::get<0>(q_.back()).key) {
            q_.pop_back();
        }
    }

    q_.push_back(std::make_pair(seed_in, num_added_keys_));
    num_added_keys_ += 1;

    if (num_added_keys_ <= window_len_) {
        return 2;
    }

    return 0;
}

}  // namespace mindex

```

- Headers for the high-level functions to generate the minimizers:
```

    /*
     * For a given sequence, generates and returns a set of minimizers.
     * The `minimizers` vector is never cleared but only inserted to, so minimizers
     * can be accumulated through multiple calls to this method.
     * The GenerateMinimizersGeneric_ method also splits the `seq` at non-[ACTG] bases
     * and generates minimizers only for the sequences which do not contain them.
     * It calls either GenerateMinimizersWithRMQ_ or GenerateMinimizersWithQueue_ to
     * actually generate the minimizers on those clean subsequences.
     */
    static int GenerateMinimizersGeneric_(std::vector<mindex128_t>& minimizers, const int8_t* seq,
                                          ind_t seq_len, indid_t seq_id, int32_t k, int32_t w,
                                          bool use_rc, bool homopolymer_suppression,
                                          int32_t max_homopolymer_run, ind_t seq_start,
                                          ind_t seq_end);

    /*
     * Expects that seq contains only [ACTG] bases.
     * Uses a deque to add/remove seeds and yield minimizers.
     */
    static int GenerateMinimizersWithQueue_(std::vector<mindex128_t>& minimizers, const int8_t* seq,
                                            ind_t seq_len, ind_t seq_offset, indid_t seq_id,
                                            int32_t k, int32_t w, bool use_rc,
                                            bool homopolymer_suppression,
                                            int32_t max_homopolymer_run);

```

- Code for the high-level functions to generate the minimizers:
```
int MinimizerIndex::GenerateMinimizersWithQueue_(std::vector<mindex128_t>& minimizers,
                                                 const int8_t* seq, ind_t seq_len, ind_t seq_offset,
                                                 indid_t seq_id, int32_t k, int32_t w, bool use_rc,
                                                 bool homopolymer_suppression,
                                                 int32_t max_homopolymer_len) {
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

    // Preallocate local memory for the minimizers.
    size_t num_minimizers = 0;
    std::vector<mindex128_t> local_minimizers;
    local_minimizers.reserve(seq_len);

    // This will keep track of the keys which were already used and enable minimizer compression.
    std::vector<bool> is_seed_used(seq_len, false);

    const minkey_t mask = (((uint64_t)1) << (2 * k)) - 1;
    minkey_t buffer = 0x0;     // Holds the current 2-bit seed representation.
    minkey_t buffer_rc = 0x0;  // Holds the reverse complement 2-bit seed at the same position.
    int32_t shift = 0;         // We keep track of the amount to shift after each 'N' base.

    auto mg = createMinimizerGenerator(w);

    std::deque<ind_t> kmer_starts;

    ind_t last_i = 0;
    // Initialize the first k bases.
    for (ind_t i = 0; i < seq_len; i++) {
        int8_t b = seq[i];

        // If enabled, skip same bases.
        if (homopolymer_suppression) {
            int32_t hp_len = 0;
            if ((i + 1) < seq_len && seq[i + 1] == b) {
                for (hp_len = 1; hp_len < max_homopolymer_len && (i + hp_len) < seq_len; ++hp_len) {
                    if (seq[i + hp_len] != b) {
                        break;
                    }
                }
                i += hp_len - 1;
                if (hp_len == max_homopolymer_len) {
                    continue;
                }
            }
        }

        // Add the base to the buffer.
        buffer = (buffer << 2) | ((((uint64_t)nuc_to_twobit[b])));
        buffer_rc = (buffer_rc >> 2) | ((((uint64_t)nuc_to_twobit_complement[b])) << (k * 2 - 2));
        kmer_starts.push_back(i);
        last_i = i + 1;
        ++shift;
        if (shift >= (k - 1)) {
            break;
        }
    }

    // Place to store all the minimizers of a window.
    std::vector<mindex::Minimizer> seeds_out;
    mindex::Minimizer seed_out(0, 0, 0, false);

    // std::cerr << "last_i = " << last_i << ", seq_len = " << seq_len << std::endl;

    for (ind_t i = (last_i); i < seq_len; i++) {
        int8_t b = seq[i];

        // If enabled, skip same bases.
        if (homopolymer_suppression) {
            int32_t hp_len = 0;
            if ((i + 1) < seq_len && seq[i + 1] == b) {
                for (hp_len = 1; hp_len < max_homopolymer_len && (i + hp_len) < seq_len; ++hp_len) {
                    if (seq[i + hp_len] != b) {
                        break;
                    }
                }
                i += hp_len - 1;
                if (hp_len == max_homopolymer_len) {
                    continue;
                }
            }
        }

        // Add the base to the buffer.
        buffer = (buffer << 2) | ((((uint64_t)nuc_to_twobit[b])));
        buffer_rc = (buffer_rc >> 2) | ((((uint64_t)nuc_to_twobit_complement[b])) << (k * 2 - 2));
        kmer_starts.push_back(i);

        // Calculate the seed key.
        minkey_t key = buffer & mask;
        minkey_t rev_key = buffer_rc & mask;

        // key = InvertibleHash_(key);
        // rev_key = InvertibleHash_(rev_key);

        int8_t flag = MINIMIZER_FLAG_DEFAULT_FWD;

        if (use_rc && rev_key < key) {
            key = rev_key;
            flag = MINIMIZER_FLAG_IS_REV;
        }

        // Get the minimizers.
        ind_t kmer_start = kmer_starts.front();
        kmer_starts.pop_front();
        int ret_val_mg = mg->Yield(mindex::Minimizer(key, seq_id, kmer_start, flag), seeds_out);

        // If something went wrong, skip adding the minimizers.
        if (ret_val_mg) {
            // std::cerr << "Warning: i = " << i << ", mg->Yield returned with value: " << ret_val_mg << std::endl;
            continue;
        }

        // Add each minimizer of the window, but skip duplicates.
        for (auto& seed_out : seeds_out) {
            if (is_seed_used[seed_out.pos]) {
                continue;
            }

            is_seed_used[seed_out.pos] = true;

            seed_out.pos += seq_offset;
            local_minimizers.emplace_back(seed_out.to_128t());
            ++num_minimizers;
        }
    }

    // Get the last key out. Separated outside the loop to save from
    // another branching. The values to Yield are dummy in this cass,
    // because we only need to get the seed_out
    {
        // Get the minimizers.
        int ret_val_mg = mg->Yield(mindex::Minimizer(0, 0, 0, 0), seeds_out);

        if (!ret_val_mg) {
            // Add each minimizer of the window, but skip duplicates.
            for (auto& seed_out : seeds_out) {
                if (is_seed_used[seed_out.pos]) {
                    continue;
                }

                is_seed_used[seed_out.pos] = true;

                seed_out.pos += seq_offset;
                local_minimizers.emplace_back(seed_out.to_128t());
                ++num_minimizers;
            }
        }
    }

    // Manually reserve so that STL doesn't double.
    if ((minimizers.size() + num_minimizers) > minimizers.capacity()) {
        minimizers.reserve((minimizers.size() + num_minimizers) * 1.1);
    }

    // Insert the local copy of the minimizers.
    minimizers.insert(minimizers.end(), local_minimizers.begin(),
                      local_minimizers.begin() + num_minimizers);

    return 0;
}

int MinimizerIndex::GenerateMinimizersGeneric_(std::vector<mindex128_t>& minimizers,
                                               const int8_t* seq, ind_t seq_len,
                                               indid_t seq_id, int32_t k, int32_t w, bool use_rc,
                                               bool homopolymer_suppression,
                                               int32_t max_homopolymer_run,
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

    // The seq will be split in parts separated by non-[ACTG] bases.
    // Parts shorter than seed_len will be skipped.
    std::vector<ind_t> split_start;
    std::vector<ind_t> split_len;
    ind_t start = seq_start;
    for (ind_t i = seq_start; i < seq_end; ++i) {
        if (!is_nuc[seq[i]]) {
            ind_t curr_len = i - start;
            if (curr_len >= k) {
                split_start.emplace_back(start);
                split_len.emplace_back(curr_len);
            }
            start = i + 1;
        }
    }
    if (start < seq_end and (seq_end - start) >= k) {
        split_start.emplace_back(start);
        split_len.emplace_back(seq_end - start);
    }

    // Process each chunk separately, and give it the start offset.
    for (size_t split_id = 0; split_id < split_start.size(); ++split_id) {
        GenerateMinimizersWithQueue_(minimizers, seq + split_start[split_id], split_len[split_id],
                                     split_start[split_id], seq_id, k, w, use_rc,
                                     homopolymer_suppression, max_homopolymer_run);
    }

    return 0;
}
```