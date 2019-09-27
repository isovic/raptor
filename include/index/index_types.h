/*
 * index_types.h
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX_2_TYPES_H_
#define SRC_MINIMIZER_INDEX_2_TYPES_H_

#include <cstdint>

typedef unsigned __int128 mindex128_t;

// Seed index int. 32 bit for now, but in case we decide to switch to 64 bits.
// We want to keep it low by default to reduce memory usage, and for larger seqs,
// we'll just block.
typedef int32_t ind_t;

// Type used for the ID of a particular sequence in the index.
typedef int32_t indid_t;

// Seed key type.
typedef uint64_t minkey_t;

static const mindex128_t MINIMIZER_64bit_MASK = (((mindex128_t)0x0FFFFFFFFFFFFFFFF));
static const mindex128_t MINIMIZER_32bit_MASK =
    (((mindex128_t)0x0000000007FFFFFFF));  // Mask the sign bit.
static const mindex128_t MINIMIZER_32bit_MASK_FULL = (((mindex128_t)0x000000000FFFFFFFF));

struct SeedSpan {
    SeedSpan() : start(0), end(0) {}
    SeedSpan(size_t _start, ind_t _end) : start(_start), end(_end) {}
    size_t start = 0;
    ind_t end = 0;
};

/*
 * Experimental code test performance with different hash maps.
 */
#define MINIMIZER_INDEX2_USING_DENSEHASH
// #define MINIMIZER_INDEX2_USING_SPARSEHASH
// #define MINIMIZER_INDEX2_USING_UNORDERED_MAP

#ifdef MINIMIZER_INDEX2_USING_UNORDERED_MAP
#include <unordered_map>
typedef std::unordered_map<minkey_t, int64_t, std::hash<minkey_t> > SeedHashType2;
#endif

#ifdef MINIMIZER_INDEX2_USING_DENSEHASH
#include <sparsehash/dense_hash_map>
using google::dense_hash_map;  // namespace where class lives by default
// typedef dense_hash_map<minkey_t, int64_t, std::hash<minkey_t> > SeedHashType2;
using SeedHashType2 = dense_hash_map<minkey_t, int64_t, std::hash<minkey_t>>;
#endif

#ifdef MINIMIZER_INDEX2_USING_SPARSEHASH
#include <sparsehash/sparse_hash_map>
using google::sparse_hash_map;  // namespace where class lives by default
// typedef sparse_hash_map<minkey_t, int64_t, std::hash<minkey_t> > SeedHashType2;
using SeedHashType2 = sparse_hash_map<minkey_t, int64_t, std::hash<minkey_t>>;
#endif

const minkey_t MINIMIZER_INDEX_EMPTY_HASH_KEY = (minkey_t)0xFFFFFFFFFFFFFFFF;

namespace mindex {

enum class IndexType {
    Minimizer,
    Dense,
    Undefined
};

}

#endif
