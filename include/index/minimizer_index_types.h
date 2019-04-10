/*
 * minimizer_index_types.h
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MINIMIZER_INDEX_2_TYPES_H_
#define SRC_MINIMIZER_INDEX_2_TYPES_H_

#include <cstdint>

typedef unsigned __int128 mindex128_t;

// Minimizer index int. 32 bit for now, but in case we decide to switch to 64 bits.
// We want to keep it low by default to reduce memory usage, and for larger seqs,
// we'll just block.
typedef int32_t ind_t;

// Type used for the ID of a particular sequence in the index.
typedef int32_t indid_t;

// Minimizer key type.
typedef uint64_t minkey_t;

static const mindex128_t MINIMIZER_64bit_MASK = (((mindex128_t)0x0FFFFFFFFFFFFFFFF));
static const mindex128_t MINIMIZER_32bit_MASK =
    (((mindex128_t)0x0000000007FFFFFFF));  // Mask the sign bit.
static const mindex128_t MINIMIZER_32bit_MASK_FULL = (((mindex128_t)0x000000000FFFFFFFF));

#endif
