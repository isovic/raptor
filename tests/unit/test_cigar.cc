#include <gtest/gtest.h>

#include <aligner/cigar.h>

TEST(Cigar, ScoreCigarAlignment1) {
    std::vector<raptor::CigarOp> cigar = {};
    int32_t match = 1, mismatch = 1, gap_open = 1, gap_ext = 1;
    int64_t result = raptor::ScoreCigarAlignment(cigar, match, mismatch, gap_open, gap_ext);
    int64_t expected = 0;
    ASSERT_EQ(result, expected);
}

TEST(Cigar, ScoreCigarAlignment2) {
    std::vector<raptor::CigarOp> cigar = { {'=', 10} };
    int32_t match = 1, mismatch = 1, gap_open = 1, gap_ext = 1;
    int64_t result = raptor::ScoreCigarAlignment(cigar, match, mismatch, gap_open, gap_ext);
    int64_t expected = 10;
    ASSERT_EQ(result, expected);
}

TEST(Cigar, ScoreCigarAlignment3) {
    std::vector<raptor::CigarOp> cigar = { {'=', 10}, {'I', 1}, {'X', 2}, {'D', 1}, {'=', 5} };
    int32_t match = 2, mismatch = -4, gap_open = -2, gap_ext = -4;
    int64_t result = raptor::ScoreCigarAlignment(cigar, match, mismatch, gap_open, gap_ext);
    int64_t expected = 20 - 2 - 8 - 2 + 10;
    ASSERT_EQ(result, expected);
}

TEST(Cigar, ScoreCigarAlignment4) {
    std::vector<raptor::CigarOp> cigar = { {'=', 10}, {'I', 5}, {'D', 7}, {'X', 2}, {'=', 5} };
    int32_t match = 2, mismatch = -4, gap_open = -2, gap_ext = -4;
    int64_t result = raptor::ScoreCigarAlignment(cigar, match, mismatch, gap_open, gap_ext);
    int64_t expected = 20 + (-2 - 4 * 4) + (-2 - 6 * 4) - 8 + 10;
    ASSERT_EQ(result, expected);
}

TEST(Cigar, ScoreCigarAlignment5) {
    // Only extended CIGAR operations are defined: '=', 'X', 'I', 'D'.
    // This tests occurrences of other characters, which will simply be ignored.
    std::vector<raptor::CigarOp> cigar = { {'M', 10}, {'I', 5}, {'D', 7}, {'X', 2}, {'N', 5} };
    int32_t match = 2, mismatch = -4, gap_open = -2, gap_ext = -4;
    int64_t result = raptor::ScoreCigarAlignment(cigar, match, mismatch, gap_open, gap_ext);
    int64_t expected = 0 + (-2 - 4 * 4) + (-2 - 6 * 4) - 8 + 0;
    ASSERT_EQ(result, expected);
}
