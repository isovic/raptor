#include <gtest/gtest.h>

#include <index/index_params.h>
#include <index/dense_index.h>
#include <index/seed.hpp>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <log/log_system.h>

#include <tests/unit/test_index_common.h>

TEST(DenseIndexTest, AddSequenceTest1) {
    auto index_params = mindex::createIndexParams();
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {"AAAAA", "CCCCC", "TTTTTG", "GGGGG"};
    std::vector<std::string> headers = {"header1", "header2", "header3", "header4"};

    // Add them to the index. This only adds the data, but does not build the index.
    for (size_t i = 0; i < seqs.size(); i++) {
        index->AddSequence(seqs[i], headers[i]);
    }

    ASSERT_EQ(index->params(), index_params);

    // Check that all sequences were added.
    ASSERT_EQ((size_t) index->num_seqs(), seqs.size());

    // Attempt to access all the info for each sequence.
    // FetchSeqAsString also indirectly tests the data_ vector.
    int64_t expected_total_len = 0;
    for (size_t i = 0; i < index->num_seqs(); i++) {
        ASSERT_EQ(index->header(i), headers[i]);
        ASSERT_EQ(index->len(i), seqs[i].size());
        ASSERT_EQ(index->FetchSeqAsString(i, false), seqs[i]);
        expected_total_len += seqs[i].size();
    }

    ASSERT_EQ((size_t) index->total_len(), expected_total_len);
}

TEST(DenseIndexTest, AddSequenceTest2) {
    auto index_params = mindex::createIndexParams();
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {"AAAAA", "CCCCC", "TTTTTG", "GGGGG"};
    std::vector<std::string> headers = {"header1", "header2", "header3", "header4"};

    // Add them to the index. This only adds the data, but does not build the index.
    for (size_t i = 0; i < seqs.size(); i++) {
        index->AddSequence((const int8_t *) seqs[i].c_str(), seqs[i].size(), (const char *) headers[i].c_str(), headers[i].size());
    }

    ASSERT_EQ(index->params(), index_params);

    // Check that all sequences were added.
    ASSERT_EQ((size_t) index->num_seqs(), seqs.size());

    // Attempt to access all the info for each sequence.
    // FetchSeqAsString also indirectly tests the data_ vector.
    int64_t expected_total_len = 0;
    for (size_t i = 0; i < index->num_seqs(); i++) {
        ASSERT_EQ(index->header(i), headers[i]);
        ASSERT_EQ(index->len(i), seqs[i].size());
        ASSERT_EQ(index->FetchSeqAsString(i, false), seqs[i]);
        expected_total_len += seqs[i].size();
    }

    ASSERT_EQ((size_t) index->total_len(), expected_total_len);
}

TEST(DenseIndexTest, AddSequencesTest3) {
    auto index_params = mindex::createIndexParams();
    auto index1 = mindex::createDenseIndex(index_params);
    auto index2 = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {"AAAAA", "CCCCC", "TTTTTG", "GGGGG"};
    std::vector<std::string> headers = {"header1", "header2", "header3", "header4"};

    // Add the sequences to the index in a way that was already extensively tested above.
    for (size_t i = 0; i < seqs.size(); i++) {
        index1->AddSequence(seqs[i], headers[i]);
    }

    // Add the sequences to the other index using the wrapper method.
    index2->AddSequences(seqs, headers);

    // Compare that the end result is the same.
    ASSERT_EQ(index2->num_seqs(), index1->num_seqs());

    for (size_t i = 0; i < index1->num_seqs(); ++i) {
        const auto& seq1 = index1->seqs()->GetSeqByID(i);
        const auto& seq2 = index2->seqs()->GetSeqByID(i);
        ASSERT_EQ(seq2->data(), seq1->data());
    }
}

TEST(DenseIndexTest, FetchSeqAsStringTest1) {
    auto index_params = mindex::createIndexParams();
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {"AAAAAAAAAACCCCCTTTTTGGGGGGGGGG", "ACTGGTCA"};
    std::vector<std::string> headers = {"header1", "header2"};

    // Add them to the index. This only adds the data, but does not build the index.
    for (size_t i = 0; i < seqs.size(); i++) {
        index->AddSequence(seqs[i], headers[i]);
    }

    // Test fetching the full sequence.
    ASSERT_EQ(index->FetchSeqAsString(0, false), seqs[0]);

    {
        // Test reverse complement.
        std::string expected("CCCCCCCCCCAAAAAGGGGGTTTTTTTTTT");
        ASSERT_EQ(index->FetchSeqAsString(0, true), expected);
    }

    {
        // Test fetching between an arbitrary start and end position.
        size_t seq_id = 0, start = 0, end = 10;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), seqs[seq_id].substr(start, (end - start)));
    }

    {
        // Test fetching between an arbitrary start and end position.
        size_t seq_id = 0, start = 0, end = 11;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), seqs[seq_id].substr(start, (end - start)));
    }

    {
        // Test fetching between an arbitrary start and end position of the second sequence.
        size_t seq_id = 1, start = 2, end = 4;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), seqs[seq_id].substr(start, (end - start)));
    }

    {
        // If start is smaller than the end, the method should return an empty string.
        size_t seq_id = 0, start = 11, end = 0;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), std::string());
    }

    {
        // Out of range seq_id should return an empty string.
        size_t seq_id = 10, start = 0, end = 11;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), std::string());
    }

    {
        // Out of range seq_id should return an empty string.
        size_t seq_id = -1, start = 0, end = 11;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), std::string());
    }

    {
        // Out of range start and end should return an empty string.
        size_t seq_id = 10, start = 0, end = 1000;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), std::string());
    }

    {
        // Out of range start and end should return an empty string.
        size_t seq_id = 10, start = 1000, end = 10000;
        ASSERT_EQ(index->FetchSeqAsString(seq_id, start, end, false), std::string());
    }
}

TEST(DenseIndexTest, FetchRawSeqTest1) {
    auto index_params = mindex::createIndexParams();
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {"AAAAA", "CCCCC", "TTTTTG", "GGGGG"};
    std::vector<std::string> headers = {"header1", "header2", "header3", "header4"};

    // Add them to the index. This only adds the data, but does not build the index.
    for (size_t i = 0; i < seqs.size(); i++) {
        index->AddSequence(seqs[i], headers[i]);
    }

    // Check if all raw data can be fetched properly.
    for (size_t i = 0; i < index->num_seqs(); i++) {
        ASSERT_EQ(strncmp((const char *) index->FetchRawSeq(i), seqs[i].c_str(), seqs[i].size()), 0);
    }
}

TEST(DenseIndexTest, TestBuild1) {
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AAAAAAAAAA",
                                    };
    std::vector<std::string> headers = {
                                    "fake_seq1",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(0, 0, 0, false),
                                        mindex::Seed::Encode(0, 0, 1, false),
                                        mindex::Seed::Encode(0, 0, 2, false),
                                        mindex::Seed::Encode(0, 0, 3, false),
                                        mindex::Seed::Encode(0, 0, 4, false),
                                        mindex::Seed::Encode(0, 0, 5, false),
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild2) {
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    // "AGCTTTTCAT",
                                    "CAAAAAAAAA",
                                    // "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    // "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTA"
                                    };
    std::vector<std::string> headers = {
                                    "fake_seq1",
                                    // "ecoli:0-70",
                                    // "lambda:0-80"
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{   mindex::Seed::Encode(0, 0, 1, false),
                                        mindex::Seed::Encode(0, 0, 2, false),
                                        mindex::Seed::Encode(0, 0, 3, false),
                                        mindex::Seed::Encode(0, 0, 4, false),
                                        mindex::Seed::Encode(0, 0, 5, false),
                                        mindex::Seed::Encode(256, 0, 0, false)
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild3) {
    /*
     * Tests building on a more complex sequence.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCA",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli:0-20",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(2, 0, 2, true),
                                        mindex::Seed::Encode(9, 0, 1, true),
                                        mindex::Seed::Encode(39, 0, 0, true),
                                        mindex::Seed::Encode(56, 0, 6, true),
                                        mindex::Seed::Encode(121, 0, 14, false),
                                        mindex::Seed::Encode(131, 0, 8, true),
                                        mindex::Seed::Encode(180, 0, 12, true),
                                        mindex::Seed::Encode(224, 0, 5, true),
                                        mindex::Seed::Encode(288, 0, 9, true),
                                        mindex::Seed::Encode(301, 0, 13, true),
                                        mindex::Seed::Encode(317, 0, 7, false),
                                        mindex::Seed::Encode(481, 0, 11, false),
                                        mindex::Seed::Encode(484, 0, 15, false),
                                        mindex::Seed::Encode(512, 0, 3, true),
                                        mindex::Seed::Encode(840, 0, 10, true),
                                        mindex::Seed::Encode(896, 0, 4, true),
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild4) {
    /*
     * Test building the index from only a selected region.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = true;
    index_params->region_rname = "ecoli";
    index_params->region_rstart = 10;
    index_params->region_rend = 20;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // "TCTGACTGCA" This should be the region.
    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(121, 0, 14, false),
                                        mindex::Seed::Encode(180, 0, 12, true),
                                        mindex::Seed::Encode(301, 0, 13, true),
                                        mindex::Seed::Encode(481, 0, 11, false),
                                        mindex::Seed::Encode(484, 0, 15, false),
                                        mindex::Seed::Encode(840, 0, 10, true),
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild5) {
    /*
     * Test building the index from only a selected region,
     * but the region is out of bounds (hanging off the end of the reference).
     * If the specified region reaches beyond the end of the sequence,
     * the entire suffix of the sequence should be used, and no error should occur.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = true;
    index_params->region_rname = "ecoli";
    index_params->region_rstart = 60;
    index_params->region_rend = 80;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(146, 0, 64, false),
                                        mindex::Seed::Encode(201, 0, 62, false),
                                        mindex::Seed::Encode(461, 0, 61, true),
                                        mindex::Seed::Encode(585, 0, 65, false),
                                        mindex::Seed::Encode(804, 0, 63, false),
                                        mindex::Seed::Encode(820, 0, 60, true),
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild6) {
    /*
     * Test building the index from only a selected region,
     * but the region is out of bounds (hanging off the beginning of the reference).
     * If the specified region starts with a value < 0, the 0 coordinate will
     * be automatically used instead.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = true;
    index_params->region_rname = "ecoli";
    index_params->region_rstart = -10;
    index_params->region_rend = 10;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(2, 0, 2, true),
                                        mindex::Seed::Encode(9, 0, 1, true),
                                        mindex::Seed::Encode(39, 0, 0, true),
                                        mindex::Seed::Encode(224, 0, 5, true),
                                        mindex::Seed::Encode(512, 0, 3, true),
                                        mindex::Seed::Encode(896, 0, 4, true),
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild7) {
    /*
     * Test building the index from only a selected region,
     * but the entire region is out of bounds for the reference.
     * The sequence should not be indexed.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = true;
    index_params->region_rname = "ecoli";
    index_params->region_rstart = 100;
    index_params->region_rend = 200;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild8) {
    /*
     * Test building the index from only a selected region,
     * but the reference specified by the region does not exist.
     * No sequences should be indexed in this scenario.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = true;
    index_params->region_rname = "bsub";
    index_params->region_rstart = 10;
    index_params->region_rend = 20;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                    };

    // raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);

}

TEST(DenseIndexTest, TestBuild9) {
    /*
     * Test a larger k-mer size and a minimizer window.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 15;
    index_params->w = 5;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli:0-70",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(9156492, 0, 50, false),
                                        mindex::Seed::Encode(111424439, 0, 20, false),
                                        mindex::Seed::Encode(167726968, 0, 0, false),
                                        mindex::Seed::Encode(222578820, 0, 30, true),
                                        mindex::Seed::Encode(281992686, 0, 25, false),
                                        mindex::Seed::Encode(364792648, 0, 10, true),
                                        mindex::Seed::Encode(507619596, 0, 15, false),
                                        mindex::Seed::Encode(518950656, 0, 35, false),
                                        mindex::Seed::Encode(555614204, 0, 45, true),
                                        mindex::Seed::Encode(664588817, 0, 55, true),
                                        mindex::Seed::Encode(939520212, 0, 40, true),
                                        mindex::Seed::Encode(959258848, 0, 5, true),
                                    };

    raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild10) {
    /*
     * Tests building on a small sequence with w > 1.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 4;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCA",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    ////////////////////////////////
    /// NOTE: The k-mer keys are still he same as in the minimizer case, because
    /// the dense index takes the minimum of the (fwd, rev) k-mer hash, and discards
    // the other one.
    ////////////////////////////////

    // These are all seeds for a w = 1:
    // Sorted by key:
        // mindex::Seed::Encode(2, 0, 2, true),
        // mindex::Seed::Encode(9, 0, 1, true),
        // mindex::Seed::Encode(39, 0, 0, true),
        // mindex::Seed::Encode(56, 0, 6, true),
        // mindex::Seed::Encode(121, 0, 14, false),
        // mindex::Seed::Encode(131, 0, 8, true),
        // mindex::Seed::Encode(180, 0, 12, true),
        // mindex::Seed::Encode(224, 0, 5, true),
        // mindex::Seed::Encode(288, 0, 9, true),
        // mindex::Seed::Encode(301, 0, 13, true),
        // mindex::Seed::Encode(317, 0, 7, false),
        // mindex::Seed::Encode(481, 0, 11, false),
        // mindex::Seed::Encode(484, 0, 15, false),
        // mindex::Seed::Encode(512, 0, 3, true),
        // mindex::Seed::Encode(840, 0, 10, true),
        // mindex::Seed::Encode(896, 0, 4, true),

    // Sorted by start position. Numbers on the right are the
    // IDs of the windows in which the particular seed is a minimizer.
        // mindex::Seed::Encode(39, 0, 0, true),
        // mindex::Seed::Encode(9, 0, 1, true),
        // mindex::Seed::Encode(2, 0, 2, true),        (0, 1, 2)
        // mindex::Seed::Encode(512, 0, 3, true),
        // mindex::Seed::Encode(896, 0, 4, true),
        // mindex::Seed::Encode(224, 0, 5, true),
        // mindex::Seed::Encode(56, 0, 6, true),      (3, 4, 5, 6)
        // mindex::Seed::Encode(317, 0, 7, false),
        // mindex::Seed::Encode(131, 0, 8, true),     (7, 8)
        // mindex::Seed::Encode(288, 0, 9, true),
        // mindex::Seed::Encode(840, 0, 10, true),
        // mindex::Seed::Encode(481, 0, 11, false),
        // mindex::Seed::Encode(180, 0, 12, true),     (9, 10)
        // mindex::Seed::Encode(301, 0, 13, true),
        // mindex::Seed::Encode(121, 0, 14, false),    (11, 12)
        // mindex::Seed::Encode(484, 0, 15, false),


    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(39, 0, 0, true),
                                        mindex::Seed::Encode(131, 0, 8, true),
                                        mindex::Seed::Encode(180, 0, 12, true),
                                        mindex::Seed::Encode(896, 0, 4, true),
                                    };

    // for (int32_t i = 0; i < index->seeds().size(); i++) {
    //     // std::cerr << "[" << i << "] " << mindex::Seed(index->seeds()[i]).Verbose() << std::endl;
    //     auto minimizer = mindex::Seed(index->seeds()[i]);
    //     std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    //     // std::cerr << "[" << i << "] " << mindex::Seed(expected[i]).Verbose() << std::endl;
    // }

    raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuild11) {
    /*
     * Tests building on a small sequence with w > 1.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 4;
    index_params->homopolymer_suppression = false;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCANNNACTNNNNNAGCTTTTCATTCTGACTGCA",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    // This is different compared to the minimizer index, because we are keeping only
    // those seeds whose position is (pos % w) == 0.
    std::vector<mindex128_t> expected{
                                        mindex::Seed::Encode(9, 0, 32, true),
                                        mindex::Seed::Encode(39, 0, 0, true),
                                        mindex::Seed::Encode(131, 0, 8, true),
                                        mindex::Seed::Encode(180, 0, 12, true),
                                        mindex::Seed::Encode(224, 0, 36, true),
                                        mindex::Seed::Encode(288, 0, 40, true),
                                        mindex::Seed::Encode(301, 0, 44, true),
                                        mindex::Seed::Encode(896, 0, 4, true),
                                    };

    // for (int32_t i = 0; i < index->seeds().size(); i++) {
    //     // std::cerr << "[" << i << "] " << mindex::Seed(index->seeds()[i]).Verbose() << std::endl;
    //     auto minimizer = mindex::Seed(index->seeds()[i]);
    //     std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    //     // std::cerr << "[" << i << "] " << mindex::Seed(expected[i]).Verbose() << std::endl;
    // }

    raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestCollectHits1) {
    /*
     * Test an empty query.
    */

    // Create the parameters.
    auto index_params = mindex::createIndexParams();

    // Create the index.
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Define the query sequence.
    std::string query("");

    std::vector<mindex::SeedHitPacked> hits = index->CollectHits(query);

    ASSERT_EQ(hits.size(), 0);
}

TEST(DenseIndexTest, TestCollectHits2) {
    /*
     * Test an empty index but a non-empty query.
    */

    // Create the parameters.
    auto index_params = mindex::createIndexParams();

    // Create the index.
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    };
    std::vector<std::string> headers = {
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Define the query sequence.
    std::string query("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC");

    std::vector<mindex::SeedHitPacked> hits = index->CollectHits(query);

    ASSERT_EQ(hits.size(), 0);
}

TEST(DenseIndexTest, TestCollectHits3) {
    // Test the collection of hits with the default parameters (k = 15, w = 5).
    // This is like a realistic scenario, but the query and the target are the same.
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    // Create the parameters.
    auto index_params = mindex::createIndexParams();

    // Create the index.
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Define the query sequence.
    std::string query = seqs[0];

    // for (int32_t i = 0; i < index->seeds().size(); i++) {
    //     auto minimizer = mindex::Seed(index->seeds()[i]);
    //     std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    // }

    std::vector<mindex::SeedHitPacked> hits = index->CollectHits(query);

    ASSERT_EQ(hits.size(), index->seeds().size());

    for (size_t i = 0; i < hits.size(); i++) {
        auto& hit = hits[i];
        // std::cerr << "[hit " << i << "] " << hit.Verbose() << std::endl;
        ASSERT_EQ(hit.QueryPos(), hit.TargetPos());
    }
}

TEST(DenseIndexTest, TestCollectHits4) {
    /*
     * Test a scenario which results in 0 hits, but neither the
     * index nor the query are empty.
    */

    // Create the parameters.
    auto index_params = mindex::createIndexParams();

    // Create the index.
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Define the query sequence.
    std::string query("AAAAAAAAAAAAAAAAAAAA");

    // for (int32_t i = 0; i < index->seeds().size(); i++) {
    //     auto minimizer = mindex::Seed(index->seeds()[i]);
    //     std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    // }

    std::vector<mindex::SeedHitPacked> hits = index->CollectHits(query);

    ASSERT_EQ(hits.size(), 0);

    // ASSERT_EQ(hits.size(), index->seeds().size());

    // for (size_t i = 0; i < hits.size(); i++) {
    //     auto& hit = hits[i];
    //     // std::cerr << "[hit " << i << "] " << hit.Verbose() << std::endl;
    //     ASSERT_EQ(hit.QueryPos(), hit.TargetPos());
    // }

}

TEST(DenseIndexTest, TestCollectHits5) {
    /*
     * Test a scenario where query is not identical to the index.
    */

    // Create the parameters.
    auto index_params = mindex::createIndexParams();

    // Create the index.
    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Define the query sequence.
    //               ("AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC");
    //                                           X
    //               ("AGCTTTTCATTCTGACTGCAACGGGCGATATGTCTCTGTGTGG-TTAAAAAAAGAGTGTCTGATAGCAGC");
    std::string query("AGCTTTTCATTCTGACTGCAACGGGCGATATGTCTCTGTGTGGTTAAAAAAAGAGTGTCTGATAGCAGC");

    std::vector<mindex::SeedHitPacked> hits = index->CollectHits(query);

    std::vector<mindex::SeedHitPacked> expected {
                                                mindex::SeedHitPacked(0, false, 0, 0, 0),
                                                mindex::SeedHitPacked(0, false, 5, 0, 5),
                                                mindex::SeedHitPacked(0, false, 10, 0, 10),
                                                mindex::SeedHitPacked(0, false, 45, 0, 44),
                                                mindex::SeedHitPacked(0, false, 50, 0, 49),
                                                mindex::SeedHitPacked(0, false, 55, 0, 54),
                                                       };

    // for (int32_t i = 0; i < index->seeds().size(); i++) {
    //     auto minimizer = mindex::Seed(index->seeds()[i]);
    //     std::cerr << "mindex::Seed::Encode(" << minimizer.key << ", " << minimizer.seq_id << ", " << minimizer.pos << ", " << ((minimizer.flag) ? "true" : "false") << ")," << std::endl;
    // }
    // for (size_t i = 0; i < hits.size(); i++) {
    //     auto& hit = hits[i];
    //     std::cerr << "mindex::SeedHitPacked::PackTo128t(" << hit.TargetId() << ", " << hit.TargetPos() << ", " << hit.QueryPos() << ")," << std::endl;
    // }

    ASSERT_EQ(hits, expected);
}

TEST(DenseIndexTest, TestBuildHomopolymerSuppression1) {
    /*
     * Tests building on a more complex sequence.
    */
    auto index_params = mindex::createIndexParams();

    index_params->k = 5;
    index_params->w = 1;
    index_params->homopolymer_suppression = true;
    index_params->max_homopolymer_len = 5;
    index_params->freq_percentil = 0.0;
    index_params->min_occ_cutoff = 0;
    index_params->is_region_specified = false;
    index_params->region_rname = "";
    index_params->region_rstart = 0;
    index_params->region_rend = 0;
    index_params->index_only_fwd_strand = false;

    auto index = mindex::createDenseIndex(index_params);

    // Define the reference sequences and their headers.
    std::vector<std::string> seqs = {
                                    "ACGTTTTG",
                                    };
    std::vector<std::string> headers = {
                                    "ecoli:0-20",
                                    };

    // Add them to the index. This only adds the data, but does not build the index.
    index->AddSequences(seqs, headers);

    // Build the index from the added sequences.
    index->BuildIndex();

    // Seeds are in a sorted order by their key.
    std::vector<mindex128_t> expected{
                                    mindex::Seed::Encode(110, 0, 0, false),
                                    };

    raptor::unit::VerboseSeeds(index->seeds(), expected);

    ASSERT_EQ(index->seeds(), expected);
}

TEST(DenseIndexTest, TestBuildHomopolymerSuppression2) {
    /*
     * Tests building a more complex sequence with homopolymer suppression.
     * The index should correspond to the index same sequence without any homopolymers.
    */

    // Index params and index for homopolymer suppresed sequence (unit under test).
    auto index_params_hp = mindex::createIndexParams();
    index_params_hp->k = 5;
    index_params_hp->w = 1;
    index_params_hp->homopolymer_suppression = true;
    index_params_hp->max_homopolymer_len = 5;
    index_params_hp->freq_percentil = 0.0;
    index_params_hp->min_occ_cutoff = 0;
    index_params_hp->is_region_specified = false;
    index_params_hp->region_rname = "";
    index_params_hp->region_rstart = 0;
    index_params_hp->region_rend = 0;
    index_params_hp->index_only_fwd_strand = false;
    auto index_hp = mindex::createDenseIndex(index_params_hp);
    // Define the reference sequences and their headers.
    std::vector<std::string> seqs_for_hp = {"AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"};
    std::vector<std::string> headers_for_hp = {"ecoli"};
    // Add them to the index. This only adds the data, but does not build the index.
    index_hp->AddSequences(seqs_for_hp, headers_for_hp);
    // Build the index from the added sequences.
    index_hp->BuildIndex();

    // Index params for deactivated homopolymer suppression and reference sequences
    // which have no homopolymers. This should result in the same index as the above.
    auto index_params_no_hp = mindex::createIndexParams();
    index_params_no_hp->k = 5;
    index_params_no_hp->w = 1;
    index_params_no_hp->homopolymer_suppression = false;
    index_params_no_hp->max_homopolymer_len = 5;
    index_params_no_hp->freq_percentil = 0.0;
    index_params_no_hp->min_occ_cutoff = 0;
    index_params_no_hp->is_region_specified = false;
    index_params_no_hp->region_rname = "";
    index_params_no_hp->region_rstart = 0;
    index_params_no_hp->region_rend = 0;
    index_params_no_hp->index_only_fwd_strand = false;
    auto index_no_hp = mindex::createDenseIndex(index_params_no_hp);
    // Define the reference sequences and their headers.
    std::vector<std::string> seqs_with_no_hp = {"AGCTCATCTGACTGCACGCATATGTCTCTGTGTGATAGAGTGTCTGATAGCAGC"};
    std::vector<std::string> headers_with_no_hp = {"ecoli"};
    // Add them to the index. This only adds the data, but does not build the index.
    index_no_hp->AddSequences(seqs_with_no_hp, headers_with_no_hp);
    // Build the index from the added sequences.
    index_no_hp->BuildIndex();

    // We need to manually compare the seed keys, because the positions will have changed.
    ASSERT_EQ(index_hp->seeds().size(), index_no_hp->seeds().size());

    raptor::unit::VerboseSeeds(index_hp->seeds(), index_no_hp->seeds());

    for (int32_t i = 0; i < index_hp->seeds().size(); i++) {
        auto minimizer_hp = mindex::Seed(index_hp->seeds()[i]);
        auto minimizer_no_hp = mindex::Seed(index_no_hp->seeds()[i]);
        ASSERT_EQ(minimizer_hp.key, minimizer_no_hp.key);
        ASSERT_EQ(minimizer_hp.flag, minimizer_no_hp.flag);
    }
}

TEST(DenseIndexTest, TestLoad1) {
    // ASSERT_EQ(0, 1);
}

TEST(DenseIndexTest, TestStore1) {
    // ASSERT_EQ(0, 1);
}
