#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <sequences/sequence_file_parser_gfa1.h>

TEST(SequenceFileParserGfa1, YieldSequence) {
    // ASSERT_EQ(result, expected);
    auto parser = mindex::createSequenceFileParserGfa1("test-data/sequence-parser/test3.gfa1");
    mindex::SequencePtr seq;

    std::vector<std::string> result_headers;
    std::vector<std::string> result_seqs;
    std::vector<std::string> result_quals;

    while ((seq = parser->YieldSequence(true)) != nullptr) {
        result_headers.emplace_back(seq->header());
        result_seqs.emplace_back(seq->GetSequenceAsString());
        result_quals.emplace_back(seq->GetQualityAsString());
    }

    std::vector<std::string> expected_headers = {
        "seq1", "seq2", "seq3", "seq4"
    };
    std::vector<std::string> expected_seqs = {
        "AAAC", "GGGC", "TTTTTAAAAATTTT", "AACTGACTGACTGAATGAACGAACT"
    };
    std::vector<std::string> expected_quals = {
        "", "", "", ""
    };

    ASSERT_EQ(result_seqs.size(), 4);
    ASSERT_EQ(result_headers, expected_headers);
    ASSERT_EQ(result_seqs, expected_seqs);
}

TEST(SequenceFileParserGfa1, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected;

    auto parser = mindex::createSequenceFileParserGfa1(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserGfa1, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserGfa1(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserGfa1, IsOpen1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserGfa1(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

// TEST(SequenceFileParserGfa1, GetFileOffset1) {
//     // This functionality is integratively tested in FileSeek1.
// }

TEST(SequenceFileParserGfa1, FileSeek1) {
    auto parser = mindex::createSequenceFileParserGfa1("test-data/sequence-parser/test3.gfa1");

    ASSERT_NE(parser, nullptr);
    ASSERT_EQ(parser->IsOpen(), true);

    std::vector<int64_t> offsets;
    std::vector<mindex::SequencePtr> seqs;

    mindex::SequencePtr seq;
    int64_t prev_offset = parser->GetFileOffset();

    // Create a mini index of the sequences.
    while ((seq = parser->YieldSequence(true)) != nullptr) {
        seqs.emplace_back(std::move(seq));
        offsets.emplace_back(prev_offset);
        prev_offset = parser->GetFileOffset();
    }

    // Create a reversed permutation of sequence IDs for reading.
    std::vector<int64_t> permutation;
    for (int64_t i = 0; i < static_cast<int64_t>(seqs.size()); ++i) {
        permutation.emplace_back(i);
    }
    std::reverse(permutation.begin(), permutation.end());

    // Lookup all the same sequences, in a random order.
    for (const auto& seq_id: permutation) {
        parser->FileSeek(offsets[seq_id]);
        seq = parser->YieldSequence(true);
        ASSERT_NE(seq, nullptr);
        ASSERT_EQ(seq->header(), seqs[seq_id]->header());
        ASSERT_EQ(seq->GetSequenceAsString(), seqs[seq_id]->GetSequenceAsString());
        ASSERT_EQ(seq->GetQualityAsString(), seqs[seq_id]->GetQualityAsString());
    }
}

//    int64_t GetFileOffset() const;
//    bool FileSeek(int64_t pos);
