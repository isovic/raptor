#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <sequences/sequence_file_parser_sam.h>

TEST(SequenceFileParserSam, YieldSequence) {
    auto parser = mindex::createSequenceFileParserSam("test-data/sequence-parser/test5.sam");
    mindex::SequencePtr seq;

    std::vector<std::string> result_headers;
    std::vector<std::string> result_seqs;
    std::vector<std::string> result_quals;

    while ((seq = parser->YieldSequence(true)) != nullptr) {
        result_headers.emplace_back(seq->header());
        result_seqs.emplace_back(seq->GetSequenceAsString());
        result_quals.emplace_back(seq->GetQualityAsString());
    }

    std::string expected_file_header = {
        "@HD\tVN:1.5\tSO:unknown\n"
        "@SQ\tSN:Chr01\tLN:1000\n"
    };
    std::vector<std::string> expected_headers = {
        "seq1", "seq2", "seq3", "seq4", "seq5", "seq6"
    };
    std::vector<std::string> expected_seqs = {
        "AAAC", "GGGC", "TTTTTAAAAATTTT", "AACTGACTGACTGAATGAACGAACT", "AAAC", "AAAC"
    };
    std::vector<std::string> expected_quals = {
        "", "", "", "!!!!!11111222223333344444", "1234", "1234"
    };

    ASSERT_EQ(result_seqs.size(), expected_seqs.size());
    ASSERT_EQ(result_headers, expected_headers);
    ASSERT_EQ(result_seqs, expected_seqs);
    ASSERT_EQ(parser->GetFileHeaderAsString(), expected_file_header);
}

TEST(SequenceFileParserSam, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test5.sam";
    std::string expected = {
        "@HD\tVN:1.5\tSO:unknown\n"
        "@SQ\tSN:Chr01\tLN:1000\n"
    };

    auto parser = mindex::createSequenceFileParserSam(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserSam, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test5.sam";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserSam(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserSam, IsOpen1) {
    std::string in_path = "test-data/sequence-parser/test5.sam";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserSam(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

// TEST(SequenceFileParserSam, GetFileOffset1) {
//     // This functionality is integratively tested in FileSeek1.
// }

TEST(SequenceFileParserSam, FileSeek1) {
    auto parser = mindex::createSequenceFileParserSam("test-data/sequence-parser/test5.sam");

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
