#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <sequences/sequence_file_parser_gfa2.h>

TEST(SequenceFileParserGfa2, YieldSequence) {
    // ASSERT_EQ(result, expected);
    auto parser = mindex::createSequenceFileParserGfa2("test-data/sequence-parser/test4.gfa2");
    mindex::SequencePtr seq;

    std::vector<std::string> result_headers;
    std::vector<std::string> result_seqs;
    std::vector<std::string> result_quals;

    while ((seq = parser->YieldSequence()) != nullptr) {
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

TEST(SequenceFileParserGfa2, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected;

    auto parser = mindex::createSequenceFileParserGfa2(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserGfa2, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserGfa2(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserGfa2, IsOpen1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserGfa2(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserGfa2, GetFileOffset1) {

}

TEST(SequenceFileParserGfa2, FileSeek1) {

}

//    int64_t GetFileOffset() const;
//    bool FileSeek(int64_t pos);
