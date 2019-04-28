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

    while ((seq = parser->YieldSequence()) != nullptr) {
        result_headers.emplace_back(seq->header());
        result_seqs.emplace_back(seq->GetSequenceAsString());
        result_quals.emplace_back(seq->GetQualityAsString());
    }

    std::string expected_file_header = {
        "@HD\tVN:1.5\tSO:unknown"
        "@SQ\tSN:seq1\tLN:4"
        "@SQ\tSN:seq2\tLN:4"
        "@SQ\tSN:seq3\tLN:14"
        "@SQ\tSN:seq4\tLN:25"
    };
    std::vector<std::string> expected_headers = {
        "seq1", "seq2", "seq3", "seq4"
    };
    std::vector<std::string> expected_seqs = {
        "AAAC", "GGGC", "TTTTTAAAAATTTT", "AACTGACTGACTGAATGAACGAACT"
    };
    std::vector<std::string> expected_quals = {
        "", "", "", "!!!!!11111222223333344444"
    };

    ASSERT_EQ(result_seqs.size(), 4);
    ASSERT_EQ(result_headers, expected_headers);
    ASSERT_EQ(result_seqs, expected_seqs);
    ASSERT_EQ(parser->GetFileHeaderAsString(), expected_file_header);
}

TEST(SequenceFileParserSam, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected;

    auto parser = mindex::createSequenceFileParserSam(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserSam, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserSam(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserSam, IsOpen1) {
    std::string in_path = "test-data/sequence-parser/test3.gfa1";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserSam(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserSam, GetFileOffset1) {

}

TEST(SequenceFileParserSam, FileSeek1) {

}

//    int64_t GetFileOffset() const;
//    bool FileSeek(int64_t pos);
