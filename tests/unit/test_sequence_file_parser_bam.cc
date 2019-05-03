#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <sequences/sequence_file_parser_bam.h>

TEST(SequenceFileParserBam, YieldSequence) {
    auto parser = mindex::createSequenceFileParserBam("test-data/sequence-parser/test6.bam");
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
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
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

TEST(SequenceFileParserBam, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test6.bam";
    std::string expected_file_header = {
        "@HD\tVN:1.5\tSO:unknown\tpb:3.0.7\n"
        "@SQ\tSN:Chr01\tLN:1000\n"
    };

    auto parser = mindex::createSequenceFileParserBam(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected_file_header);
}

TEST(SequenceFileParserBam, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test6.bam";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserBam(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

// TEST(SequenceFileParserBam, IsOpen1) {
//     std::string in_path = "test-data/sequence-parser/test6.bam";
//     bool expected = true;

//     auto parser = mindex::createSequenceFileParserBam(in_path);
//     bool result = parser->IsOpen();

//     ASSERT_EQ(result, expected);
// }

TEST(SequenceFileParserBam, GetFileOffset1) {

}

TEST(SequenceFileParserBam, FileSeek1) {

}

//    int64_t GetFileOffset() const;
//    bool FileSeek(int64_t pos);
