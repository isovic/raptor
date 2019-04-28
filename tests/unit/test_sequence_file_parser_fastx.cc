#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <sequences/sequence_file_parser_fastx.h>

TEST(SequenceFileParserFastx, YieldSequence1) {
    // ASSERT_EQ(result, expected);
    auto parser = mindex::createSequenceFileParserFastx("test-data/sequence-parser/test1.fasta");
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
        "seq1", "seq2", "seq3 additional_value", "seq4"
    };
    std::vector<std::string> expected_seqs = {
        "AAAC", "GGGC", "TTTTTAAAAATTTT", "AACTGACTGACTGAATGAACGAACT"
    };
    std::vector<std::string> expected_quals = {
        "", "", "", ""
    };

    ASSERT_EQ(result_headers, expected_headers);
    ASSERT_EQ(result_seqs, expected_seqs);
}

TEST(SequenceFileParserFastx, YieldSequence2) {
    // ASSERT_EQ(result, expected);
    auto parser = mindex::createSequenceFileParserFastx("test-data/sequence-parser/test2.fastq");

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
        "seq1", "seq2", "seq3 additional_value", "seq4"
    };
    std::vector<std::string> expected_seqs = {
        "AAAC", "GGGC", "TTTTTAAAAATTTT", "AACTGACTGACTGAATGAACGAACT"
    };
    std::vector<std::string> expected_quals = {
        "!!!!", "AAAA", "CCCCCCCCCCCCCC", "!!!!!11111222223333344444"
    };

    ASSERT_EQ(result_headers, expected_headers);
    ASSERT_EQ(result_seqs, expected_seqs);
    ASSERT_EQ(result_quals, expected_quals);
}

TEST(SequenceFileParserFastx, GetFileHeaderAsString1) {
    std::string in_path = "test-data/sequence-parser/test1.fasta";
    std::string expected;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, GetFileHeaderAsString2) {
    std::string in_path = "test-data/sequence-parser/test2.fastq";
    std::string expected;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    std::string result = parser->GetFileHeaderAsString();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, GetFilePath1) {
    std::string in_path = "test-data/sequence-parser/test1.fasta";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, GetFilePath2) {
    std::string in_path = "test-data/sequence-parser/test2.fastq";
    std::string expected = in_path;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    std::string result = parser->GetFilePath();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, IsOpen1) {
    std::string in_path = "test-data/sequence-parser/test1.fasta";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, IsOpen2) {
    std::string in_path = "test-data/sequence-parser/test2.fastq";
    bool expected = true;

    auto parser = mindex::createSequenceFileParserFastx(in_path);
    bool result = parser->IsOpen();

    ASSERT_EQ(result, expected);
}

TEST(SequenceFileParserFastx, GetFileOffset1) {

}

TEST(SequenceFileParserFastx, FileSeek1) {

}
