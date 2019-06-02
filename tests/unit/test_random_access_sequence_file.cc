#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <sequences/sequence_file.h>
#include <sequences/sequence_file_composite_fofn.h>
#include <sequences/random_access_sequence_file.h>
#include <sequences/sequence_serializer.h>
#include <utility/memtime.h>
#include <raptor_fetch/params_raptor_fetch.h>

TEST(RandomAccessSequenceFile, FetchSequenceUsingId) {
    // The RaptorDB was generated using: "$ raptor-reshape -r test-data/ecoli-small/reads.6x.fwd.fasta -o test-data/raptordb-fetch/test-1 --split-blocks --block-size 0.1".
    std::string in_path("test-data/raptordb-fetch/test-1.rdb");
    std::string in_exp_path("test-data/ecoli-small/reads.6x.fwd.fasta");

    int64_t max_open_files = 50;
	auto random_seq_file = mindex::createRandomAccessSequenceFile(in_path, max_open_files);

	auto seq_file_parser = mindex::createSequenceFileCompositeFofn({in_exp_path}, mindex::SequenceFormat::Auto);
	auto seq_file = seq_file_parser->YieldAll();

    for (const auto& seq: seq_file->seqs()) {
        auto result_seq = random_seq_file->FetchSequence(seq->abs_id());

        // std::cerr << "Testing: abs_id = " << seq->abs_id() << ", name: '" << seq->header() << "'\n";

        ASSERT_NE(result_seq, nullptr);
        ASSERT_EQ(result_seq->header(), seq->header());
        ASSERT_EQ(result_seq->GetSequenceAsString(), seq->GetSequenceAsString());
        ASSERT_EQ(result_seq->GetQualityAsString(), seq->GetQualityAsString());
    }
}

TEST(RandomAccessSequenceFile, FetchSequenceUsingQname) {
    std::string in_path("test-data/raptordb-fetch/test-1.rdb");
    std::string in_exp_path("test-data/ecoli-small/reads.6x.fwd.fasta");

    int64_t max_open_files = 50;
	auto random_seq_file = mindex::createRandomAccessSequenceFile(in_path, max_open_files);

	auto seq_file_parser = mindex::createSequenceFileCompositeFofn({in_exp_path}, mindex::SequenceFormat::Auto);
	auto seq_file = seq_file_parser->YieldAll();

    for (const auto& seq: seq_file->seqs()) {
        auto result_seq = random_seq_file->FetchSequence(seq->header());

        ASSERT_NE(result_seq, nullptr);
        ASSERT_EQ(result_seq->header(), seq->header());
        ASSERT_EQ(result_seq->GetSequenceAsString(), seq->GetSequenceAsString());
        ASSERT_EQ(result_seq->GetQualityAsString(), seq->GetQualityAsString());
    }
}

TEST(RandomAccessSequenceFile, FetchSequenceUsingQnameSingleFile) {
    std::string in_path("test-data/raptordb-fetch/test-2.rdb");
    std::string in_exp_path("test-data/ecoli-small/reads.6x.fwd.fasta");

    int64_t max_open_files = 50;
	auto random_seq_file = mindex::createRandomAccessSequenceFile(in_path, max_open_files);

	auto seq_file_parser = mindex::createSequenceFileCompositeFofn({in_exp_path}, mindex::SequenceFormat::Auto);
	auto seq_file = seq_file_parser->YieldAll();

    for (const auto& seq: seq_file->seqs()) {
        auto result_seq = random_seq_file->FetchSequence(seq->header());

        ASSERT_NE(result_seq, nullptr);
        ASSERT_EQ(result_seq->header(), seq->header());
        ASSERT_EQ(result_seq->GetSequenceAsString(), seq->GetSequenceAsString());
        ASSERT_EQ(result_seq->GetQualityAsString(), seq->GetQualityAsString());
    }
}
