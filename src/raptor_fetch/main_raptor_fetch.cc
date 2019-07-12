#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <raptor/args_raptor.h>
#include <limits>
#include <sequences/sequence_file.h>
#include <sequences/random_access_sequence_file.h>
#include <sequences/sequence_serializer.h>
#include <utility/memtime.h>
#include <utility/stringutil.h>
#include <raptor_fetch/params_raptor_fetch.h>
#include <raptor_fetch/overlaps/overlap_file.h>
#include <raptor_fetch/args_raptor_fetch.h>

namespace raptor {
namespace fetch {

enum class OutFormat {
    FASTA,
    ERC
};

inline void WriteSequence(std::ostream& os, const mindex::SequencePtr& seq, bool use_id_for_output, const OutFormat& out_fmt);
void WriteSubsequence(std::ostream& os, const mindex::SequencePtr& seq, bool use_id_for_output, const OutFormat& out_fmt, int64_t start, int64_t end);


inline void WriteSequence(std::ostream& os, const mindex::SequencePtr& seq, bool use_id_for_output, const OutFormat& out_fmt) {
    WriteSubsequence(os, seq, use_id_for_output, out_fmt, 0, seq->len());
}

void WriteSubsequence(std::ostream& os, const mindex::SequencePtr& seq, bool use_id_for_output, const OutFormat& out_fmt, int64_t start, int64_t end) {
    if (out_fmt == OutFormat::FASTA) {
        os << ">";
    }
    if (use_id_for_output) {
        os << std::setw(9) << std::setfill('0') << seq->abs_id();
    } else {
        os << raptor::TrimToFirstWhiteSpace(seq->header());
    }
    if (out_fmt == OutFormat::FASTA) {
        os << std::endl;
    } else if (out_fmt == OutFormat::ERC) {
        os << "\t";
    } else {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Unknown output format.");
    }
    os << seq->GetSubsequenceAsString(start, end)<< "\n";
}

void RunToolERC(std::shared_ptr<raptor::ParamsRaptorFetch> params, std::ostream& os) {
    auto seq_file = mindex::createRandomAccessSequenceFile(params->rdb_path, 50);

    // Create a reader.
    auto overlap_file = raptor::createOverlapFile(params->in_paths[0]);

    // Load the entire file into memory.
    overlap_file->LoadNextBatch(0, false, params->min_len, params->min_score, params->min_idt);

    auto& batch = overlap_file->batch();

    int64_t span = 0, span_start = 0;
    while ((span = raptor::FindQuerySpan(batch, span_start))) {

        if (span >= params->min_cov) {
            bool a_written = false;
            for (int64_t oid = span_start; oid < (span_start + span); ++oid) {
                auto& ovl = batch[oid];
                auto& qda = overlap_file->GetQueryDataById(ovl->a_id());
                auto& qdb = overlap_file->GetQueryDataById(ovl->b_id());
                if (qda == nullptr || qdb == nullptr) {
                    continue;
                }
                if (qda->IsFiltered() || qdb->IsFiltered() || ovl->IsFiltered()) {
                    continue;
                }

                if (a_written == false) {
                    auto seq_a = seq_file->FetchSequence(qda->name);
                    if (seq_a == nullptr) {
                        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Could not find seq_a, name: '%s'. Skipping all seq_a overlaps.", qda->name.c_str());
                        break;
                    }
                    if (ovl->a_rev()) {
                        seq_a->ReverseComplement();
                    }
                    raptor::fetch::WriteSequence(os, seq_a, params->use_id_for_output, raptor::fetch::OutFormat::ERC);
                    a_written = true;
                }
                auto seq_b = seq_file->FetchSequence(qdb->name);
                if (seq_b == nullptr) {
                    WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Could not find seq_b, name: '%s'. Skipping all seq_a overlaps.", qdb->name.c_str());
                    break;
                }
                if (ovl->b_rev()) {
                    seq_b->ReverseComplement();
                }

                raptor::fetch::WriteSequence(os, seq_b, params->use_id_for_output, raptor::fetch::OutFormat::ERC);
            }

            os << "+ +\n";
        }

        span_start += span;
    }
}

void RunToolClip(std::shared_ptr<raptor::ParamsRaptorFetch> params, std::ostream& os) {
    auto seq_file = mindex::createRandomAccessSequenceFile(params->rdb_path, 50);
    std::ifstream ifs(params->in_paths[0]);
    std::string line;
    while(std::getline(ifs, line)) {
        std::vector<std::string> tokens = raptor::TokenizeToWhitespaces(line);
        if (tokens.size() < 5) {
            continue;
        }
        if (tokens[4] != "0") {
            // Check if the query was filtered. 0 means all good.
            continue;
        }
        int32_t start = std::stoi(tokens[1]);
        int32_t end = std::stoi(tokens[2]);
        auto seq = seq_file->FetchSequence(tokens[0]);
        if (seq == nullptr) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Could not find find sequence: '%s' in the RaptorDB. Skipping.", tokens[0].c_str());
            continue;
        }

        raptor::fetch::WriteSequence(os, seq, params->use_id_for_output, raptor::fetch::OutFormat::FASTA);
    }
}

void RunToolFetch(std::shared_ptr<raptor::ParamsRaptorFetch> params, std::ostream& os) {
    auto seq_file = mindex::createRandomAccessSequenceFile(params->rdb_path, 50);
    std::ifstream ifs(params->in_paths[0]);
    std::string line;
    while(std::getline(ifs, line)) {
        std::vector<std::string> tokens = raptor::TokenizeToWhitespaces(line);
        if (tokens.size() == 0) {
            continue;
        }
        auto seq = seq_file->FetchSequence(tokens[0]);
        if (seq == nullptr) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Could not find find sequence: '%s' in the RaptorDB. Skipping.", tokens[0].c_str());
            continue;
        }
        auto header_tokens = raptor::TokenizeToWhitespaces(seq->header());

        raptor::fetch::WriteSequence(os, seq, params->use_id_for_output, raptor::fetch::OutFormat::FASTA);
    }
}

}
}

void RunTool(std::shared_ptr<raptor::ParamsRaptorFetch> params, std::ostream& os) {
	LOG_ALL("Raptor-fetch starts processing.\n");

    if (params->in_paths.size() == 0) {
        FATAL_REPORT(ERR_OPENING_FILE, "No input files are specified. Exiting.");
    }

    if (params->job_str == "erc") {
        raptor::fetch::RunToolERC(params, os);

    } else if (params->job_str == "clip") {
        raptor::fetch::RunToolClip(params, os);

    } else if (params->job_str == "fetch") {
        raptor::fetch::RunToolFetch(params, os);

    } else {
		ERROR_REPORT(ERR_UNEXPECTED_VALUE, "Unknown job_str: '%s'.", params->job_str.c_str());
    }

	LOG_ALL("Memory usage: %.2f GB\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
	LOG_ALL("Done!\n");
}

int Run(int argc, char **argv) {

	std::string program_name(argv[0]);

	auto parameters = raptor::createParamsRaptorFetch();

	if (raptor::ProcessArgsRaptorFetch(argc, argv, parameters)) {
		return 1;
	}

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters->verbose_level);
	if (parameters->verbose_level == 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
	} else if (parameters->verbose_level > 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
	}
	fflush(stdout);

	RunTool(parameters, std::cout);

	return 0;
}

int main(int argc, char **argv) {
  	return Run(argc, argv);
}

#endif
