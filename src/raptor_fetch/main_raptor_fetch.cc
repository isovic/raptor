#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
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

void RunTool(std::shared_ptr<raptor::ParamsRaptorFetch> params) {
	LOG_ALL("Raptor-fetch starts processing.\n");

    if (params->in_paths.size() == 0) {
        FATAL_REPORT(ERR_OPENING_FILE, "No input files are specified. Exiting.");
    }

    auto seq_file = mindex::createRandomAccessSequenceFile(params->rdb_path, 50);

    if (params->job_str == "erc") {
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
                            std::cerr << "Warning could not find seq_a, name: '" << qda->name << "'. Skipping all seq_a overlaps.\n";
                            break;
                        }
                        if (ovl->a_rev()) {
                            seq_a->ReverseComplement();
                        }
                        fprintf (stdout, "%09lld\t%s\n", seq_a->abs_id(), seq_a->GetSequenceAsString().c_str());
                        a_written = true;
                    }
                    auto seq_b = seq_file->FetchSequence(qdb->name);
                    if (seq_b == nullptr) {
                        std::cerr << "Warning could not find seq_b, name: '" << qdb->name << "'. Skipping the overlap.\n";
                        break;
                    }
                    if (ovl->b_rev()) {
                        seq_b->ReverseComplement();
                    }
                    fprintf (stdout, "%09lld\t%s\n", seq_b->abs_id(), seq_b->GetSequenceAsString().c_str());
                }

                std::cout << "+ +\n";
            }

            span_start += span;
        }
    } else if (params->job_str == "clip") {
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
                std::cerr << "Warning: could not find sequence: '" << tokens[0] << "' in the RaptorDB. Skipping.\n";
                continue;
            }
            fprintf (stdout, ">%09lld\n%s\n", seq->abs_id(), seq->GetSubsequenceAsString(start, end).c_str());

        }
    } else if (params->job_str == "fetch") {
        std::ifstream ifs(params->in_paths[0]);
        std::string line;
        while(std::getline(ifs, line)) {
            std::vector<std::string> tokens = raptor::TokenizeToWhitespaces(line);
            if (tokens.size() == 0) {
                continue;
            }
            auto seq = seq_file->FetchSequence(tokens[0]);
            if (seq == nullptr) {
                std::cerr << "Warning: could not find sequence: '" << tokens[0] << "' in the RaptorDB. Skipping.\n";
                continue;
            }
            auto header_tokens = raptor::TokenizeToWhitespaces(seq->header());
            fprintf (stdout, ">%s\n%s\n", header_tokens[0].c_str(), seq->GetSequenceAsString().c_str());
        }

    } else {
        std::cerr << "ERROR: Unknown job_str: '" << params->job_str << "'. Exiting.\n";
        return;
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

	RunTool(parameters);

	return 0;
}

int main(int argc, char **argv) {
  	return Run(argc, argv);
}

#endif
