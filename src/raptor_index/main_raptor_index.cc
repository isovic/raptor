#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <params/params_raptor_index.h>
#include <raptor/index_factory.h>
#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
#include <fstream>
#include <raptor/args_raptor.h>
#include <sequences/sequence_file.h>
#include <utility/memtime.h>

void RunRaptorIndex(std::shared_ptr<raptor::ParamsRaptorIndex> parameters) {
	/*
	* Create the index.
	*/
	LOG_ALL("Creating the index.\n");

    auto index = mindex::createMinimizerIndex(parameters->index_params);

	int64_t batch_size_int64 = static_cast<int64_t>(parameters->batch_size);
	const int64_t batch_size_in_bytes = static_cast<int64_t>(parameters->batch_size * 1024 * 1024);

	for (auto& ref_path: parameters->ref_paths) {
		auto seqs = mindex::createSequenceFile(ref_path, mindex::SequenceFormat::Auto);

		while (seqs->LoadBatchMB(batch_size_int64)) {
			for (const auto& seq: seqs->seqs()) {
				int64_t index_size = index->total_len();
				if ((index_size + seq->len()) > batch_size_in_bytes) {
					// Store the batch to disk.
					LOG_ALL("Storing a batch index to disk.\n");
				    index = mindex::createMinimizerIndex(parameters->index_params);
					index->BuildIndex();
					LOG_ALL("Reloaded an empty index.\n");
				}
				index->AddSequence(seq->data(), seq->header());
				LOG_ALL("index->total_len() = %ld, batch_size_in_bytes = %ld, batch_size_int64 = %ld\n", index->total_len(), batch_size_in_bytes, batch_size_int64);
			}

			LOG_ALL("Batch done! Memory usage: %.2f GB\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
			fflush(stdout);
		}
	}

	if (index != nullptr && index->total_len() > 0) {
		// Store the batch to disk.
		index->BuildIndex();
		LOG_ALL("Storing the final batch index to disk.\n");
	}

	LOG_ALL("Memory usage: %.2f GB\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
	LOG_ALL("Done!\n");
}

int Run(int argc, char **argv) {

	std::string program_name(argv[0]);

	auto parameters = raptor::createParamsRaptorIndex();

	if (raptor::ProcessArgsRaptorIndex(argc, argv, parameters)) {
		return 1;
	}

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters->verbose_level);
	if (parameters->verbose_level == 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
	} else if (parameters->verbose_level > 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
	}
	fflush(stdout);

	RunRaptorIndex(parameters);

	return 0;
}

int main(int argc, char **argv) {
  	return Run(argc, argv);
}

#endif
