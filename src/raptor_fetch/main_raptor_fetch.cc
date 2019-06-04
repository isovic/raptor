#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
#include <fstream>
#include <args_raptor.h>
#include <limits>
#include <sequences/sequence_file.h>
#include <index/random_access_sequence_file.h>
#include <sequences/sequence_serializer.h>
#include <utility/memtime.h>
#include <raptor_fetch/params_raptor_fetch.h>

void RunTool(std::shared_ptr<raptor::ParamsRaptorFetch> parameters) {
	LOG_ALL("Raptor-fetch starts processing.\n");

	auto seq_file = mindex::createRandomAccessSequenceFile(parameters->in_path, 50);

	seq_file->FetchSequence(0);
	seq_file->FetchSequence(15000);
	seq_file->FetchSequence(3);
	seq_file->FetchSequence(12000);
	seq_file->FetchSequence("m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/95/0_11468");

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
