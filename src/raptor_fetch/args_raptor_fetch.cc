/*
 * args.cc
 *
 *  Created on: Jun 04, 2019
 *      Author: isovic
 */

#include <raptor_fetch/args_raptor_fetch.h>
#include <raptor/args_raptor.h>
#include <lib/argparser.h>
#include <utility/stringutil.h>
#include <version.h>
#include <string>
#include <thread>
#include <iostream>
#include <utility/files.hpp>

namespace raptor {

int ProcessArgsRaptorFetch(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorFetch> params) {
    bool help = false;
    bool version = false;

    ArgumentParser argparser;

    argparser.AddArgument(&params->rdb_path, VALUE_TYPE_STRING, "r", "rdb", "",
                          "Path to the RaptorDB.", 0,
                          "Input/Output options");
    argparser.AddArgument(&params->in_paths, VALUE_TYPE_STRING_LIST, "i", "input", "",
                          "One or more input overlap paths.", 0,
                          "Input/Output options");
    argparser.AddArgument(&params->out_prefix, VALUE_TYPE_STRING, "o", "out-prefix", "",
                          "Prefix of the output files to generate.", 0,
                          "Input/Output options");
    argparser.AddArgument(&params->job_str, VALUE_TYPE_STRING, "", "job", "erc",
                          "Options:"
                          "\n fetch - The input is a list of qnames, one per line. Fetches the sequences from the rdb."
                          "\n clip  - The input is a BED-3 file of sequences and regions to extract."
                          "\n erc   - The input is an M4 file of overlaps. Output is a pile of sequences for Falcon consensus."
                          , 0,
                          "Input/Output options");
    argparser.AddArgument(&params->min_cov, VALUE_TYPE_INT32, "", "min-cov", "0",
                          "Minimum coverage of a sequence to retain it.",
                          0, "Filtering");
    argparser.AddArgument(&params->min_len, VALUE_TYPE_INT64, "", "min-len", "0",
                          "Minimum length of a sequence or a clipped span to retain it.",
                          0, "Filtering");
    argparser.AddArgument(&params->min_score, VALUE_TYPE_DOUBLE, "", "min-score", "0",
                          "Overlaps with score below this value will not be considered for processing.",
                          0, "Filtering");
    argparser.AddArgument(&params->min_idt, VALUE_TYPE_DOUBLE, "", "min-idt", "0.0",
                          "Overlaps with the identity value below this will not be considered for processing.",
                          0, "Filtering");
    argparser.AddArgument(
        &params->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5",
        "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0,
        "Other options");

    argparser.AddArgument(&version, VALUE_TYPE_BOOL, "", "version", "0", "Output the version and exit.", 0,
                          "Other options");

    argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0,
                          "Other options");

    argparser.ProcessArguments(argc, argv);

    // Store the command line in params.
    std::stringstream ss_command_line;
    for (int i = 0; i < argc; i++) {
        if (i > 0) ss_command_line << " ";
        ss_command_line << argv[i];
    }
    params->command_line = ss_command_line.str();

    // Check if help was triggered.
    if (argparser.GetArgumentByLongName("help")->is_set == true) {
        VerboseShortHelp(argc, argv);
        fprintf(stderr, "%s\n", argparser.VerboseUsage().c_str());
        exit(0);
    }

    if (version) {
        std::stringstream ss;
        ss << RAPTOR_VERSION_STRING << "\n";
        fprintf(stdout, "%s", ss.str().c_str());
        exit(0);
    }

    // Sanity check for the RaptorDB.
    if (argparser.GetArgumentByLongName("rdb")->is_set == false) {
        fprintf(stderr, "Please specify the path to the RaptorDB file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpAndExit(argc, argv);
    }
    if (!raptor::FileExists(params->rdb_path.c_str())) {
        fprintf(stderr, "RaptorDB does not exist: '%s'\n\n", params->rdb_path.c_str());
        VerboseShortHelpAndExit(argc, argv);
    }

    // Sanity check for the input paths.
    if (params->in_paths.size() == 0) {
        fprintf(stderr, "Please specify the path to the input file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpAndExit(argc, argv);
    }
    for (auto& in_path: params->in_paths) {
        if (!raptor::FileExists(in_path.c_str())) {
            fprintf(stderr, "Reference does not exist: '%s'\n\n", in_path.c_str());
            VerboseShortHelpAndExit(argc, argv);
        }
    }

    // Write this out for every debug verbose level.
    if (params->verbose_level > 5) {
        fprintf(stderr, "%s\n", argparser.VerboseArguments().c_str());
    }

    return 0;
}

} /* namespace raptor */
