/*
 * args.cc
 *
 *  Created on: Jun 04, 2019
 *      Author: isovic
 */

#include <args_raptor_fetch.h>
#include <args_raptor.h>
#include <lib/argparser.h>
#include <utility/stringutil.h>
#include <version.h>
#include <string>
#include <thread>
#include <iostream>
#include <utility/files.hpp>

namespace raptor {

int ProcessArgsRaptorFetch(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorFetch> parameters) {
    bool help = false;
    bool version = false;

    ArgumentParser argparser;

    argparser.AddArgument(&parameters->in_path, VALUE_TYPE_STRING, "i", "input", "",
                          "Input Raptor DB file.", 0,
                          "Input/Output options");

    argparser.AddArgument(&version, VALUE_TYPE_BOOL, "", "version", "0", "Output the version and exit.", 0,
                          "Other options");

    argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0,
                          "Other options");

    argparser.ProcessArguments(argc, argv);

    // Store the command line in parameters.
    std::stringstream ss_command_line;
    for (int i = 0; i < argc; i++) {
        if (i > 0) ss_command_line << " ";
        ss_command_line << argv[i];
    }
    parameters->command_line = ss_command_line.str();

    // Check if help was triggered.
    if (argparser.GetArgumentByLongName("help")->is_set == true) {
        std::stringstream ss;
        ss << SOFTWARE_NAME
           << " - Data fetching tool for Raptor - A very accurate and sensitive long-read, high error-rate sequence mapper\n"
           << SOFTWARE_NAME << " ";
        ss << "Version: " << RAPTOR_VERSION_STRING << "\n";
        ss << "Build date: " << std::string(RELEASE_DATE).c_str() << "\n";
        ss << "\n";
        ss << LICENCE_INFORMATION << "\n";
        ss << "\n";
        ss << "Usage:\n";
        ss << "\traptor-reshape [options] -r <reference_file>";
        ss << "\n";

        fprintf(stderr, "%s\n", ss.str().c_str());
        fprintf(stderr, "%s\n", argparser.VerboseUsage().c_str());
        exit(0);
    }

    if (version) {
        std::stringstream ss;
        ss << RAPTOR_VERSION_STRING << "\n";
        fprintf(stdout, "%s", ss.str().c_str());
        exit(0);
    }

    // Sanity check for the input path.
    if (!raptor::FileExists(parameters->in_path.c_str())) {
        fprintf(stderr, "Input file does not exist: '%s'\n\n", parameters->in_path.c_str());
        VerboseShortHelpAndExit(argc, argv);
    }

    // Write this out for every debug verbose level.
    if (parameters->verbose_level > 5) {
        fprintf(stderr, "%s\n", argparser.VerboseArguments().c_str());
    }

    return 0;
}

} /* namespace raptor */
