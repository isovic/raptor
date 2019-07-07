/*
 * args.cc
 *
 *  Created on: Mar 07, 2018
 *      Author: isovic
 */

#include <raptor_reshape/args_raptor_reshape.h>
#include <raptor/args_raptor.h>
#include <params/params_raptor.h>
#include <lib/argparser.h>
#include <utility/stringutil.h>
#include <version.h>
#include <string>
#include <thread>
#include <aligner/global_margins.h>
#include <aligner/pairwise_penalties.h>
#include <iostream>
#include <utility/files.hpp>
#include <utility/fofn.h>

namespace raptor {

int ProcessArgsRaptorReshape(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorReshape> parameters) {
    bool help = false;
    bool version = false;

    ArgumentParser argparser;
    std::string in_fmt, out_fmt;

    argparser.AddArgument(&parameters->in_paths, VALUE_TYPE_STRING_LIST, "i", "in", "",
                          "Path to the sequences (fastq or fasta). Can be specified multiple times.", 0,
                          "Input/Output options");
    argparser.AddArgument(&parameters->out_prefix, VALUE_TYPE_STRING, "o", "out-prefix", "",
                          "Prefix of the output files that will be generated.", 0,
                          "Input/Output options");
    argparser.AddArgument(&in_fmt, VALUE_TYPE_STRING, "", "in-fmt", "auto",
                          "Format of the input sequence file. Options are:"
                          "\n auto  - Determines the format automatically from file extension."
                          "\n fasta - Loads FASTQ or FASTA files."
                          "\n fastq - Loads FASTQ or FASTA files."
                          "\n gfa   - Graphical Fragment Assembly format."
                          "\n gfa1  - GFA-1 format."
                          "\n gfa2  - GFA-2 format."
                          "\n sam   - Sequence Alignment/Mapping format."
#ifdef RAPTOR_COMPILED_WITH_PBBAM
                          "\n bam   - Binary Sequence Alignment/Mapping format."
                          "\n xml   - PacBio Dataset format."
#endif
                          "\n fofn  - File Of File Names. Format is determined automatically.",
                          0, "Input/Output options");

    argparser.AddArgument(&out_fmt, VALUE_TYPE_STRING, "", "out-fmt", "fasta",
                          "Format of the output file(s), unless '--symlink' is used. Options are:\n fastq, fasta.",
                        //   "Format of the output file(s), unless '--symlink' is used. Options are:\n fastq, fasta, gfa1, gfa2, sam.",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->in_batch_size, VALUE_TYPE_DOUBLE, "", "in-batch", "400",
                          "Queries will be loaded in batches of the size specified in megabytes. "
                          "Value <= 0 loads the entire file.",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->block_size, VALUE_TYPE_DOUBLE, "", "block-size", "400",
                          "Partition sequences into blocks. "
                          "Value <= 0 loads the entire file.",
                          0, "Input/Output options");
    argparser.AddArgument(
        &parameters->split_blocks, VALUE_TYPE_BOOL, "", "split-blocks", "0",
        "If true, each block into a separate file.", 0,
        "Input/Output options");
    argparser.AddArgument(
        &parameters->rename_seqs, VALUE_TYPE_BOOL, "", "rename", "0",
        "Sequences will be renamed in order of appearance.", 0,
        "Input/Output options");
    argparser.AddArgument(
        &parameters->keep_lowercase, VALUE_TYPE_BOOL, "", "keep-lowercase", "0",
        "If this parameter is not specified, lowercase bases will be converted to uppercase.", 0,
        "Other options");
    argparser.AddArgument(
        &parameters->symlink_files, VALUE_TYPE_BOOL, "", "symlink", "0",
        "Sequences will not be written to block files, but only indexed, while still pointing to their original files.", 0,
        "Other options");
    argparser.AddArgument(
        &parameters->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5",
        "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0,
        "Other options");
    // argparser.AddArgument(
    //     &parameters->filter_pb_max, VALUE_TYPE_BOOL, "", "filt-pb-max", "0",
    //     "Keep only maximum length subread per ZMW.", 0,
    //     "Other options");

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
           << " - Data reshaping tool for Raptor - A very accurate and sensitive long-read, high error-rate sequence mapper\n"
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

    // Sanity check for the reference path.
    if (parameters->in_paths.size() == 0) {
        fprintf(stderr, "Please specify the path to the reference file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpAndExit(argc, argv);
    }

    /////////////////////////////////////////////////
    /// Parsing and validating the input formats. ///
    /////////////////////////////////////////////////
    // First check that the specified format is valid.
    parameters->in_fmt = mindex::SequenceFormatFromString(in_fmt);
    if (parameters->in_fmt == mindex::SequenceFormat::Unknown) {
        fprintf(stderr, "Unsupported input format '%s'!\n\n", in_fmt.c_str());
        VerboseShortHelpAndExit(argc, argv);
    }
    // Parse the input format.
    parameters->in_fmt = mindex::SequenceFormatFromString(in_fmt);
    if (parameters->in_fmt == mindex::SequenceFormat::Unknown) {
        fprintf(stderr, "Unknown input format '%s'!\n\n", in_fmt.c_str());
        VerboseShortHelpAndExit(argc, argv);
    }
    // Collect all FOFN files.
    parameters->in_paths = ExpandPathList(parameters->in_fmt, in_fmt, parameters->in_paths);
    // Validate the input files and formats.
    ValidateInputFiles(argc, argv, parameters->in_fmt, parameters->in_paths);
    // In case the input was RaptorDB, modify the infmt for future use in the index factory.
    if (IsInputFormatRaptorDB(parameters->in_fmt, parameters->in_paths)) {
        parameters->in_fmt = mindex::SequenceFormat::RaptorDB;
    }
    parameters->in_fmt = (parameters->in_fmt == mindex::SequenceFormat::FOFN) ? mindex::SequenceFormat::Auto : parameters->in_fmt;
    /////////////////////////////////////////////////

    // Parse the output format.
    parameters->out_fmt = mindex::SequenceFormatFromString(out_fmt);
    if (parameters->out_fmt == mindex::SequenceFormat::Unknown || parameters->out_fmt == mindex::SequenceFormat::Auto || parameters->out_fmt == mindex::SequenceFormat::GFA) {
        fprintf(stderr, "Unsupported output format '%s'!\n\n", out_fmt.c_str());
        VerboseShortHelpAndExit(argc, argv);
    }

    if (parameters->symlink_files) {
        if (parameters->rename_seqs) {
            fprintf(stderr, "Option '--rename' cannot be used in combination with '--symlink'.\n\n");
            VerboseShortHelpAndExit(argc, argv);
        }
        if (parameters->keep_lowercase) {
            fprintf(stderr, "Option '--keep-lowercase' cannot be used in combination with '--symlink'.\n\n");
            VerboseShortHelpAndExit(argc, argv);
        }
        if (parameters->split_blocks) {
            fprintf(stderr, "Option '--split-blocks' cannot be used in combination with '--symlink'.\n\n");
            VerboseShortHelpAndExit(argc, argv);
        }
    }

    // Write this out for every debug verbose level.
    if (parameters->verbose_level > 5) {
        fprintf(stderr, "%s\n", argparser.VerboseArguments().c_str());
    }

    return 0;
}

} /* namespace raptor */
