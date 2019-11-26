/*
 * args_raptor_sova.cc
 *
 *  Created on: Nov 12, 2019
 *      Author: isovic
 */

#include <raptor_sova/args_raptor_sova.h>
#include <raptor_sova/params_raptor_sova.h>
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
#include <index/index_types.h>
#include <sequences/sequence_file_utils.h>

namespace raptor {
namespace sova {

void VerboseShortHelpRaptorSova(int argc, char **argv) {
    fprintf(stderr, "For detailed help, please run with -h option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example usage:\n");
    fprintf(stderr, "  ./raptor-sova -r reads.fa -q reads.fa -o out.paf         # Overlapping.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s\n", LICENCE_INFORMATION.c_str());
    fprintf(stderr, "Version: %d.%d.%d-%s\n", RAPTOR_VERSION_MAJOR, RAPTOR_VERSION_MINOR,
            RAPTOR_VERSION_PATCH, RAPTOR_VERSION_COMMIT.c_str());
    fprintf(stderr, "Build date: %s\n", std::string(RELEASE_DATE).c_str());
    fprintf(stderr, "\n");
}

void VerboseShortHelpRaptorSovaAndExit(int argc, char **argv, int ret_val = 1) {
    VerboseShortHelpRaptorSova(argc, argv);
    exit(ret_val);
}

int ProcessArgsRaptorSova(int argc, char **argv, std::shared_ptr<raptor::sova::ParamsRaptorSova> parameters) {
    bool help = false;
    bool version = false;

    ArgumentParser argparser;

    std::string infmt, ref_fmt, graph_fmt, outfmt;
    std::string index_type_string;
    std::string index_region_string;
    std::string batch_size_str("0");

    // Define the composite options which can be expanded internally.
    std::string composite_str_ovl_hifi("--overlap-skip-self --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 --one-hit-per-target -k 25 -w 10 --end-bonus 200 --diff --flank-ext-len 100 --no-sezs");
    argparser.AddCompositeArgument("ovl", composite_str_ovl_hifi);

    argparser.AddArgument(&parameters->ref_paths, VALUE_TYPE_STRING_LIST, "r", "ref", "",
                          "Path to the reference sequence (fastq or fasta). Can be specified multiple times.", 0,
                          "Input/Output options");
    argparser.AddArgument(&parameters->query_paths, VALUE_TYPE_STRING_LIST, "q", "query", "",
                          "Path to the reads file. Can be specified multiple times.", 0, "Input/Output options");
    argparser.AddArgument(&parameters->out_path, VALUE_TYPE_STRING, "o", "out", "",
                          "Path to the output file that will be generated.", 0,
                          "Input/Output options");
    argparser.AddArgument(&ref_fmt, VALUE_TYPE_STRING, "", "ref-fmt", "auto",
                          "Format of the reference sequence file. Options are:"
                          "\n auto  - Determines the format automatically from file extension."
                          "\n fastq - Loads FASTQ or FASTA files.\n fasta - Loads FASTQ or FASTA files."
                          "\n gfa   - Graphical Fragment Assembly format (any)."
                          "\n gfa1  - GFA-1 format.\n gfa2  - GFA-2 format."
                          "\n sam   - Sequence Alignment/Mapping format."
#ifdef RAPTOR_COMPILED_WITH_PBBAM
                          "\n bam   - Binary Sequence Alignment/Mapping format."
#endif
                          "\n fofn  - File Of File Names."
                          "\n rdb   - RaptorDB format."
                          "\n ",
                          0, "Input/Output options");
    argparser.AddArgument(&infmt, VALUE_TYPE_STRING, "", "in-fmt", "auto",
                          "Format in which to input reads. Options are:"
                          "\n auto  - Determines the format automatically from file extension."
                          "\n fastq - Loads FASTQ or FASTA files."
                          "\n fasta - Loads FASTQ or FASTA files."
                          "\n gfa   - Graphical Fragment Assembly format."
                          "\n gfa1  - GFA-1 format."
                          "\n gfa2  - GFA-2 format"
                          "\n sam   - Sequence Alignment/Mapping format."
#ifdef RAPTOR_COMPILED_WITH_PBBAM
                          "\n bam   - Binary Sequence Alignment/Mapping format."
                          "\n xml   - PacBio Dataset format."
#endif
                          "\n fofn  - File Of File Names."
                          "\n rdb   - RaptorDB format."
                          "\n ",
                          0, "Input/Output options");
    argparser.AddArgument(&outfmt, VALUE_TYPE_STRING, "", "out-fmt", "paf",
                          "Format in which to output results. Options are:"
                          "\n sam  - Standard SAM output."
#ifdef RAPTOR_COMPILED_WITH_PBBAM
                          "\n bam   - Binary Sequence Alignment/Mapping format."
#endif
                          "\n paf  - PAF format useful for overlaps."
                          "\n gfa2 - GFA2 format."
                          "\n mhap - MHAP format useful for overlaps."
                          "\n m4   - BLASR M4 format."
                          "\n ",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->strict_format, VALUE_TYPE_BOOL, "", "strict-fmt", "0",
                          "If specified, all non-mandatory fields will be omitted from the output. Useful for tests.",
                          0, "Input/Output options");
    argparser.AddArgument(&batch_size_str, VALUE_TYPE_STRING, "B", "batch", "500",
                          "Queries will be loaded in batches of the size specified in megabytes. If there is a "
                          "trailing 'x', then the batch size is in terms of coverage of reference genome. "
                          "Value <= 0 loads the entire file.",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->rdb_block_ref, VALUE_TYPE_INT64, "", "rdb-block-ref", "-1",
                          "If the references are in the RaptorDB format, this specifies the block to index. "
                          "Value < 0 loads the entire file.",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->rdb_block_query, VALUE_TYPE_INT64, "", "rdb-block-query", "-1",
                          "If the queries are in the RaptorDB format, this specifies a particular block to process. "
                          "Value < 0 processes each block as a batch, sequentially.",
                          0, "Input/Output options");

    argparser.AddArgument(&(index_type_string), VALUE_TYPE_STRING, "", "index-type", "minimizer",
                          "Type of index to use:"
                          "\n  minimizer"
                          "\n  dense"
                          "\n ",
                          0, "Index options");
    argparser.AddArgument(&(index_region_string), VALUE_TYPE_STRING, "", "region", "",
                          "Indexes only the specified region. Region format: chr:start-end. Start "
                          "is 0-based, and end is not inclusive. If end is <= 0, the entire suffix "
                          "of the reference will be indexed.",
                          0, "Index options");
    argparser.AddArgument(
        &(parameters->index_params->index_only_fwd_strand), VALUE_TYPE_BOOL, "", "fwd", "0",
        "If specified, only the fwd strand of the reference will be indexed.", 0, "Index options");
    argparser.AddArgument(&(parameters->index_params->k), VALUE_TYPE_INT32, "k", "seed-len", "15",
                          "Length of the seeds used for hashing and lookup.", 0, "Index options");
    argparser.AddArgument(&(parameters->index_params->w), VALUE_TYPE_INT32, "w", "minimizer-window",
                          "5",
                          "Length of the window to select a minimizer from if index type is 'minimizer'. Step size if index type is 'dense'.",
                          0, "Index options");
    // argparser.AddArgument(&(parameters->index_params->base_spacing), VALUE_TYPE_INT32, "", "space",
    //                       "1",
    //                       "Spaced seeds. Every <val> base will be used to form a seed.",
    //                       0, "Index options");
    argparser.AddArgument(&(parameters->index_params->freq_percentil), VALUE_TYPE_DOUBLE, "",
                          "freq-percentile", "0.0002",
                          "Filer this fraction of most frequent seeds in the lookup process.", 0,
                          "Index options");
    argparser.AddArgument(&(parameters->index_params->min_occ_cutoff), VALUE_TYPE_INT32, "",
                          "min-occ-cutoff", "0",
                          "If the freq-percentil value is less than this count, this will be used "
                          "instead for seed count filtering.",
                          0, "Index options");
    argparser.AddArgument(
        &(parameters->index_params->homopolymer_suppression), VALUE_TYPE_BOOL, "", "hp-suppress",
        "0", "If specified, homopolymer runs will be suppressed into a single seed point.", 0,
        "Index options");
    argparser.AddArgument(&(parameters->index_params->max_homopolymer_len), VALUE_TYPE_INT32, "",
                          "max-hp-len", "5", "Compress only homopolymer runs below this length.", 0,
                          "Index options");
    argparser.AddArgument(
        &parameters->keep_lowercase, VALUE_TYPE_BOOL, "", "keep-lowercase", "0",
        "If this parameter is not specified, lowercase bases will be converted to uppercase.", 0,
        "Index options");
    argparser.AddArgument(&parameters->index_params->min_tlen, VALUE_TYPE_INT64, "", "min-tlen",
                          "0",
                          "Sequences shorter than this will not be indexed.",
                          0, "Index options");

    // Mapping options.
    argparser.AddArgument(&parameters->mapper_params->min_qlen, VALUE_TYPE_INT32, "", "min-qlen",
                          "50",
                          "If a query is shorter than this, it will be marked as unmapped. This "
                          "value can be lowered if the reads are known to be accurate.",
                          0, "Mapping options");

    argparser.AddArgument(&parameters->mapper_params->seed_join_dist, VALUE_TYPE_INT32, "",
                          "max-seed-dist", "5000",
                          "Maximum distance between two seeds to combine in an anchor after LIS.",
                          0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->min_num_seeds, VALUE_TYPE_INT32, "",
                          "min-num-seeds", "3",
                          "Minimum number of seeds in an anchor. Less important as it is trumped "
                          "by min-cov-bases which is more specific.",
                          0, "Mapping options");
    argparser.AddArgument(
        &parameters->mapper_params->min_cov_bases, VALUE_TYPE_INT32, "", "min-cov-bases", "30",
        "Minimum number of covered bases in an anchor to keep it.", 0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->chain_min_span, VALUE_TYPE_INT32, "", "min-chain-span",
                          "1000",
                          "Minimum chain span to ratain it.",
                          0, "Mapping options");

    argparser.AddArgument(&parameters->mapper_params->chain_bandwidth, VALUE_TYPE_INT32, "", "chain-bw",
                          "100",
                          "Diagonal bandwidth to merge seeds into chains.",
                          0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->align_bandwidth, VALUE_TYPE_DOUBLE, "", "aln-bw",
                          "0.01",
                          "Bandwidth for alignment, fraction of the query span.",
                          0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->align_max_diff, VALUE_TYPE_DOUBLE, "", "aln-maxd",
                          "0.03",
                          "Expected maximum diff rate between sequences.",
                          0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->min_identity, VALUE_TYPE_DOUBLE, "",
                          "min-idt", "98",
                          "Minimum percent alignment identity allowed to report the alignment.",
                          0, "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->min_map_len, VALUE_TYPE_INT64, "", "min-map-len", "1000",
                          "Output only alignments/mappings/overlaps above this length.", 0,
                          "Mapping options");

    argparser.AddArgument(
        &parameters->mapper_params->overlap_skip_self_hits, VALUE_TYPE_BOOL, "", "overlap-skip-self", "0",
        "A mapping is discarded if the query and the target headers are the same. "
        "Applicable only if list of reference sequences is the same as list of query sequences (order is important).",
        0, "Overlapping options");

    argparser.AddArgument(
            &parameters->mapper_params->overlap_single_arc, VALUE_TYPE_BOOL, "", "overlap-single-arc", "0",
        "Output overlaps only in one strand of a query-target pair. "
        "Applicable only if list of reference sequences is the same as list of query sequences (order is important).",
        0, "Overlapping options");

    // Filtering options.
    argparser.AddArgument(&parameters->bestn, VALUE_TYPE_INT64, "", "bestn", "0",
                          "Output best N alignments/mappings/overlaps. If <= 0 all mappings within "
                          "bestn-threshold from best score will be output.",
                          0, "Filtering options");
    argparser.AddArgument(&parameters->one_hit_per_target, VALUE_TYPE_BOOL, "",
                          "one-hit-per-target", "0",
                          "Only one query-target pair will be reported. Useful for overlapping.",
                          0, "Filtering options");



    argparser.AddArgument(&parameters->num_threads, VALUE_TYPE_INT64, "t", "threads", "-1",
                          "Number of threads to use. If '-1', number of threads will be equal to "
                          "min(24, num_cores/2).",
                          0, "Other options");
    argparser.AddArgument(
        &parameters->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5",
        "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0,
        "Other options");
    argparser.AddArgument(&parameters->start_read, VALUE_TYPE_INT64, "s", "start", "0",
                          "Ordinal number of the read from which to start processing data.", 0,
                          "Other options");
    argparser.AddArgument(
        &parameters->num_reads_to_process, VALUE_TYPE_INT64, "n", "num-reads", "-1",
        "Number of reads to process per batch. Value of '-1' processes all reads.", 0,
        "Other options");

    argparser.AddArgument(&parameters->debug_qid, VALUE_TYPE_INT64, "y", "debug-read", "-1",
                          "ID of the read to give the detailed verbose output.", 0,
                          "Debug options");
    argparser.AddArgument(&parameters->debug_qname, VALUE_TYPE_STRING, "Y", "debug-qname", "",
                          "QNAME of the read to give the detailed verbose output. Has precedence "
                          "over -y. Use quotes to specify.",
                          0, "Debug options");

    argparser.AddArgument(&parameters->composite, VALUE_TYPE_COMPOSITE, "x", "composite", "",
        "Pre-set parameters for different use cases. Valid options are:\n"
        " \"ovl\"   - Equivalent to: '" + composite_str_ovl_hifi + "'\n"
        );

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
        VerboseShortHelpRaptorSova(argc, argv);
        fprintf(stderr, "%s\n", argparser.VerboseUsage().c_str());
        exit(0);
    }

    if (version) {
        std::stringstream ss;
        ss << RAPTOR_VERSION_STRING << "\n";
        fprintf(stdout, "%s", ss.str().c_str());
        exit(0);
    }

    // In case debug_qname was set, and it was specified using quote signs on the command line,
    // remove the quotes.
    if (parameters->debug_qname.size() > 2 && ((parameters->debug_qname.front())) == '"' &&
        ((parameters->debug_qname.back())) == '"') {
        parameters->debug_qname =
            parameters->debug_qname.substr(1, (parameters->debug_qname.size() - 2));
    }

    // Sanity check for the reference path.
    if (parameters->ref_paths.size() == 0) {
        fprintf(stderr, "Please specify the path to the reference file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    // Sanity check for the reads path.
    if (parameters->query_paths.size() == 0 && parameters->calc_only_index == false) {
        fprintf(stderr, "Please specify the path to the reads file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    /////////////////////////////////////////////////
    /// Parsing and validating the input formats. ///
    /////////////////////////////////////////////////
    // Parse the input format.
    parameters->infmt = mindex::SequenceFormatFromString(infmt);
    if (parameters->infmt == mindex::SequenceFormat::Unknown) {
        fprintf(stderr, "Unknown input format '%s'!\n\n", infmt.c_str());
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }
    // Collect all FOFN files.
    parameters->query_paths = ExpandPathList(parameters->infmt, infmt, parameters->query_paths);
    // Validate the input files and formats.
    bool validate_rv1 = ValidateInputFiles(parameters->infmt, parameters->query_paths);
    if (validate_rv1 == false) {
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }
    // In case the input was RaptorDB, modify the infmt for future use in the index factory.
    if (IsInputFormatRaptorDB(parameters->infmt, parameters->query_paths)) {
        parameters->infmt = mindex::SequenceFormat::RaptorDB;
    }
    #ifdef RAPTOR_COMPILED_WITH_PBBAM
        if (IsInputFormatXML(parameters->infmt, parameters->query_paths)) {
            parameters->infmt = mindex::SequenceFormat::XML;
        }
    #endif
    /////////////////////////////////////////////////

    /////////////////////////////////////////////////////
    /// Parsing and validating the reference formats. ///
    /////////////////////////////////////////////////////
    // Parse the reference format.
    parameters->ref_fmt = mindex::SequenceFormatFromString(ref_fmt);
    if (parameters->ref_fmt == mindex::SequenceFormat::Unknown) {
        fprintf(stderr, "Unknown reference format '%s'!\n\n", ref_fmt.c_str());
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }
    // Collect all FOFN files.
    parameters->ref_paths = ExpandPathList(parameters->ref_fmt, infmt, parameters->ref_paths);
    // Validate the input files and formats.
    bool validate_rv2 = ValidateInputFiles(parameters->ref_fmt, parameters->ref_paths);
    if (validate_rv2 == false) {
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }
    // In case the input was RaptorDB, modify the infmt for future use in the index factory.
    if (IsInputFormatRaptorDB(parameters->ref_fmt, parameters->ref_paths)) {
        parameters->ref_fmt = mindex::SequenceFormat::RaptorDB;
    }
    #ifdef RAPTOR_COMPILED_WITH_PBBAM
        if (IsInputFormatXML(parameters->ref_fmt, parameters->ref_paths)) {
            parameters->ref_fmt = mindex::SequenceFormat::XML;
        }
    #endif
    /////////////////////////////////////////////////////

    // Parse the output format.
    parameters->outfmt = raptor::OutputFormatFromString(outfmt);
    if (parameters->outfmt == OutputFormat::Unknown) {
        fprintf(stderr, "Unknown output format '%s'!\n\n", outfmt.c_str());
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    if (parameters->num_threads < 0) {
        parameters->num_threads = std::min(24, ((int)std::thread::hardware_concurrency()) / 2);
    }

    if (parameters->index_params->w <= 0) {
        fprintf(stderr, "Minimizer window length cannot be <= 0!\n");
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    if (parameters->index_params->freq_percentil < 0.0 ||
        parameters->index_params->freq_percentil > 1.0) {
        fprintf(stderr, "Frequency percentil should be in range [0.0, 1.0].\n");
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    // Write this out for every debug verbose level.
    if (parameters->verbose_level > 5) {
        fprintf(stderr, "%s\n", argparser.VerboseArguments().c_str());
    }

    parameters->batch_size = 0.0;
    if (batch_size_str.size() > 0) {
        if (batch_size_str.back() == 'x') {
            parameters->batch_size = std::stod(batch_size_str.substr(0, batch_size_str.size() - 1));
            parameters->batch_type = mindex::BatchLoadType::Coverage;
        } else {
            parameters->batch_size = std::stod(batch_size_str);
        }
    }

    parameters->ref_and_reads_path_same = false;
    if (parameters->ref_paths == parameters->query_paths) {
        parameters->ref_and_reads_path_same = true;
    }
    parameters->mapper_params->ref_and_reads_path_same = parameters->ref_and_reads_path_same;

    parameters->index_params->is_region_specified = false;
    if (index_region_string.size() > 0) {
        parameters->index_params->is_region_specified = true;

        size_t rname_end = index_region_string.find_last_of(":");
        if (rname_end == std::string::npos) {
            parameters->index_params->region_rname = index_region_string;
            parameters->index_params->region_rstart = 0;
            parameters->index_params->region_rend = -1;
        } else {
            parameters->index_params->region_rname = index_region_string.substr(0, rname_end);
            size_t dash_pos = index_region_string.find_first_of("-", (rname_end + 1));
            if (dash_pos == std::string::npos) {
                fprintf(stderr,
                        "Region format is incorrect. Two coordinates need to be provided, "
                        "separated by '-', in the format: chr:start-end.\n");
                VerboseShortHelpRaptorSovaAndExit(argc, argv);
            }
            parameters->index_params->region_rstart =
                atoi(index_region_string.substr((rname_end + 1), dash_pos).c_str());
            parameters->index_params->region_rend =
                atoi(index_region_string.substr((dash_pos + 1)).c_str());
        }
    }

    // Decode the index type.
    if (index_type_string == "minimizer") {
        parameters->index_type = mindex::IndexType::Minimizer;
    } else if (index_type_string == "dense") {
        parameters->index_type = mindex::IndexType::Dense;
    } else {
        parameters->index_type = mindex::IndexType::Undefined;
        fprintf(stderr, "Unknown index type.\n");
        VerboseShortHelpRaptorSovaAndExit(argc, argv);
    }

    parameters->mapper_params->verbose_level = parameters->verbose_level;
    parameters->mapper_params->debug_qid = parameters->debug_qid;
    parameters->mapper_params->debug_qname = parameters->debug_qname;

#ifdef RAPTOR_TESTING_MODE
    if (parameters->verbose_level >= 5) {
        fprintf (stderr, "Compiled with RAPTOR_TESTING_MODE.\n");
    }
#endif

    return 0;
}

} /* namespace sova */
} /* namespace raptor */
