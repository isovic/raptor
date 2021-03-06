/*
 * args_raptor.cc
 *
 *  Created on: Mar 07, 2018
 *      Author: isovic
 */

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
#include <index/index_types.h>
#include <sequences/sequence_file_utils.h>

namespace raptor {

void VerboseShortHelpRaptor(int argc, char **argv) {
    fprintf(stderr, "For detailed help, please run with -h option.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example usage:\n");
    fprintf(stderr, "  ./raptor -r escherichia_coli.fa -q reads.fastq -o out.paf         # Mapping.\n");
    fprintf(stderr, "  ./raptor -r escherichia_coli.fa -q reads.fastq --align -o out.paf # Align.\n");
    fprintf(stderr, "  ./raptor -g graph.gfa2 -r ref.fa -q reads.fasta -o out.paf        # Graph mapping.\n");
    fprintf(stderr, "  ./raptor-reshape -i reads.fasta --block-size 400 -o reads         # Build a RaptorDB.\n");
    fprintf(stderr, "  ./raptor -x overlap -r reads.rdb -q reads.rdb -o overlaps.paf     # Map/overlap to the RDB.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s\n", LICENCE_INFORMATION.c_str());
    fprintf(stderr, "Version: %d.%d.%d-%s\n", RAPTOR_VERSION_MAJOR, RAPTOR_VERSION_MINOR,
            RAPTOR_VERSION_PATCH, RAPTOR_VERSION_COMMIT.c_str());
    fprintf(stderr, "Build date: %s\n", std::string(RELEASE_DATE).c_str());
    fprintf(stderr, "\n");
}

void VerboseShortHelpRaptorAndExit(int argc, char **argv, int ret_val) {
    VerboseShortHelpRaptor(argc, argv);
    exit(ret_val);
}

int ProcessArgsRaptor(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptor> parameters) {
    bool help = false;
    bool version = false;

    ArgumentParser argparser;

    std::string infmt, ref_fmt, graph_fmt, outfmt;
    std::string index_type_string;
    std::string index_region_string;
    std::string batch_size_str("0");
    bool no_relabel_secondary_supp = false;

    // Define the composite options which can be expanded internally.
    // std::string composite_str_overlap("--overlap-skip-self --overlap-single-arc -B 1000 --min-map-len 1000 --bestn 0 --bestn-threshold 1.0");
    std::string composite_str_ovl_raw("--overlap-skip-self --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 --sec-ratio 0.0 --one-hit-per-target -k 15 -w 10 --no-relabel");
    argparser.AddCompositeArgument("ovl-raw", composite_str_ovl_raw);

    std::string composite_str_ovl_hifi("--overlap-skip-self --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 --sec-ratio 0.0 --one-hit-per-target -k 25 -w 10 --end-bonus 200 --diff --flank-ext-len 100 --no-sezs --no-relabel");
    argparser.AddCompositeArgument("ovl-hifi", composite_str_ovl_hifi);

    std::string composite_str_ovl_miniasm("--overlap-skip-self --overlap-single-arc --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 --sec-ratio 0.0 --one-hit-per-target --no-relabel");
    argparser.AddCompositeArgument("ovl-miniasm", composite_str_ovl_miniasm);

    // std::string composite_str_ovl_ipa("--overlap-skip-self --min-map-len 1000 --bestn 0 --bestn-threshold 1.0 --one-hit-per-target "
    //                             "--match 1 --mismatch 1 --gap-open 1 --gap-ext 1 "
    //                             "--zdrop 400 --end-bonus 500 --no-sezs --flank-ext-len 500"); // --sezs
    // argparser.AddCompositeArgument("ipa", composite_str_ovl_ipa);

    argparser.AddArgument(&parameters->ref_paths, VALUE_TYPE_STRING_LIST, "r", "ref", "",
                          "Path to the reference sequence (fastq or fasta). Can be specified multiple times.", 0,
                          "Input/Output options");
    argparser.AddArgument(&parameters->query_paths, VALUE_TYPE_STRING_LIST, "q", "query", "",
                          "Path to the reads file. Can be specified multiple times.", 0, "Input/Output options");
    argparser.AddArgument(&parameters->graph_path, VALUE_TYPE_STRING, "g", "agraph", "",
                          "Path to the graph file. Graph is assumed to be asymmetric in the sense "
                          "of reverse complement, and revcmp edges will be added internally. This "
                          "is useful for mapping to assembly graphs.",
                          0, "Input/Output options");
    argparser.AddArgument(&parameters->graph_path, VALUE_TYPE_STRING, "G", "graph", "",
                          "Path to the graph file (GFA-1 or GFA-2, automatically determined from the header).", 0, "Input/Output options");
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
    argparser.AddArgument(&graph_fmt, VALUE_TYPE_STRING, "", "graph-fmt", "gfa",
                          "Format of the input graph. Options are:\n gfa - Graphical Fragment "
                          "Assembly format, version 1 or 2.",
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

    // argparser.AddArgument(&parameters->is_rna, VALUE_TYPE_BOOL, "", "rna", "0", "Query sequences
    // are RNAseq/IsoSeq. Impacts chaining process.", 0, "Input/Output options");

    // argparser.AddArgument(&parameters->rebuild_index, VALUE_TYPE_BOOL, "", "rebuild-index", "0",
    //                       "Always rebuild index even if it already exists in given path.", 0,
    //                       "Index options");
    // argparser.AddArgument(
    //     &parameters->auto_rebuild_index, VALUE_TYPE_BOOL, "", "auto-rebuild-index", "1",
    //     "Rebuild index only if an existing index is of an older version or corrupt.", 0,
    //     "Index options");
    // argparser.AddArgument(&parameters->index_on_the_fly, VALUE_TYPE_BOOL, "", "fly-index", "0",
    //                       "Index will be constructed on the fly, without storing it to disk. If it "
    //                       "already exists on disk, it will be loaded unless --rebuild-index is "
    //                       "specified.",
    //                       0, "Index options");
    // argparser.AddArgument(&parameters->calc_only_index, VALUE_TYPE_BOOL, "I", "index-only", "0",
    //                       "Build only the index from the given reference and exit. If not "
    //                       "specified, index will automatically be built if it does not exist, or "
    //                       "loaded from file otherwise.",
    //                       0, "Index options");
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

    // argparser.AddArgument(
    //     &(parameters->mapper_params->max_hits), VALUE_TYPE_INT32, "", "max-hits", "-1",
    //     "If a query has more hits than this, set it as unmapped. If -1, no threshold is set.", 0,
    //     "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->min_qlen, VALUE_TYPE_INT32, "", "min-qlen",
                          "50",
                          "If a query is shorter than this, it will be marked as unmapped. This "
                          "value can be lowered if the reads are known to be accurate.",
                          0, "Mapping options");

    argparser.AddArgument(
        &parameters->mapper_params->diag_margin, VALUE_TYPE_INT32, "", "max-diag-margin", "500",
        "Maximum diagonal shift between two seeds to not include them in the same LIS frame.", 0,
        "Mapping options");
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
    argparser.AddArgument(
        &parameters->mapper_params->min_dp_score, VALUE_TYPE_INT32, "", "min-dp-score", "60",
        "Minimum number of covered bases in an anchor to keep it.", 0, "Mapping options");
    argparser.AddArgument(
        &parameters->mapper_params->chain_max_skip, VALUE_TYPE_INT32, "",
        "chain-max-skip", "25",
        "Heuristic to reduce the number of predecessors to inspect during the chaining DP.", 0,
        "Mapping options");
    argparser.AddArgument(
        &parameters->mapper_params->chain_max_predecessors, VALUE_TYPE_INT32, "",
        "chain-max-pred", "500",
        "Heuristic to reduce the number of predecessors to inspect during the chaining DP.", 0,
        "Mapping options");
    argparser.AddArgument(
        &parameters->mapper_params->chain_max_dist, VALUE_TYPE_INT32, "", "chain-max-dist", "10000",
        "Maximum distance between two anchors to chain them. "
        "This parameter is also used as a predecessor lookup window in graph-based mapping.", 0, "Mapping options");
    argparser.AddArgument(
        &parameters->mapper_params->chain_max_bandwidth, VALUE_TYPE_INT32, "",
        "chain-max-bandwidth", "3000",
        "Maximum allowed gap between anchors to chain them. Use a larger value for larger SVs.", 0,
        "Mapping options");
    argparser.AddArgument(&parameters->mapper_params->graph_allowed_score_diff_frac, VALUE_TYPE_DOUBLE, "",
                          "graph-score-diff-frac", "0.95",
                          "Allowed fraction from the maximum backtrack edge, or leaf node, "
                          "to retain an alternate path.",
                          0, "Graph Mapping options");
    argparser.AddArgument(&parameters->mapper_params->graph_max_path_edges, VALUE_TYPE_INT32, "",
                          "graph-max-path-edges", "100",
                          "Maximum number of edges to follow during local graph construction. "
                          "Value <= 0 will follow all available edges.",
                          0, "Graph Mapping options");
    argparser.AddArgument(&parameters->mapper_params->score_anchors, VALUE_TYPE_BOOL, "",
                          "score-anchors", "0",
                          "If specified, anchors will be scored via alignment during linear mapping. This can produce "
                          "more accurate chaining.",
                          0, "Graph Mapping options");
    argparser.AddArgument(&parameters->mapper_params->graph_score_anchors, VALUE_TYPE_BOOL, "",
                          "graph-score-anchors", "0",
                          "If specified, anchors will be scored via alignment during graph-based mapping. This can produce "
                          "more accurate chaining.",
                          0, "Graph Mapping options");
    argparser.AddArgument(&parameters->mapper_params->no_graph_mapping, VALUE_TYPE_BOOL, "",
                          "no-gm", "0",
                          "If specified, graph mapping will not be applied. Instead, only linear mapping anchors will "
                          "be used. Useful for faster overlapping.",
                          0, "Graph Mapping options");

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

    argparser.AddArgument(
        &parameters->mapper_params->flank_ext_len, VALUE_TYPE_INT32, "", "flank-ext-len", "0",
        "If != 0, front and back of a path will be extended using DP. Value == 0 is a no-op. "
        "Value > 0 sets a limit to the maximum flank length for extension. Value < 0 will force extend "
        "flanks of any length (essentially - end-to-end alignment).",
        0, "Overlapping options");

    // argparser.AddArgument(&parameters->mapper_params->graph_chain_inversions, VALUE_TYPE_BOOL, "",
    //                       "chain-inversions", "0",
    //                       "Attempt to chain reverse complement strand of the read as well, to "
    //                       "produce higher scoring paths.",
    //                       0, "Mapping options");



    int32_t match = 5;
    int32_t mismatch = 4;
    std::string gap_open_penalty_str;
    std::string gap_ext_penalty_str;
    std::string aligner_str;
    int32_t bandwidth = -1;
    int32_t zbandwidth = -1;
    int32_t end_bonus = 0;
    bool no_stop_ext_on_zero_score = false;
    double min_secondary_to_primary_ratio = 0.80;
    int32_t allowed_suppl_overlap = 0;

    argparser.AddArgument(&parameters->do_align, VALUE_TYPE_BOOL, "", "align", "0",
                          "If selected, alignment will be produced for the mappings.", 0,
                          "Alignment options");
    argparser.AddArgument(&aligner_str, VALUE_TYPE_STRING, "", "aligner", "edlib",
                          "Specifies which library should be used for alignment. Options are:\n"
                          "    edlib       - Edit distance alignment, bit-vector approach. Fastest.\n"
                          "    ksw2-single - KSW2 library, single affine gap penalty.\n"
                          "    ksw2-double - KSW2 library, double affine gap penalty.",
                          0, "Alignment options");
    argparser.AddArgument(&parameters->do_diff, VALUE_TYPE_BOOL, "", "diff", "0",
                          "Align anchors without traceback. Useful for overlapping. Complementary to '--align'. The '--aligner' option does not affect this one.", 0,
                          "Alignment options");
    argparser.AddArgument(
                          &match, VALUE_TYPE_INT32, "", "match", "2",
                          "Match score for the DP alignment. Ignored when Edlib is used.",
                          0, "Alignment options");
    argparser.AddArgument(&mismatch, VALUE_TYPE_INT32, "", "mismatch", "4",
                          "Mismatch penalty for the DP alignment. Ignored when Edlib is used.",
                          0, "Alignment options");
    argparser.AddArgument(&gap_open_penalty_str, VALUE_TYPE_STRING, "", "gap-open", "2,1",
                          "Gap open penalty for the DP alignment. Ignored when Edlib is used.",
                          0, "Alignment options");
    argparser.AddArgument(&gap_ext_penalty_str, VALUE_TYPE_STRING, "", "gap-ext", "4,13",
                          "Gap extend penalty for the DP alignment. Ignored when Edlib is used.",
                          0, "Alignment options");
    argparser.AddArgument(&parameters->aligner_params->use_basic_cigar, VALUE_TYPE_BOOL, "",
                          "basic-cigar", "0", "Use the basic CIGAR format for output alignments.",
                          0, "Alignment options");
    argparser.AddArgument(&bandwidth, VALUE_TYPE_INT32, "",
                          "bandwidth", "-1",
                          "Bandwidth to use for alignment. -1 uses dynamic banding for Edlib, "
                          "and switches off banded alignment for other aligners. If not specified,"
                          "-1 will be used for Edlib, and 500 for ksw2-*.",
                          0, "Alignment options");
    argparser.AddArgument(&parameters->aligner_params->no_extend_alignment, VALUE_TYPE_BOOL, "",
                          "no-ext", "0",
                          "If specified, alignment will not be extended beyond the mapping boundary.",
                          0, "Alignment options");

    argparser.AddArgument(&zbandwidth, VALUE_TYPE_INT32, "",
                          "zbandwidth", "500",
                          "Bandwidth to use for extension alignment.",
                          0, "Alignment options");
    argparser.AddArgument(&end_bonus, VALUE_TYPE_INT32, "",
                          "end-bonus", "0",
                          "An extra end bonus score for KSW2-based alignment.",
                          0, "Alignment options");

    argparser.AddArgument(&parameters->aligner_params->aligner_opt.zdrop, VALUE_TYPE_INT32, "",
                          "zdrop", "400",
                          "The Z-drop heuristic for extension alignment.",
                          0, "Alignment options");
    argparser.AddArgument(&no_stop_ext_on_zero_score, VALUE_TYPE_BOOL, "",
                          "no-sezs", "0",
                          "By default, extension alignment will stop when the alignment score "
                          "drops below 0 (like local alignment). With this flag, even scores < 0 "
                          "will be considered with the zdrop break threshold. "
                          "Only valid for Edlib-based extension alignment.",
                          0, "Alignment options");

    // argparser.AddArgument(&parameters->max_evalue, VALUE_TYPE_DOUBLE, "",
    //                       "evalue", "1e0",
    //                       "Threshold for E-value. If E-value > FLT, read will be called unmapped. "
    //                       "If FLT < 0.0, threshold not applied.",
    //                       0, "Filtering options");
    argparser.AddArgument(&parameters->min_identity, VALUE_TYPE_DOUBLE, "",
                          "min-idt", "65.0",
                          "Minimum percent alignment identity allowed to report the alignment.",
                          0, "Filtering options");
    argparser.AddArgument(&parameters->bestn, VALUE_TYPE_INT64, "", "bestn", "5",
                          "Output best N alignments/mappings/overlaps. If <= 0 all secondary "
                          "mappings will be output. If 1, only primary will be output.",
                          0, "Filtering options");
    argparser.AddArgument(&parameters->bestn_threshold, VALUE_TYPE_DOUBLE, "",
                          "bestn-threshold", "0.01",
                          "If bestn != 1, all mappings with a score within this fraction from the "
                          "best score will be retained.",
                          0, "Filtering options");
    argparser.AddArgument(&min_secondary_to_primary_ratio, VALUE_TYPE_DOUBLE, "",
                          "sec-ratio", "0.80",
                          "Minimum secondary/primary score ratio to retain the secondary "
                          "mapping/alignment.",
                          0, "Filtering options");
    argparser.AddArgument(&allowed_suppl_overlap, VALUE_TYPE_INT32, "",
                          "suppl-ovl", "50",
                          "Allowed supplementary/primary overlap in bp. If 0, even a single bp"
                          "overlap will mark a region secondary instead of supplementary.",
                          0, "Filtering options");
    argparser.AddArgument(&parameters->min_map_len, VALUE_TYPE_INT64, "", "min-map-len", "0",
                          "Output only alignments/mappings/overlaps above this length.", 0,
                          "Filtering options");
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
    argparser.AddArgument(&no_relabel_secondary_supp, VALUE_TYPE_BOOL, "", "no-relabel", "0",
                          "Do not to relabel supplementary alignments.", 0,
                          "Other options");

    argparser.AddArgument(&parameters->debug_qid, VALUE_TYPE_INT64, "y", "debug-read", "-1",
                          "ID of the read to give the detailed verbose output.", 0,
                          "Debug options");
    argparser.AddArgument(&parameters->debug_qname, VALUE_TYPE_STRING, "Y", "debug-qname", "",
                          "QNAME of the read to give the detailed verbose output. Has precedence "
                          "over -y. Use quotes to specify.",
                          0, "Debug options");
    argparser.AddArgument(
        &parameters->verbose_output, VALUE_TYPE_INT64, "b", "verbose-out", "0",
        "Helpful debug comments can be placed in SAM output lines (at the end). Comments can be "
        "turned off by setting this parameter to 0. Different values increase/decrease verbosity "
        "level.\n0 - verbose off\n1 - server mode, command line will be omitted to obfuscate "
        "paths.\n2 - umm this one was skipped by accident. The same as 0.\n>=3 - detailed verbose "
        "is added for each alignment, including timing measurements and other.\n4 - qnames and "
        "rnames will not be trimmed to the first space.\n5 - QVs will be omitted (if available).",
        0, "Debug options");

    argparser.AddArgument(&parameters->composite, VALUE_TYPE_COMPOSITE, "x", "composite", "",
        "Pre-set parameters for different use cases. Valid options are:\n"
        " \"ovl-raw\"   - Equivalent to: '" + composite_str_ovl_raw + "'\n" +
        " \"ovl-hifi\"   - Equivalent to: '" + composite_str_ovl_hifi + "'\n"
        // + " \"ovl-miniasm\"   - Equivalent to: '" + composite_str_ovl_miniasm + "'\n", 0, "Composite options"
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
        VerboseShortHelpRaptor(argc, argv);
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
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    // Sanity check for the reads path.
    if (parameters->query_paths.size() == 0 && parameters->calc_only_index == false) {
        fprintf(stderr, "Please specify the path to the reads file.\n");
        fprintf(stderr, "\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    /////////////////////////////////////////////////
    /// Parsing and validating the input formats. ///
    /////////////////////////////////////////////////
    // Parse the input format.
    parameters->infmt = mindex::SequenceFormatFromString(infmt);
    if (parameters->infmt == mindex::SequenceFormat::Unknown) {
        fprintf(stderr, "Unknown input format '%s'!\n\n", infmt.c_str());
        VerboseShortHelpRaptorAndExit(argc, argv);
    }
    // Collect all FOFN files.
    parameters->query_paths = ExpandPathList(parameters->infmt, infmt, parameters->query_paths);
    // Validate the input files and formats.
    bool validate_rv1 = ValidateInputFiles(parameters->infmt, parameters->query_paths);
    if (validate_rv1 == false) {
        VerboseShortHelpRaptorAndExit(argc, argv);
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
        VerboseShortHelpRaptorAndExit(argc, argv);
    }
    // Collect all FOFN files.
    parameters->ref_paths = ExpandPathList(parameters->ref_fmt, infmt, parameters->ref_paths);
    // Validate the input files and formats.
    bool validate_rv = ValidateInputFiles(parameters->ref_fmt, parameters->ref_paths);
    if (validate_rv1 == false) {
        VerboseShortHelpRaptorAndExit(argc, argv);
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

    // Parse the input graph format.
    parameters->graph_fmt = raptor::GraphFormatFromString(graph_fmt);
    if (parameters->graph_fmt == GraphFormat::Unknown) {
        fprintf(stderr, "Unknown input graph format '%s'!\n\n", graph_fmt.c_str());
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    // Parse the output format.
    parameters->outfmt = raptor::OutputFormatFromString(outfmt);
    if (parameters->outfmt == OutputFormat::Unknown) {
        fprintf(stderr, "Unknown output format '%s'!\n\n", outfmt.c_str());
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    if (parameters->num_threads < 0) {
        parameters->num_threads = std::min(24, ((int)std::thread::hardware_concurrency()) / 2);
    }

    if (parameters->index_params->w <= 0) {
        fprintf(stderr, "Minimizer window length cannot be <= 0!\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    if (parameters->index_params->freq_percentil < 0.0 ||
        parameters->index_params->freq_percentil > 1.0) {
        fprintf(stderr, "Frequency percentil should be in range [0.0, 1.0].\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
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

    if (parameters->mapper_params->is_rna == true) {
        fprintf(stderr, "The RNA feature is unfortunately not implemented yet.\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    if (argparser.GetArgumentByLongName("graph")->is_set == true &&
        argparser.GetArgumentByLongName("agraph")->is_set == true) {
        fprintf(stderr,
                "Two types of input graph are specified. Use either '--graph' or '--agraph'.\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    parameters->add_symmetric_arcs = false;
    if (argparser.GetArgumentByLongName("agraph")->is_set == true) {
        parameters->add_symmetric_arcs = true;
    }

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
                VerboseShortHelpRaptorAndExit(argc, argv);
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
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    if (parameters->do_align && parameters->do_diff) {
        fprintf (stderr, "Both '--align' and '--diff' options are specified. These are mutually exclusive.\n");
        VerboseShortHelpRaptorAndExit(argc, argv);
    }

    // Parse the alignment parameters.
    {
        auto gap_open_segments = raptor::Tokenize(gap_open_penalty_str, ',');
        auto gap_ext_segments = raptor::Tokenize(gap_ext_penalty_str, ',');
        if (gap_open_segments.size() != gap_ext_segments.size()) {
            fprintf (stderr, "The gap open and gap ext parameters need to have same number of segments.\n");
            VerboseShortHelpRaptorAndExit(argc, argv);
        }
        if (gap_open_segments.size() == 0 || gap_ext_segments.size() == 0) {
            fprintf (stderr, "At least one gap penalty segment needs to be specified.\n");
            VerboseShortHelpRaptorAndExit(argc, argv);
        }

        std::string::size_type sz;
        std::vector<raptor::AffinePiece> affine_pieces;
        for (size_t i = 0; i < gap_open_segments.size(); i++) {
            float go = std::stof(gap_open_segments[i], &sz);
            float ge = std::stof(gap_ext_segments[i], &sz);
            affine_pieces.emplace_back(raptor::AffinePiece(-go, -ge));
        }

        raptor::PiecewisePenalties p(match, -mismatch, affine_pieces);
        parameters->aligner_params->aligner_opt.p = p;
        parameters->aligner_params->aligner_opt.bandwidth= bandwidth;
        parameters->aligner_params->aligner_opt.zbandwidth = zbandwidth;
        parameters->aligner_params->aligner_opt.end_bonus = end_bonus;
        parameters->aligner_params->aligner_opt.do_traceback = true;
        parameters->aligner_params->aligner_opt.gm = raptor::GlobalMargins();
        // Use the positive logic downstream. The negative logic is only for the command line parsing.
        parameters->aligner_params->aligner_opt.stop_ext_on_zero_score = !no_stop_ext_on_zero_score;

        parameters->aligner_params->aligner_type =  (aligner_str == "edlib") ? raptor::AlignerType::Edlib :
                                                    (aligner_str == "ksw2-single") ? raptor::AlignerType::KSW2Single :
                                                    (aligner_str == "ksw2-double") ? raptor::AlignerType::KSW2Double :
                                                    raptor::AlignerType::Unknown;

        if (parameters->aligner_params->aligner_type == raptor::AlignerType::Unknown) {
            fprintf (stderr, "Unknown aligner.\n");
            VerboseShortHelpRaptorAndExit(argc, argv);
        }

        // KSW2 tends to segfault on larger sequences if bandwidth is not specified
        // (most likely it has quadratic memory complexity). That's why we specify the
        // different default values for Edlib and KSW2 if the bandwidth was not
        // explicitly specified.
        if (argparser.GetArgumentByLongName("bandwidth")->is_set == false) {
            if (parameters->aligner_params->aligner_type == raptor::AlignerType::Edlib) {
                parameters->aligner_params->aligner_opt.bandwidth = -1;
            } else if (parameters->aligner_params->aligner_type == raptor::AlignerType::KSW2Single) {
                parameters->aligner_params->aligner_opt.bandwidth = 500;
            } else if (parameters->aligner_params->aligner_type == raptor::AlignerType::KSW2Double) {
                parameters->aligner_params->aligner_opt.bandwidth = 500;
            }
        }

        parameters->aligner_params->min_identity = parameters->min_identity;
        parameters->aligner_params->max_evalue = parameters->max_evalue;
    }

    parameters->mapper_params->verbose_level = parameters->verbose_level;
    parameters->mapper_params->debug_qid = parameters->debug_qid;
    parameters->mapper_params->debug_qname = parameters->debug_qname;

    parameters->aligner_params->verbose_level = parameters->verbose_level;
    parameters->aligner_params->debug_qid = parameters->debug_qid;
    parameters->aligner_params->debug_qname = parameters->debug_qname;
    parameters->aligner_params->is_rna = parameters->mapper_params->is_rna;

    parameters->relabel_secondary_supp = !no_relabel_secondary_supp;
    parameters->mapper_params->relabel_secondary_supp = !no_relabel_secondary_supp;
    parameters->aligner_params->relabel_secondary_supp = !no_relabel_secondary_supp;

    parameters->mapper_params->min_secondary_to_primary_ratio = min_secondary_to_primary_ratio;
    parameters->aligner_params->min_secondary_to_primary_ratio = min_secondary_to_primary_ratio;

    parameters->mapper_params->allowed_suppl_overlap = allowed_suppl_overlap;
    parameters->aligner_params->allowed_suppl_overlap = allowed_suppl_overlap;

#ifdef RAPTOR_TESTING_MODE
    if (parameters->verbose_level >= 5) {
        fprintf (stderr, "Compiled with RAPTOR_TESTING_MODE.\n");
    }
#endif

    return 0;
}

// int ProcessArgsRaptorIndex(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorIndex> parameters) {
//     bool help = false;
//     bool version = false;

//     ArgumentParser argparser;

//     std::string infmt, ref_fmt, graph_fmt, outfmt;
//     std::string index_region_string;

//     argparser.AddArgument(&parameters->ref_paths, VALUE_TYPE_STRING_LIST, "r", "ref", "",
//                           "Path to the sequences (fastq or fasta). Can be specified multiple times.", 0,
//                           "Input/Output options");
//     argparser.AddArgument(&parameters->out_path, VALUE_TYPE_STRING, "o", "out", "",
//                           "Path to the output file that will be generated.", 0,
//                           "Input/Output options");
//     // argparser.AddArgument(&ref_fmt, VALUE_TYPE_STRING, "", "ref-fmt", "auto",
//     //                       "Format of the reference sequence file. Options are:\n auto  - Determines the "
//     //                       "format automatically from file extension.\n fastq - Loads FASTQ or "
//     //                       "FASTA files.\n fasta - Loads FASTQ or FASTA files.\n gfa   - Graphical "
//     //                       "Fragment Assembly format.\n gfa1  - GFA-1 format.\n gfa2  - GFA-2 format\n sam   - Sequence Alignment/Mapping format.",
//     //                       0, "Input/Output options");
//     argparser.AddArgument(&parameters->batch_size, VALUE_TYPE_DOUBLE, "B", "batch", "200",
//                           "Queries will be loaded in batches of the size specified in megabytes. If there is a "
//                           "trailing 'x', then the batch size is in terms of coverage of reference genome. "
//                           "Value <= 0 loads the entire file.",
//                           0, "Input/Output options");
//     argparser.AddArgument(&(index_region_string), VALUE_TYPE_STRING, "", "region", "",
//                           "Indexes only the specified region. Region format: chr:start-end. Start "
//                           "is 0-based, and end is not inclusive. If end is <= 0, the entire suffix "
//                           "of the reference will be indexed.",
//                           0, "Index options");
//     argparser.AddArgument(
//         &(parameters->index_params->index_only_fwd_strand), VALUE_TYPE_BOOL, "", "fwd", "0",
//         "If specified, only the fwd strand of the reference will be indexed.", 0, "Index options");
//     argparser.AddArgument(&(parameters->index_params->k), VALUE_TYPE_INT32, "k", "seed-len", "15",
//                           "Length of the seeds used for hashing and lookup.", 0, "Index options");
//     argparser.AddArgument(&(parameters->index_params->w), VALUE_TYPE_INT32, "w", "minimizer-window",
//                           "5",
//                           "Length of the window to select a minimizer from. If equal to 1, "
//                           "minimizers will be turned off.",
//                           0, "Index options");
//     argparser.AddArgument(&(parameters->index_params->freq_percentil), VALUE_TYPE_DOUBLE, "",
//                           "freq-percentile", "0.0002",
//                           "Filer this fraction of most frequent seeds in the lookup process.", 0,
//                           "Index options");
//     argparser.AddArgument(&(parameters->index_params->min_occ_cutoff), VALUE_TYPE_INT32, "",
//                           "min-occ-cutoff", "0",
//                           "If the freq-percentil value is less than this count, this will be used "
//                           "instead for seed count filtering.",
//                           0, "Index options");
//     argparser.AddArgument(
//         &(parameters->index_params->homopolymer_suppression), VALUE_TYPE_BOOL, "", "hp-suppress",
//         "0", "If specified, homopolymer runs will be suppressed into a single seed point.", 0,
//         "Index options");
//     argparser.AddArgument(&(parameters->index_params->max_homopolymer_len), VALUE_TYPE_INT32, "",
//                           "max-hp-len", "5", "Compress only homopolymer runs below this length.", 0,
//                           "Index options");
//     argparser.AddArgument(
//         &parameters->keep_lowercase, VALUE_TYPE_BOOL, "", "keep-lowercase", "0",
//         "If this parameter is not specified, lowercase bases will be converted to uppercase.", 0,
//         "Index options");
//     argparser.AddArgument(&parameters->num_threads, VALUE_TYPE_INT64, "t", "threads", "-1",
//                           "Number of threads to use. If '-1', number of threads will be equal to "
//                           "min(24, num_cores/2).",
//                           0, "Other options");
//     argparser.AddArgument(
//         &parameters->verbose_level, VALUE_TYPE_INT64, "v", "verbose", "5",
//         "Verbose level. If equal to 0 nothing except strict output will be placed on stdout.", 0,
//         "Other options");

//     argparser.AddArgument(&version, VALUE_TYPE_BOOL, "", "version", "0", "Output the version and exit.", 0,
//                           "Other options");

//     argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0,
//                           "Other options");

//     argparser.ProcessArguments(argc, argv);

//     // Store the command line in parameters.
//     std::stringstream ss_command_line;
//     for (int i = 0; i < argc; i++) {
//         if (i > 0) ss_command_line << " ";
//         ss_command_line << argv[i];
//     }
//     parameters->command_line = ss_command_line.str();

//     // Check if help was triggered.
//     if (argparser.GetArgumentByLongName("help")->is_set == true) {
//         std::stringstream ss;
//         ss << SOFTWARE_NAME
//            << " - Indexing tool for Raptor - A very accurate and sensitive long-read, high error-rate sequence mapper\n"
//            << SOFTWARE_NAME << " ";
//         ss << "Version: " << RAPTOR_VERSION_STRING << "\n";
//         ss << "Build date: " << std::string(RELEASE_DATE).c_str() << "\n";
//         ss << "\n";
//         ss << LICENCE_INFORMATION << "\n";
//         ss << "\n";
//         ss << "Usage:\n";
//         ss << "\traptor [options] -r <reference_file> -q <reads_file> -o <output_sam_path>\n";
//         ss << "\n";

//         fprintf(stderr, "%s\n", ss.str().c_str());
//         fprintf(stderr, "%s\n", argparser.VerboseUsage().c_str());
//         exit(0);
//     }

//     if (version) {
//         std::stringstream ss;
//         ss << RAPTOR_VERSION_STRING << "\n";
//         fprintf(stdout, "%s", ss.str().c_str());
//         exit(0);
//     }

//     // Sanity check for the reference path.
//     if (parameters->ref_paths.size() == 0) {
//         fprintf(stderr, "Please specify the path to the reference file.\n");
//         fprintf(stderr, "\n");
//         VerboseShortHelpRaptorAndExit(argc, argv);
//     }
//     for (auto& ref_path: parameters->ref_paths) {
//         if (!raptor::FileExists(ref_path.c_str())) {
//             fprintf(stderr, "Reference does not exist: '%s'\n\n", ref_path.c_str());
//             VerboseShortHelpRaptorAndExit(argc, argv);
//         }
//     }

//     // // Parse the reference format.
//     // parameters->ref_fmt =
//     //     (ref_fmt == "auto")
//     //         ? mindex::SequenceFormat::Auto
//     //         : (ref_fmt == "fasta")
//     //             ? mindex::SequenceFormat::Fasta
//     //             : (ref_fmt == "fastq") ? mindex::SequenceFormat::Fastq
//     //             : (ref_fmt == "gfa") ? mindex::SequenceFormat::GFA
//     //             : (ref_fmt == "gfa1") ? mindex::SequenceFormat::GFA1
//     //             : (ref_fmt == "gfa2") ? mindex::SequenceFormat::GFA2
//     //             : (ref_fmt == "sam") ? mindex::SequenceFormat::SAM
//     //             : mindex::SequenceFormat::Unknown;
//     // if (parameters->ref_fmt == mindex::SequenceFormat::Unknown) {
//     //     fprintf(stderr, "Unknown reference format '%s'!\n\n", ref_fmt.c_str());
//     //     VerboseShortHelpRaptorAndExit(argc, argv);
//     // }

//     if (parameters->num_threads < 0) {
//         parameters->num_threads = std::min(24, ((int)std::thread::hardware_concurrency()) / 2);
//     }

//     if (parameters->index_params->w <= 0) {
//         fprintf(stderr, "Minimizer window length cannot be <= 0!\n");
//         VerboseShortHelpRaptorAndExit(argc, argv);
//     }

//     if (parameters->index_params->freq_percentil < 0.0 ||
//         parameters->index_params->freq_percentil > 1.0) {
//         fprintf(stderr, "Frequency percentil should be in range [0.0, 1.0].\n");
//         VerboseShortHelpRaptorAndExit(argc, argv);
//     }

//     // Write this out for every debug verbose level.
//     if (parameters->verbose_level > 5) {
//         fprintf(stderr, "%s\n", argparser.VerboseArguments().c_str());
//     }

//     parameters->index_params->is_region_specified = false;
//     if (index_region_string.size() > 0) {
//         parameters->index_params->is_region_specified = true;

//         size_t rname_end = index_region_string.find_last_of(":");
//         if (rname_end == std::string::npos) {
//             parameters->index_params->region_rname = index_region_string;
//             parameters->index_params->region_rstart = 0;
//             parameters->index_params->region_rend = -1;
//         } else {
//             parameters->index_params->region_rname = index_region_string.substr(0, rname_end);
//             size_t dash_pos = index_region_string.find_first_of("-", (rname_end + 1));
//             if (dash_pos == std::string::npos) {
//                 fprintf(stderr,
//                         "Region format is incorrect. Two coordinates need to be provided, "
//                         "separated by '-', in the format: chr:start-end.\n");
//                 VerboseShortHelpRaptorAndExit(argc, argv);
//             }
//             parameters->index_params->region_rstart =
//                 atoi(index_region_string.substr((rname_end + 1), dash_pos).c_str());
//             parameters->index_params->region_rend =
//                 atoi(index_region_string.substr((dash_pos + 1)).c_str());
//         }
//     }

//     return 0;
// }

} /* namespace raptor */
