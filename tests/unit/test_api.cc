#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <lib/argparser.h>
#include <index/sequence_file.h>
// #include <utility/utility_general.h>
#include <writer/raptor_results_writer.h>
#include <graph/segment_graph_parser.h>

#include <raptor/index_factory.h>
#include <params/params_raptor.h>
#include <raptor/raptor.h>
#include <raptor/mapper.h>
#include <raptor/graph_mapper.h>

#include <raptor/raptor_aligner.h>
#include <raptor/path_aligner.h>

#include <log/log_system.h>

///////////////////////////////////////////////
// Sample test data                         ///
///////////////////////////////////////////////
std::vector<std::string> refs_1 = {
                            "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
                            "TGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGG"
                            "TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTAC"
                            "ACAACATCCATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT"
                            };
std::vector<std::string> rnames_1 = {
                                "ecoli-0:240",
                                };
std::vector<std::string> queries_1 = {
                                "TCACTAAATACTTTAACCAATATAGGCATAGCGCACAGACAGATAAAAATTACAGAGTAC"
                                };
std::vector<std::string> qnames_1 = {
                                "ecoli-120:180",
                                };
///////////////////////////////////////////////

TEST(APIExamples, APITest1) {
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    // Create the index parameters;
    auto index_params = mindex::createIndexParams();
    index_params->k = 5;
    index_params->w = 1;

    // Create the index.
    mindex::IndexPtr index = std::move(mindex::createMinimizerIndex(index_params));
    index->AddSequences(refs_1, rnames_1);
    index->Build();

    // Create an empty graph.
    raptor::GraphPtr graph = raptor::createSegmentGraph();

	// Split the segment graph on edges.
	raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(graph);

    // Create a mapper.
    auto mapper_params = raptor::createParamsMapper();
    auto mapper = raptor::createMapper(index, mapper_params);

    // Create a graph mapper from the graph.
    auto graph_mapper = raptor::createGraphMapper(index, graph, ssg, mapper_params);

    // Create a writer object for output.
    std::ostringstream oss;
    auto writer = raptor::createRaptorResultsWriter(oss, index, raptor::OutputFormat::PAF);

    auto reads = mindex::createSequenceFile();
    reads->LoadAll(qnames_1, queries_1);

    for (auto& seq: reads->seqs()) {
        raptor::RaptorResults results;

        // Linear mapping.
        results.mapping_result = mapper->Map(seq);

        // Mapping results can at this point be used downstream. However,
        // to pull them through the writer, we need to first run graph mapping
        // (otherwise, the writter will not generate output).
        // The graph can be empty.
        results.graph_mapping_result = graph_mapper->Map(seq, results.mapping_result);

        results.regions = results.graph_mapping_result->CollectRegions(true);

        // results.mapping_result->Filter(1, 0.01, -1, false);
        writer->WriteSingleResult(reads, results, false, true, false);
    }

    std::ostringstream expected;
    expected << "ecoli-120:180\t60\t0\t56\t+\tecoli-0:240\t240\t120\t176\t60\t56\t60\tcm:i:12\tnm:i:-1\tas:i:60\tpi:i:0\tpj:i:0\tpn:i:1\tps:i:0\tcg:Z:*\n";

    ASSERT_EQ(oss.str(), expected.str());
}

TEST(APIExamples, FullPipelineWithValueAccessAPI) {
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ////////////////////////////
    /// Sample input data.   ///
    ////////////////////////////
    // The Index will be loaded from refs_1 and rnames_1.
    auto reads = mindex::createSequenceFile();
    reads->LoadAll(qnames_1, queries_1);

    ////////////////////////////
    /// Construct the Index. ///
    ////////////////////////////
    // Create the index parameters;
    auto index_params = mindex::createIndexParams();
    index_params->k = 5;
    index_params->w = 1;
    // Create the index.
    mindex::IndexPtr index = std::move(mindex::createMinimizerIndex(index_params));
    index->AddSequences(refs_1, rnames_1);
    index->Build();

    /////////////////////////////////
    /// Initialize an empty graph ///
    /////////////////////////////////
    raptor::GraphPtr graph = raptor::createSegmentGraph();

	// Split the segment graph on edges.
	raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(graph);

    /////////////////////////////////
    /// Create a linear mapper.   ///
    /////////////////////////////////
    auto mapper_params = raptor::createParamsMapper();
    auto mapper = raptor::createMapper(index, mapper_params);

    /////////////////////////////////
    /// Create a graph mapper.    ///
    /////////////////////////////////
    auto graph_mapper = raptor::createGraphMapper(index, graph, ssg, mapper_params);

    /////////////////////////////////
    /// Create an aligner.        ///
    /////////////////////////////////
    std::shared_ptr<raptor::ParamsAligner> aligner_params = raptor::createParamsAligner();
    // This is a global sequence aligner.
    auto aligner = raptor::createAligner(aligner_params->aligner_type, aligner_params->aligner_opt);
    // This is used to extend the alignments beyond the global matching boundary.
    auto ext_aligner = raptor::createAligner(raptor::AlignerType::KSW2Single, aligner_params->aligner_opt);
    // Create the Raptor aligner object that knows how to align paths.
    auto raptor_aligner = createRaptorAligner(index, aligner_params, aligner, aligner, ext_aligner);

    // Output stream.
    std::ostringstream oss;

    for (auto& seq: reads->seqs()) {
        raptor::RaptorResults results;

        ///////////////////////////////////////////////////////////////////////
        ///// Example: Linear alignment.                                  /////
        ///////////////////////////////////////////////////////////////////////
        // Linear mapping.
        results.mapping_result = mapper->Map(seq);

        // Write out all the mapped regions.
        std::vector<std::shared_ptr<raptor::RegionBase>> mapped_regions = results.mapping_result->CollectRegions(false);
        for (const auto& aln: mapped_regions) {
            oss << reads->GetSeqByAbsID(aln->QueryID())->header() << "\t"
                << aln->QueryLen() << "\t"
                << aln->QueryStart() << "\t"
                << aln->QueryEnd() << "\t"
                << (aln->TargetRev() ? "-" : "+") << "\t"
                << index->header(aln->TargetID()) << "\t"
                << aln->TargetLen() << "\t"
                << aln->TargetFwdStart() << "\t"
                << aln->TargetFwdEnd() << "\t"
                << aln->CoveredBasesQuery() << "\t"
                << (aln->TargetFwdEnd() - aln->TargetFwdStart()) << "\t"
                << 60
                << "\n";
        }


        ///////////////////////////////////////////////////////////////////////
        ///// Example: Graph Mapping.                                     /////
        ///////////////////////////////////////////////////////////////////////
        // Do the Graph Mapping.
        results.graph_mapping_result = graph_mapper->Map(seq, results.mapping_result);

        // Write out all the graph-mapped regions. Interfaces are the same.
        std::vector<std::shared_ptr<raptor::RegionBase>> graph_mapped_regions = results.graph_mapping_result->CollectRegions(false);
        for (const auto& aln: graph_mapped_regions) {
            oss << reads->GetSeqByAbsID(aln->QueryID())->header() << "\t"
                << aln->QueryLen() << "\t"
                << aln->QueryStart() << "\t"
                << aln->QueryEnd() << "\t"
                << (aln->TargetRev() ? "-" : "+") << "\t"
                << index->header(aln->TargetID()) << "\t"
                << aln->TargetLen() << "\t"
                << aln->TargetFwdStart() << "\t"
                << aln->TargetFwdEnd() << "\t"
                << aln->CoveredBasesQuery() << "\t"
                << (aln->TargetFwdEnd() - aln->TargetFwdStart()) << "\t"
                << 60
                << "\n";
        }

        ///////////////////////////////////////////////////////////////////////
        ///// Example: Alignment.                                         /////
        ///////////////////////////////////////////////////////////////////////
        // Align the graph-mapped regions.
        // Filter mappings, but leave room for error and mapq calculation.
        results.graph_mapping_result->Filter(-1, 0.20, 0, 0, false);

        results.aln_result = raptor_aligner->AlignPaths(seq, results.graph_mapping_result->paths());

        results.aln_result->Filter(10,      // bestn
                                   1.0,     // bestn_threshold
                                   0,       // min_map_len,
                                   0,       // min_mapq,
                                   75,    // min_identity
                                   false);  // only sort, and ignore other filters

        std::vector<std::shared_ptr<raptor::RegionBase>> aligned_graph_mapped_regions = results.aln_result->CollectRegions(false);
        for (const auto& aln: aligned_graph_mapped_regions) {
            oss << reads->GetSeqByAbsID(aln->QueryID())->header() << "\t"
                << aln->QueryLen() << "\t"
                << aln->QueryStart() << "\t"
                << aln->QueryEnd() << "\t"
                << (aln->TargetRev() ? "-" : "+") << "\t"
                << index->header(aln->TargetID()) << "\t"
                << aln->TargetLen() << "\t"
                << aln->TargetFwdStart() << "\t"
                << aln->TargetFwdEnd() << "\t"
                << aln->CoveredBasesQuery() << "\t"
                << (aln->TargetFwdEnd() - aln->TargetFwdStart()) << "\t"
                << 60
                << "\n";
        }
        ///////////////////////////////////////////////////////////////////////

    }

    std::ostringstream expected;
    expected << "ecoli-120:180\t60\t0\t56\t+\tecoli-0:240\t240\t120\t176\t60\t56\t60\n";
    expected << "ecoli-120:180\t60\t0\t56\t+\tecoli-0:240\t240\t120\t176\t60\t56\t60\n";
    expected << "ecoli-120:180\t60\t0\t60\t+\tecoli-0:240\t240\t120\t180\t60\t60\t60\n";

    ASSERT_EQ(oss.str(), expected.str());
}
