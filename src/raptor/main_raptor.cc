#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <raptor/raptor.h>
#include <params/params_raptor.h>
#include <raptor/yield_index.h>
#include <log/log_tools.h>
#include <sequences/sequence_file.h>
#include <sequences/sequence_file_composite_fofn.h>
#include <sequences/sequence_file_composite_pbxml.h>
#include <sequences/sequence_file_composite_factory.h>
#include <sequences/random_access_sequence_file.h>
#include <lib/argparser.h>
#include <version.h>
#include <writer/raptor_results_writer_factory.h>
#include <graph/segment_graph_parser.h>
#include <graph/split_segment_graph.h>
#include <iostream>
#include <fstream>
#include <raptor/args_raptor.h>

#include <raptor/mapper.h>
#include <raptor/graph_mapper.h>

#include <utility/memtime.h>
#include <utility/tictoc.h>

void APIExample() {
    std::vector<std::string> ref_paths = {std::string("test_data/ref/ecoli_K12_MG1655_U00096.3.fasta")};
    std::string reads_path("test_data_link/ecoli-public-p6c4/m141013_011508_sherri_c100709962550000001823135904221533_s1_p0.subreads.fasta");

    // Create an index.
    auto index_params = mindex::createParamsIndex();
    auto index = raptor::YieldIndex(ref_paths, mindex::SequenceFormat::Auto, std::string("temp.rai"), true,
                                true, true, -1, mindex::IndexType::Minimizer, index_params);

    // Create a mapper.
    auto mapper_params = raptor::createParamsMapper();
    auto mapper = raptor::createMapper(index, mapper_params);

    // Create an empty graph.
    raptor::GraphPtr graph = raptor::createSegmentGraph();

	// Split the segment graph on edges.
	raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(graph);

    // Create a graph mapper from the graph.
    auto graph_mapper = raptor::createGraphMapper(index, graph, ssg, mapper_params);

    // Create a writer object for output to std::cout.
    auto writer = raptor::createRaptorResultsWriter("-", index, raptor::OutputFormat::PAF);

    // Process reads one by one in a single thread.
	auto seq_file_parser = mindex::createSequenceFileCompositeFofn({reads_path}, mindex::SequenceFormat::Fastq);
	auto reads = seq_file_parser->YieldAll();

    for (auto& seq: reads->seqs()) {
        auto mapping_result = mapper->Map((seq));

        // Mapping results can at this point be used downstream. However,
        // to pull them through the writer, we need to first run graph mapping
        // (otherwise, the writter will not generate output).
        auto graph_mapping_result = graph_mapper->Map((seq), mapping_result);

        auto regions = graph_mapping_result->CollectRegions(false);
        auto mapq = graph_mapping_result->CalcMapq();
        auto results = raptor::createRaptorResults(seq->id(), regions, {}, mapq);

        // results.mapping_result->Filter(1, 0.01, -1, false);
        writer->WriteSingleResult(reads, results, false, true, false);
    }

}

void RunRaptor(std::shared_ptr<raptor::ParamsRaptor> parameters) {
	TicToc tt_all;

	LOG_ALL("Num threads: %ld\n", parameters->num_threads);

	/*
	* Create the index.
	*/
	LOG_ALL("Creating the index.\n");
	std::string index_path = "temp.rai"; // parameters->ref_path + std::string(".rai");

	auto index = raptor::YieldIndex(parameters->ref_paths, parameters->ref_fmt, index_path, parameters->rebuild_index,
								parameters->index_on_the_fly,
								parameters->auto_rebuild_index,
								parameters->rdb_block_ref,
								parameters->index_type,
								parameters->index_params);

	LOG_ALL("Memory usage: %.2f GB. Time: %.2f/%.2f sec.\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0), tt_all.get_secs(true), tt_all.get_cpu_secs(true));

	/*
	* Construct the graph.
	*/
	// raptor::GraphPtr graph(nullptr);
	raptor::GraphPtr graph = std::move(raptor::createSegmentGraph());
	raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph();
	if (parameters->graph_path.size() > 0) {
		LOG_ALL("Loading graph from '%s'.\n", parameters->graph_path.c_str());
		graph = raptor::GraphLoader::FromGFA(parameters->graph_path, index, parameters->add_symmetric_arcs);
		// ssg = raptor::createSplitSegmentGraph(graph);

		// #ifdef RAPTOR_TESTING_MODE
		// 	if (parameters->verbose_level >= 9) {
		// 		std::cerr << "Loaded graph:" << std::endl << graph->Verbose() << std::endl << std::endl;
		// 		std::cerr << "Split segment graph:" << std::endl << ssg->Verbose() << std::endl << std::endl;

		// 		std::ofstream ofs_debug_graph_sg("temp/debug-graph/segment_graph.html");
        //         if (ofs_debug_graph_sg.is_open()) {
		// 			graph->WriteAsHTML(ofs_debug_graph_sg);
		// 		}

		// 		std::ofstream ofs_debug_graph_ssg("temp/debug-graph/split_segment_graph.html");
        //         if (ofs_debug_graph_ssg.is_open()) {
		// 			ssg->WriteAsHTML(ofs_debug_graph_ssg);
		// 		}
		// 	}
		// #endif
	}

	// Create a writer for results. If "-" or empty, writing is to stdout.
	auto writer = raptor::createRaptorResultsWriter(parameters->out_path, index, parameters->outfmt);

	int64_t total_processed = 0;

	LOG_ALL("Creating the Raptor object.\n");
	auto raptor = raptor::createRaptor(index, graph, ssg, parameters, true);

	// The RandomAccessSequenceFile has a slightly different interface compared
	// to the Composite sequence files, and that's why here we have a separate
	// branch for both options. This is becauase for a RandomAccessSequenceFile
	// batch size is equal to the block size.
	if (parameters->infmt == mindex::SequenceFormat::RaptorDB) {
		// If the input is a RaptorDB, then a batch size is one block of the input.
		const std::string& query_path = parameters->query_paths[0];
		LOG_ALL("Processing reads from: %s\n", query_path.c_str());

		// Open the input DB file.
		mindex::RandomAccessSequenceFilePtr random_seqfile = mindex::createRandomAccessSequenceFile(query_path, 50);

		int64_t num_blocks = static_cast<int64_t>(random_seqfile->db_blocks().size());
		int64_t start_block_id = 0;
		int64_t end_block_id = num_blocks;

		if (parameters->rdb_block_query >= 0) {
			start_block_id = parameters->rdb_block_query;
			end_block_id = parameters->rdb_block_query + 1;
			LOG_ALL("Only block %ld will be processed.\n", start_block_id);
		}

		if (num_blocks > 0 && start_block_id >= num_blocks) {
			FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The start_block_id >= num_blocks. start_block_id = %ld, num_blocks = %ld. Exiting.", start_block_id, num_blocks);
		}
		if (num_blocks > 0 && end_block_id > num_blocks) {	// Non-inclusive end.
			FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The end_block_id >= num_blocks. end_block_id = %ld, num_blocks = %ld. Exiting.", end_block_id, num_blocks);
		}

		// Write the header if needed (e.g. parameters->outfmt == raptor::OutputFormat::SAM).
		const auto& header_groups = random_seqfile->GetHeaderGroups();
		writer->WriteHeader(header_groups);

		// Process all blocks.
		for (int64_t block_id = start_block_id; block_id < end_block_id; ++block_id) {
			TicToc tt_batch_total;
			TicToc tt_batch_load;
			LOG_ALL("Loading block %ld from the RaptorDB.\n", block_id);
			mindex::SequenceFilePtr reads = random_seqfile->FetchBlock(block_id);
			tt_batch_load.stop();

			LOG_ALL("Mapping batch of %ld reads (%.2lf MB). Processed %ld reads so far.\n",
						reads->seqs().size(), reads->total_size() / (1048576.0), total_processed);
			raptor->Clear();

			TicToc tt_batch_align;
			raptor->Align(reads);
			tt_batch_align.stop();

			TicToc tt_batch_write;
			writer->WriteBatch(reads, raptor->results(), parameters->do_align, parameters->strict_format == false, parameters->one_hit_per_target);
			tt_batch_write.stop();

			total_processed += reads->seqs().size();
			if (parameters->num_reads_to_process > 0 &&
				total_processed > (parameters->start_read + parameters->num_reads_to_process)) {
				break;
			}
			tt_batch_total.stop();
			LOG_ALL("Batch done! Memory usage: %.2f GB. Time: a:%.2f/%.2f, l:%.2f/%.2f, w:%.2f/%.2f, t:%.2f/%.2f\n",
						((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0),
						tt_batch_align.get_secs(), tt_batch_align.get_cpu_secs(),
						tt_batch_load.get_secs(), tt_batch_load.get_cpu_secs(),
						tt_batch_write.get_secs(), tt_batch_write.get_cpu_secs(),
						tt_batch_total.get_secs(), tt_batch_total.get_cpu_secs());
			fflush(stdout);
		}

	} else {
		LOG_ALL("Processing reads.\n");
		int64_t batch_size_in_mb = static_cast<int64_t>(parameters->batch_size);
		if (parameters->batch_type == mindex::BatchLoadType::Coverage) {
			batch_size_in_mb = static_cast<int64_t>(std::ceil(index->total_len() * parameters->batch_size / (1024.0 * 1024.0)));
			LOG_ALL("Batch loading input queries by reference coverage %.4lf.\n", parameters->batch_size);
			LOG_ALL("Rough batch size: %ld MB.\n", batch_size_in_mb);
		} else {
			LOG_ALL("Batch loading input queries by size specified in MB.\n");
			LOG_ALL("Rough batch size: %ld MB.\n", batch_size_in_mb);
		}

		// Open a sequence file object which will iterate through all inputs.
		mindex::SequenceFileCompositeBasePtr reads_parser = mindex::createSequenceFileCompositeFactory(parameters->query_paths, parameters->infmt);

		// Write the header if needed (e.g. parameters->outfmt == raptor::OutputFormat::SAM).
		const auto& header_groups = reads_parser->GetHeaderGroups();
		writer->WriteHeader(header_groups);

		// Iterate through all the batches.
		TicToc tt_batch_total;
		TicToc tt_batch_load;
		mindex::SequenceFilePtr reads = mindex::createSequenceFile();
		while ((reads = reads_parser->YieldBatchMB(batch_size_in_mb))) {
			tt_batch_load.stop();

			LOG_ALL("Mapping batch of %ld reads (%.2lf MB). Processed %ld reads so far.\n",
						reads->seqs().size(), reads->total_size() / (1048576.0), total_processed);
			raptor->Clear();

			TicToc tt_batch_align;
			raptor->Align(reads);
			tt_batch_align.stop();

			TicToc tt_batch_write;
			writer->WriteBatch(reads, raptor->results(), parameters->do_align, parameters->strict_format == false, parameters->one_hit_per_target);
			tt_batch_write.stop();

			total_processed += reads->seqs().size();
			if (parameters->num_reads_to_process > 0 &&
				total_processed > (parameters->start_read + parameters->num_reads_to_process)) {
				break;
			}
			tt_batch_total.stop();
			LOG_ALL("Batch done! Memory usage: %.2f GB. Time: a:%.2f/%.2f, l:%.2f/%.2f, w:%.2f/%.2f, t:%.2f/%.2f\n",
						((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0),
						tt_batch_align.get_secs(), tt_batch_align.get_cpu_secs(),
						tt_batch_load.get_secs(), tt_batch_load.get_cpu_secs(),
						tt_batch_write.get_secs(), tt_batch_write.get_cpu_secs(),
						tt_batch_total.get_secs(), tt_batch_total.get_cpu_secs());
			fflush(stdout);
			tt_batch_total.start();
			tt_batch_load.start();
		}
	}

	tt_all.stop();
	LOG_ALL("Memory usage: %.2f GB. Time: %.2f/%.2f sec.\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0), tt_all.get_secs(), tt_all.get_cpu_secs());
	LOG_ALL("Done!\n");
}

int Run(int argc, char **argv) {

	std::string program_name(argv[0]);

	auto parameters = raptor::createParamsRaptor();

	if (raptor::ProcessArgsRaptor(argc, argv, parameters)) {
		return 1;
	}

	LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters->verbose_level);

	if (parameters->verbose_level == 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
	} else if (parameters->verbose_level > 1) {
		LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
	}
	fflush(stdout);

	RunRaptor(parameters);

	return 0;
}

int main(int argc, char **argv) {
  	return Run(argc, argv);
}

#endif
