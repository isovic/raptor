#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <raptor_sova/args_raptor_sova.h>
#include <raptor_sova/raptor_sova_worker_queue.h>
#include <utility/memtime.h>

#include <raptor/yield_index.h>
#include <writer/raptor_results_writer_factory.h>
#include <sequences/sequence_file_composite_factory.h>
#include <sequences/random_access_sequence_file.h>
#include <utility/tictoc.h>

#include <utility/stringutil.h>
void PrintAsM4(FILE *fp_out, const mindex::IndexPtr& index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping) {
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    int32_t score = mapping->Score();
    int32_t cov_bases_q = mapping->CoveredBasesQuery();
    int32_t cov_bases_t = mapping->CoveredBasesTarget();

    double qspan = q_end - q_start;
    double tspan = t_end - t_start;
    double identity_q = (qspan != 0.0) ? (static_cast<double>(cov_bases_q)) / (qspan) : -1.0;
    double identity_t = (tspan != 0.0) ? (static_cast<double>(cov_bases_t)) / (tspan) : -1.0;
    double identity = std::min(identity_q, identity_t);

    if (edit_dist >= 0) {
        double edit_dist_double = edit_dist;
        identity_q = (qspan != 0) ? ((qspan  - edit_dist_double) / qspan) : -2.0;
        identity_t = (tspan != 0) ? ((tspan - edit_dist_double) / tspan) : -2.0;
        identity = std::min(identity_q, identity_t);
    }

    if (t_is_rev) {
        // Output in the fwd strand always.
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

    // std::cout  << q_name << " "
    //     << t_name << " "
    //     << -score << " "
    //     << 100.0 * jaccard << " "
    //     << 0 << " "
    //     << q_start << " "
    //     << q_end << " "
    //     << q_len << " "
    //     << ((t_is_rev) ? "1" : "0") << " "
    //     << t_start << " "
    //     << t_end << " "
    //     << t_len
	// 	<< "\n";

	fprintf(fp_out, "%s %s %d %.2lf %d %d %d %d %d %d %d %d\n",
			q_name.c_str(), t_name.c_str(), score, 100.0 * identity,
			0, q_start, q_end, q_len,
			static_cast<int32_t>(t_is_rev), t_start, t_end, t_len);

}

void AdHocWriteBatch(const mindex::IndexPtr& index, const mindex::SequenceFilePtr seqs,
						const std::vector<std::unique_ptr<raptor::RaptorResults>>& results,
						bool write_custom_tags) {

    for (auto& result: results) {
		if (result == nullptr) {
			continue;
		}

		bool do_output = !result->regions().empty();
		// std::string timings_all = OutputFormatter::TimingMapToString(result->timings());
		// int32_t mapq = result->mapq();
		const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result->regions();
		// int64_t q_id = result->q_id();

        for (size_t i = 0; i < regions_to_write.size(); i++) {
            auto aln = regions_to_write[i];
            bool is_secondary = (aln->PathId() > 0);
            bool is_supplementary = (aln->SegmentId() > 0);
            auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());

			PrintAsM4(stdout, index, qseq, aln);

            // if (outfmt_ == raptor::OutputFormat::SAM) {
            //     *oss_ptr_ << OutputFormatter::ToSAM(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            // } else if (outfmt_ == raptor::OutputFormat::PAF) {
            //     *oss_ptr_ << OutputFormatter::ToPAF(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            // } else if (outfmt_ == raptor::OutputFormat::GFA2) {
            //     *oss_ptr_ << OutputFormatter::ToGFA2Edge(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            // } else if (outfmt_ == raptor::OutputFormat::MHAP) {
            //     *oss_ptr_ << OutputFormatter::ToMHAP(index_, qseq, aln, mapq);
            // } else if (outfmt_ == raptor::OutputFormat::M4) {
            //     *oss_ptr_ << OutputFormatter::ToM4(index_, qseq, aln, mapq);
            // }
        }
    }
}

void RunTool(std::shared_ptr<raptor::sova::ParamsRaptorSova> parameters) {
	TicToc tt_all;

	LOG_ALL("Num threads: %ld\n", parameters->num_threads);

	/*
	* Create the index.
	*/
	LOG_ALL("Creating the index.\n");
	auto index = raptor::YieldIndex(parameters->ref_paths, parameters->ref_fmt, "temp.rai", true,
								parameters->index_on_the_fly,
								true,
								parameters->rdb_block_ref,
								parameters->index_type,
								parameters->index_params);
	LOG_ALL("Memory usage: %.2f GB. Time: %.2f/%.2f sec.\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0), tt_all.get_secs(true), tt_all.get_cpu_secs(true));

	// Create a writer for results. If "-" or empty, writing is to stdout.
	auto writer = raptor::createRaptorResultsWriter(parameters->out_path, index, parameters->outfmt);

	int64_t total_processed = 0;

	LOG_ALL("Creating the Raptor object.\n");
	auto mapper = raptor::sova::createRaptorSovaWorkerQueue(index, parameters);

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
			mapper->Clear();

			TicToc tt_batch_align;
			mapper->Align(reads);
			tt_batch_align.stop();

			TicToc tt_batch_write;
			AdHocWriteBatch(index, reads, mapper->results(), true);
			// writer->WriteBatch(reads, mapper->results(), parameters->do_align, parameters->strict_format == false, parameters->one_hit_per_target);
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
			mapper->Clear();

			TicToc tt_batch_align;
			mapper->Align(reads);
			tt_batch_align.stop();

			TicToc tt_batch_write;
			AdHocWriteBatch(index, reads, mapper->results(), true);
			// writer->WriteBatch(reads, mapper->results(), parameters->do_align, parameters->strict_format == false, parameters->one_hit_per_target);
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

	auto parameters = raptor::sova::createParamsRaptorSova();

	if (raptor::sova::ProcessArgsRaptorSova(argc, argv, parameters)) {
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
