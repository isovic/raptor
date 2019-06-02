#ifndef RUN_ALL_TESTS_

#include <string>
#include <vector>

#include <log/log_tools.h>
#include <lib/argparser.h>
#include <version.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <args.h>
#include <limits>
#include <sequences/sequence_file.h>
#include <sequences/sequence_file_composite_fofn.h>
#include <sequences/sequence_serializer.h>
#include <utility/memtime.h>
#include <utility/stringutil.h>
#include <raptor_reshape/params_raptor_reshape.h>

/*
	RaptorDB v0.2 specification.
	Each line is a separate entity. Order of lines does not matter. Lines start with a keyword character,
	depending on that encode:
	- Version: "V\t<float-version_num>"
	- Files: "F\t<int64-file_id>\t<str-file_path>\t<str-file_format>"
	- Sequences: "S\t<int64-seq_id>\t<str-header>\t<int64-seq_len>\t<int64-file_id>\t<int64-start_offset_in_file>\t<int64-data_length>"
	- Blocks: "B\t<int64-block_id>\t<int64-seq_id_start>\t<int64-seq_id_end>\t<int64-num_bases_in_block>"
*/

void RunTool(std::shared_ptr<raptor::ParamsRaptorReshape> parameters) {
	LOG_ALL("Raptor-reshape starts processing.\n");

	if (parameters->in_paths.empty()) {
		FATAL_REPORT(ERR_UNEXPECTED_VALUE, "No input files specified.");
	}

	int64_t load_batch_size_int64 = static_cast<int64_t>(parameters->in_batch_size);
	int64_t write_block_size_int64 = (parameters->block_size > 0) ? (static_cast<int64_t>(parameters->block_size * 1024 * 1024)) : std::numeric_limits<int64_t>::max();

	auto seq_file_parser = mindex::createSequenceFileCompositeFofn(parameters->in_paths, parameters->in_fmt);
	auto seq_file = mindex::createSequenceFile();

	int64_t total_num_blocks = 0;
	int64_t total_out_files = 0;
	int64_t total_num_bases = 0;
	int64_t total_num_seqs = 0;
	int64_t block_num_bases = 0;
	int64_t block_start_seq_id = 0;

	std::string out_path = parameters->out_prefix + ".block." + std::to_string(total_out_files);
	std::ofstream ofs;

	if (parameters->symlink_files == false) {
		ofs.open(out_path);
	}

	std::string out_rdb_path = parameters->out_prefix + ".rdb";
	std::ofstream ofs_rdb(out_rdb_path);

	std::string out_format = mindex::SequenceFormatToString(parameters->out_fmt);

	ofs_rdb << "V\t0.2\n";

	if (parameters->symlink_files == false) {
		ofs_rdb << "F\t" << total_out_files << "\t" << out_path << "\t" << out_format << "\n";
	}

	int64_t file_before = -1;

	while ((seq_file = seq_file_parser->YieldBatchOfOne())) {
		int64_t pos_now = seq_file_parser->GetCurrentFileOffset();
		int64_t pos_last = seq_file_parser->GetCurrentFilePreviousOffset();
		int64_t file_now = seq_file_parser->GetCurrentFileId();
		int64_t file_id = (parameters->symlink_files) ? file_now : total_out_files;

		if (pos_now < 0) {
			WARNING_REPORT(ERR_UNEXPECTED_VALUE, "pos_now = %ld\n", pos_now);
		}

		if (seq_file->seqs().size() == 0) {
			break;
		}

		const auto& seq = seq_file->seqs()[0];
		int64_t curr_bases = seq->len();
		total_num_bases += curr_bases;
		block_num_bases += curr_bases;

		if (parameters->rename_seqs && parameters->symlink_files == false) {
			std::ostringstream ss_buff;
			ss_buff << std::setfill('0') << std::setw(9) << seq->abs_id();
			seq->header(ss_buff.str());
		}

		if (parameters->symlink_files == false) {	// WITH dumping of sequences to disk.
			// Overwrite the pos_before and pos_after, if we're writing the seqs to disk.
			int64_t num_written_bytes = static_cast<int64_t>(ofs.tellp());
			mindex::SequenceSerializer::SerializeSequence(ofs, seq, parameters->out_fmt);
			int64_t num_written_bytes_after = static_cast<int64_t>(ofs.tellp());

			ofs_rdb << "S\t" << total_num_seqs << "\t" << raptor::TrimToFirstWhiteSpace(seq->header()) << "\t" << seq->len()
					<< "\t" << file_id << "\t" << num_written_bytes << "\t" << (num_written_bytes_after - num_written_bytes) << "\n";

		} else {	// WITHOUT writing to disk, just pointing to existing locations.
			if (file_now != file_before) {
				// In case of symlinking, the out format has to be the same as input format.
				auto curr_out_fmt = mindex::GetSequenceFormatFromPath(parameters->in_paths[file_now]);
				std::string out_format = mindex::SequenceFormatToString(curr_out_fmt);
				ofs_rdb << "F\t" << file_now << "\t" << parameters->in_paths[file_now] << "\t" << out_format << "\n";
			}

			ofs_rdb << "S\t" << total_num_seqs << "\t" << raptor::TrimToFirstWhiteSpace(seq->header()) << "\t" << seq->len()
					<< "\t" << file_id << "\t" << pos_last << "\t" << (pos_now - pos_last) << "\n";
		}

		if (block_num_bases > write_block_size_int64) {
			LOG_ALL("Block done! Memory usage: %.2f GB\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
			fflush(stdout);

			ofs_rdb << "B\t" << total_num_blocks << "\t" << block_start_seq_id << "\t" << (total_num_seqs + 1) << "\t" << block_num_bases << "\n";
			block_start_seq_id = total_num_seqs + 1;	// One after the current. The current is overflowing the write_block_size_int64 and included in that block.

			++total_num_blocks;
			block_num_bases = 0;

			if (parameters->split_blocks && parameters->symlink_files == false) {
				++total_out_files;
				out_path = parameters->out_prefix + ".block." + std::to_string(total_out_files);
				ofs = std::ofstream(out_path);
				ofs_rdb << "F\t" << total_out_files << "\t" << out_path << "\t" << out_format << "\n";
			}
		}

		++total_num_seqs;
		file_before = file_now;
	}
	if (block_start_seq_id != total_num_seqs) {
		ofs_rdb << "B\t" << total_num_blocks << "\t" << block_start_seq_id << "\t" << total_num_seqs << "\t" << block_num_bases << "\n";
	}

	LOG_ALL("Memory usage: %.2f GB\n", ((double) raptor::getPeakRSS()) / (1024.0 * 1024.0 * 1024.0));
	LOG_ALL("Done!\n");
}

int Run(int argc, char **argv) {

	std::string program_name(argv[0]);

	auto parameters = raptor::createParamsRaptorReshape();

	if (raptor::ProcessArgsRaptorReshape(argc, argv, parameters)) {
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
