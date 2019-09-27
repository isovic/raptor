/*
 * yield_index.cc
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#include <raptor/yield_index.h>
// #include <utility/utility_general.h>
#include <sequences/sequence_file.h>
#include <log/log_tools.h>
#include <params/params_raptor.h>
#include <utility/files.hpp>
#include <sequences/random_access_sequence_file.h>
#include <sequences/sequence_file_composite_factory.h>
#include <index/index_factory.h>

namespace raptor {

std::shared_ptr<mindex::IndexBase> YieldIndex(const std::vector<std::string>& ref_paths, mindex::SequenceFormat ref_fmt,
                                                const std::string& index_path, bool rebuild_index,
                                                bool index_on_the_fly, bool auto_rebuild_index,
                                                int64_t rdb_block_id,
                                                mindex::IndexType index_type,
                                                std::shared_ptr<mindex::IndexParams> index_params) {

    auto index = mindex::createIndex(index_type, index_params);

    bool load = !rebuild_index && raptor::FileExists(index_path);
    bool store = false;
    // bool store = !index_on_the_fly; // Not available for now.

    if (load) {
        LOG_ALL("Loading index from file: '%s'.\n", index_path.c_str());

        int ret_load = index->Load(index_path);

        if (ret_load) {
            LOG_DEBUG("Problems loading index: index is of wrong version or index file corrupt.\n");

            if (auto_rebuild_index) {
                LOG_ALL("Rebuilding the index.\n");
                load = false;
            } else {
                FATAL_REPORT(
                    ERR_GENERIC,
                    "Not rebuilding the index automatically (specify --auto-rebuild-index).");
            }
        }
    }

    // Separate 'if' because there can be a fallthrough case when index needs to be rebuilt.
    if (!load) {
        LOG_ALL("Building the index for k = %d, w = %d, freq_percentile = %f\n", index_params->k,
                index_params->w, index_params->freq_percentil);

        if (ref_fmt != mindex::SequenceFormat::RaptorDB) {
            for (auto& ref_path: ref_paths) {
                LOG_ALL("Loading reference sequences from file: '%s'.\n", ref_path.c_str());
                auto refs_parser = mindex::createSequenceFileCompositeFactory({ref_path}, ref_fmt);
                auto refs = mindex::createSequenceFile();
                while ((refs = refs_parser->YieldBatchOfOne()) != nullptr) {
                    for (auto& ref : refs->seqs()) {
                        index->AddSequence(ref->data(), ref->header());
                    }
                }
            }

            LOG_ALL("Building the index.\n");
            index->BuildIndex();

            LOG_ALL("Number of sequences in the index: %ld\n", index->num_seqs());
            LOG_ALL("Finished building index.\n");

            if (store) {
                LOG_ALL("Storing the index to file: '%s'.\n", index_path.c_str());
                index->Store(index_path);
            }
        } else {
            LOG_ALL("Loading reference sequences from RaptorDB file: '%s'.\n", ref_paths[0].c_str());
            auto rasf = mindex::createRandomAccessSequenceFile(ref_paths[0], 50);
            int64_t num_blocks = static_cast<int64_t>(rasf->db_blocks().size());

            if (num_blocks > 0 && rdb_block_id >= num_blocks) {
                FATAL_REPORT(
                    ERR_UNEXPECTED_VALUE,
                    "Specified block_id is larger than the number of blocks in RaptorDB. block_id = %ld, rasf->db_blocks().size() = %lu.", rdb_block_id, rasf->db_blocks().size());
            }
            mindex::SequenceFilePtr seq_file = nullptr;
            if (rdb_block_id >= 0) {
                seq_file = rasf->FetchBlock(rdb_block_id);
            } else {
                seq_file = rasf->FetchAll();
            }

            LOG_ALL("Building the index.\n");
            index->SetSequenceFile(seq_file);
            index->BuildIndex();

            LOG_ALL("Number of sequences in the index: %ld\n", index->num_seqs());
            LOG_ALL("Finished building index.\n");

            if (store) {
                LOG_ALL("Storing the index to file: '%s'.\n", index_path.c_str());
                index->Store(index_path);
            }
        }
    }

    if (index == nullptr) {
        FATAL_REPORT(ERR_GENERIC, "No index was generated! Exiting.");
    }

    return index;
}

}  // namespace raptor
