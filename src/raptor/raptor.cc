/*
 * raptor.cc
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#include <thread>
#include <raptor/raptor.h>
#include <lib/thread_pool.hpp>
#include <log/log_tools.h>
#include <aligner/aligner_containers.h>
#include <aligner/aligner_factory.h>
#include <raptor/mapper.h>
#include <raptor/raptor_aligner.h>
#include <raptor/graph_mapper.h>
#include <raptor/path_aligner.h>
// #include <containers/single_aln_region.h>
#include <utility/range_tools.hpp>
#include <utility/stringutil.h>
#include <aligner/difflib_edlib.h>

namespace raptor {

std::unique_ptr<Raptor> createRaptor(const mindex::IndexPtr index, const raptor::GraphPtr graph,
                                     const raptor::SplitSegmentGraphPtr ssg,
                                     const std::shared_ptr<ParamsRaptor> params,
                                     bool cleanup_target_hits) {
    return std::unique_ptr<Raptor>(new Raptor(index, graph, ssg, params, cleanup_target_hits));
}

Raptor::Raptor(const mindex::IndexPtr index, const raptor::GraphPtr graph,
               const raptor::SplitSegmentGraphPtr ssg,
               const std::shared_ptr<ParamsRaptor> params,
               bool cleanup_target_hits)
    : params_(params),
      cleanup_target_hits_(cleanup_target_hits),
      index_(index),
      graph_(graph),
      ssg_(ssg),
      results_(),
      mapping_status_(RaptorReturnValue::NotRunYet) {}

Raptor::~Raptor() {}

void Raptor::Clear() {
    results_.clear();
    results_.shrink_to_fit();
    mapping_status_ = RaptorReturnValue::NotRunYet;
}

RaptorReturnValue Raptor::Align(const mindex::SequenceFilePtr reads) {
    // Sanity check on parameters.
    if (params_ == nullptr) {
        return RaptorReturnValue::ParamsNotSet;
    }
    if (params_->num_threads <= 0) {
        return RaptorReturnValue::InvalidNumThreads;
    }
    if (index_ == nullptr) {
        return RaptorReturnValue::IndexIsNull;
    }

#ifdef RAPTOR_TESTING_MODE
        LOG_ALL("Pairwise aligning reads to the reference.\n");
#endif

    auto pairwise_ret = PairwiseAlign_(reads, results_, cleanup_target_hits_);

    mapping_status_ = RaptorReturnValue::OK;

    return RaptorReturnValue::OK;
}

RaptorReturnValue Raptor::PairwiseAlign_(const mindex::SequenceFilePtr reads, std::vector<std::unique_ptr<RaptorResults>>& results, bool cleanup_target_hits) const {

    results.clear();
    results.resize(reads->seqs().size());

    // Sanity check on parameters.
    if (params_ == nullptr) {
        return RaptorReturnValue::ParamsNotSet;
    }
    if (params_->num_threads <= 0) {
        return RaptorReturnValue::InvalidNumThreads;
    }
    if (index_ == nullptr) {
        return RaptorReturnValue::IndexIsNull;
    }

    // Fill out the threads with reads.
    std::shared_ptr<thread_pool::ThreadPool> thread_pool =
        thread_pool::createThreadPool(params_->num_threads);
    std::vector<std::future<int>> thread_futures;

    int64_t start_read = (params_->debug_qid <= 0) ? params_->start_read : params_->debug_qid;

    // Initialize mapping jobs as non-temporary objects.
    std::vector<MappingJob> mapping_jobs;
    if (params_->num_threads > 1) {
        // Check if only a range of reads is to be processed.
        // This is handled by a special if case to increase
        // efficiency if all reads are to be processed (no branching
        // in the inner loop).
        if (params_->num_reads_to_process <= 0 && params_->debug_qid < 0 && params_->debug_qname.empty() == true) {
            for (int i = 0; i < reads->seqs().size(); i++) {
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }

        } else if (params_->debug_qid >= 0) {
            for (int i = 0; i < reads->seqs().size(); i++) {
                int64_t id = reads->seqs()[i]->abs_id();
                if (id != params_->debug_qid) {
                    continue;
                }
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }

        } else if (params_->debug_qname.empty() == false) {
            for (int i = 0; i < reads->seqs().size(); i++) {
                const auto& qname = raptor::TrimToFirstWhiteSpace(reads->seqs()[i]->header());
                if (qname != params_->debug_qname) {
                    continue;
                }
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }

        } else {
            for (int i = 0; i < reads->seqs().size(); i++) {
                int64_t id = reads->seqs()[i]->abs_id();
                if (id < start_read) {
                    continue;
                }
                if (id >= (start_read + params_->num_reads_to_process)) {
                    break;
                }
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }
        }
    } else {
        if (params_->num_reads_to_process <= 0 && params_->debug_qid < 0 && params_->debug_qname.empty() == true) {
            MappingJob mj;
            mj.read_range.start = 0;
            mj.read_range.end = reads->seqs().size();
            mapping_jobs.push_back(mj);

        } else if (params_->debug_qid >= 0) {
            for (int i = 0; i < reads->seqs().size(); i++) {
                int64_t id = reads->seqs()[i]->abs_id();
                if (id != params_->debug_qid) {
                    continue;
                }
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }

        } else if (params_->debug_qname.empty() == false) {
            for (int i = 0; i < reads->seqs().size(); i++) {
                const auto& qname = raptor::TrimToFirstWhiteSpace(reads->seqs()[i]->header());
                if (qname != params_->debug_qname) {
                    continue;
                }
                MappingJob mj;
                mj.read_range.start = i;
                mj.read_range.end = i + 1;
                mapping_jobs.push_back(mj);
            }

        } else {
            // Here we expect a linear increase in sequence IDs.
            int64_t id = reads->seqs()[0]->abs_id();
            int64_t num_reads = reads->seqs().size();

            if (id < (start_read + params_->num_reads_to_process) &&
                (id + num_reads) >= start_read) {
                MappingJob mj;
                mj.read_range.start = (id >= start_read) ? 0 : (start_read - id);
                mj.read_range.end =
                    ((id + num_reads) < (start_read + params_->num_reads_to_process))
                        ? (num_reads)
                        : ((start_read + params_->num_reads_to_process) - id);
                mapping_jobs.push_back(mj);
            }
        }
    }

    // Add the jobs to the threadpool.
    for (int i = 0; i < mapping_jobs.size(); i++) {
        thread_futures.emplace_back(
            thread_pool->submit_task(MappingWorker_, reads, index_, graph_, ssg_, params_,
                                     std::cref(mapping_jobs[i]), std::ref(results), cleanup_target_hits));
    }

    // Wait for threads to finish.
    for (auto& it : thread_futures) {
        it.wait();
    }

    return RaptorReturnValue::OK;
}

int Raptor::MappingWorker_(const mindex::SequenceFilePtr reads, const mindex::IndexPtr index,
                           const raptor::GraphPtr graph, const raptor::SplitSegmentGraphPtr ssg,
                           const std::shared_ptr<ParamsRaptor> params,
                           const MappingJob& mapping_job,
                           std::vector<std::unique_ptr<RaptorResults>>& results,
                           bool cleanup_target_hits // If true, target_hits will be cleared to reduce memory consumption.
                           ) {
    // The results vector needs to be large enough to store results for a particular range.
    if (results.size() < mapping_job.read_range.end) {
        return 1;
    }

    auto aligner = raptor::createAligner(params->aligner_params->aligner_type, params->aligner_params->aligner_opt);
    auto aligner_ext = raptor::createAligner(params->aligner_params->aligner_type, params->aligner_params->aligner_opt);

    // raptor::AlignmentOptions aligner_opt_flank = params->aligner_params->aligner_opt;
    // aligner_opt_flank.
    // auto aligner_ext_flank = raptor::createAligner(params->aligner_params->aligner_type, params->aligner_params->aligner_opt);

    auto graph_mapper = createGraphMapper(index, graph, ssg, params->mapper_params);
    auto raptor_aligner = createRaptorAligner(index, params->aligner_params, aligner, aligner, aligner_ext);

    std::shared_ptr<raptor::LinearMappingResult> mapping_result;
    std::shared_ptr<raptor::GraphMappingResult> graph_mapping_result;
    std::shared_ptr<raptor::AlignedMappingResult> aln_result;

    for (int64_t i = mapping_job.read_range.start; i < mapping_job.read_range.end; i++) {
#ifdef RAPTOR_TESTING_MODE
        if (params->verbose_level >= 9) {
            LOG_ALL("Processing: q_id = %ld, q_len = %ld, qname = '%s'\n", i, reads->seqs()[i]->data().size(), reads->seqs()[i]->header().c_str());
        }
#endif
        if (reads->seqs()[i] == nullptr) {
            LOG_ALL("Warning: Nullptr sequence found: q_id = %ld, q_len = %ld, qname = '%s'. Skipping.\n", i, reads->seqs()[i]->data().size(), reads->seqs()[i]->header().c_str());
            continue;
        }

        // Shorthand to the qseq for later.
        const auto& qseq = reads->seqs()[i];

        // These are the results which will be stored in a new RaptorResults object.
        std::vector<std::shared_ptr<raptor::RegionBase>> regions;
        std::unordered_map<std::string, double> timings;
        int32_t mapq = 0;
        int64_t q_id = qseq->abs_id(); // It's important that this is not the absolute ID, but the local batch ID.

        // Linear mapping.
        auto mapper = createMapper(index, params->mapper_params);
        mapping_result = mapper->Map(qseq);

        if (graph != nullptr && params->mapper_params->no_graph_mapping == false) {
            graph_mapping_result = graph_mapper->Map(qseq, mapping_result);

        } else {
            // Map to a dummy graph, where each path will be (almost) equal to a single chain,
            // without any additional edges added.

            // mapping_result->Filter(params->bestn, params->bestn_threshold,
            // params->min_map_len, false);
            graph_mapping_result = graph_mapper->DummyMap(qseq, mapping_result);
        }

        if (params->do_align) {

            DEBUG_QSEQ(params, qseq, LOG_ALL("Aligning the paths because params->do_align == true.\n"));
            DEBUG_QSEQ(params, qseq, LOG_ALL("Filtering mappings...\n"));

            // Filter mappings, but leave room for error and mapq calculation.
            // The bestn value of 2 will allow for secondary alignments to come through if their score is close enough to the primary. Important because region with lower mapping score might align better.
            // In case params->bestn == 1, this will be applied after alignment.
            int64_t prefilter_bestn = (params->bestn <= 0) ? params->bestn : std::max(params->bestn, static_cast<int64_t>(2));
            graph_mapping_result->Filter(
                                            prefilter_bestn,
                                            params->bestn_threshold,
                                            params->min_map_len / 2,  // Ideally we wouldn't filter by min_map_len here because alignment can extend it, but realistically that can cause many short mappings to be aligned and thus slow down the process.
                                            false);

            DEBUG_QSEQ(params, qseq, LOG_ALL("Done filtering mappings.\n"));

            #ifdef RAPTOR_TESTING_MODE
                if (params->debug_qid == qseq->abs_id() ||
                    params->debug_qname == std::string(qseq->header())) {
                    if (graph_mapping_result == nullptr) {
                        LOG_ALL("Warning, graph_mapping_result == nullptr! q_id = %ld, q_len = %ld, qname = '%s'. Skipping.\n", i, qseq->data().size(), qseq->header().c_str());
                        continue;
                    }
                    const std::vector<std::shared_ptr<raptor::LocalPath>>& paths = graph_mapping_result->paths();
                    LOG_NOHEADER("Filtered paths before alignment:\n");
                    for (size_t path_id = 0; path_id < paths.size(); ++path_id) {
                        const auto path = paths[path_id];
                        if (path == nullptr) {
                            LOG_NOHEADER("  [path %ld] nullptr\n\n", path_id);
                            continue;
                        }
                        LOG_NOHEADER("  [path %ld] %s\n\n", path_id, path->Verbose().c_str());
                        // LOG_NOHEADER("  [path %ld]\n\n", path_id);
                    }
                }
            #endif

            DEBUG_QSEQ(params, qseq, LOG_ALL("Aligning paths...\n"));

            aln_result = raptor_aligner->AlignPaths(
                (reads->seqs()[i]), graph_mapping_result->paths());

            DEBUG_QSEQ(params, qseq, LOG_ALL("Finished AlignPaths. Now on to filtering.\n"));

            aln_result->Filter(params->bestn, params->bestn_threshold,
                                          params->min_map_len,
                                          params->min_identity,
                                          false);

            DEBUG_QSEQ(params, qseq, LOG_ALL("Done filtering.\n"));

            regions = aln_result->CollectRegions(params->one_hit_per_target);

            DEBUG_QSEQ(params, qseq, LOG_ALL("Done collecting regions.\n"));

        } else {
            DEBUG_QSEQ(params, qseq, LOG_ALL("Not aligning the paths because params->do_align == false.\n"));

            // Extend the flanks if required.
            if (params->mapper_params->flank_ext_len != 0) {
                graph_mapping_result->Filter(params->bestn,
                                                params->bestn_threshold,
                                                params->min_map_len / 2,  // Ideally we wouldn't filter by min_map_len here because alignment can extend it, but realistically that can cause many short mappings to be aligned and thus slow down the process.
                                                false);
                DEBUG_QSEQ(params, qseq, LOG_ALL("Done filtering (1/2).\n"));

                DEBUG_QSEQ(params, qseq, LOG_ALL("Extending flanks up to length %ld.\n", params->mapper_params->flank_ext_len));
                const std::vector<std::shared_ptr<raptor::LocalPath>>& paths = graph_mapping_result->paths();
                std::vector<std::shared_ptr<raptor::LocalPath>> new_paths;
                for (size_t path_id = 0; path_id < paths.size(); ++path_id) {
                    auto path = paths[path_id];
                    if (path == nullptr) {
                        continue;
                    }
                    auto new_path = raptor::PathAligner::FlankExtend(index, qseq, aligner_ext, params->mapper_params->flank_ext_len, path);
                    new_paths.emplace_back(new_path);
                }
                graph_mapping_result->paths(new_paths);

                // Filter again with the full min_map_len, now that the flanks were extended.
                graph_mapping_result->Filter(params->bestn,
                                                params->bestn_threshold,
                                                params->min_map_len,
                                                false);
                DEBUG_QSEQ(params, qseq, LOG_ALL("Done filtering (2/2).\n"));

            } else {
                graph_mapping_result->Filter(params->bestn,
                                                params->bestn_threshold,
                                                params->min_map_len,
                                                false);
                DEBUG_QSEQ(params, qseq, LOG_ALL("Done filtering (1/1).\n"));
                DEBUG_QSEQ(params, qseq, LOG_ALL("Not extending the flanks because flank_ext_len == 0."));
            }

            regions = graph_mapping_result->CollectRegions(params->one_hit_per_target);

            // If required, align using a diff aligner.
            if (params->do_diff) {
                DEBUG_QSEQ(params, qseq, LOG_ALL("Performing diff alignment without traceback."));

                auto difflib = raptor::createDifflibEdlib();
                for (size_t region_id = 0; region_id < regions.size(); ++region_id) {
                    auto& region = regions[region_id];
                    const char* tseq_raw = (const char*) index->FetchRawSeq(region->TargetID());

                    int64_t diffs = difflib->CalcDiffs((const char *) qseq->data().data(), region->QueryStart(), region->QueryEnd(), region->QueryLen(),
                                                        tseq_raw, region->TargetRev(), region->TargetFwdStart(), region->TargetFwdEnd(), region->TargetLen());

                    region->SetEditDistance(diffs);
                }
            }

            DEBUG_QSEQ(params, qseq, LOG_ALL("Done collecting regions.\n"));
        }

        // Accumulate all timing results.
        if (mapping_result != nullptr) {
            const auto& other_timings = mapping_result->Timings();
            timings.insert(other_timings.begin(), other_timings.end());
        }
        if (graph_mapping_result != nullptr) {
            const auto& other_timings = graph_mapping_result->Timings();
            timings.insert(other_timings.begin(), other_timings.end());
        }
        if (aln_result != nullptr) {
            const auto& other_timings = aln_result->Timings();
            timings.insert(other_timings.begin(), other_timings.end());
        }

        results[i] = raptor::createRaptorResults(q_id, regions, timings);
    }

    return 0;
}

} /* namespace raptor */
