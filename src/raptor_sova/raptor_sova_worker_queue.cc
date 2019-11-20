/*
 * raptor_sova_worker_queue.cc
 *
 *  Created on: Nov 13, 2019
 *      Author: Ivan Sovic
 */

#include <thread>
#include <raptor_sova/raptor_sova_worker_queue.h>
#include <lib/thread_pool.hpp>
#include <log/log_tools.h>
#include <aligner/aligner_containers.h>
#include <aligner/aligner_factory.h>
#include <raptor_sova/sova_mapper.h>
#include <raptor/raptor_aligner.h>
#include <raptor/graph_mapper.h>
#include <raptor/path_aligner.h>
#include <utility/range_tools.hpp>
#include <utility/stringutil.h>
#include <aligner/difflib_edlib.h>

namespace raptor {
namespace sova {

std::unique_ptr<RaptorSovaWorkerQueue> createRaptorSovaWorkerQueue(const mindex::IndexPtr index,
                                     const std::shared_ptr<ParamsRaptorSova> params) {
    return std::unique_ptr<RaptorSovaWorkerQueue>(new RaptorSovaWorkerQueue(index, params));
}

RaptorSovaWorkerQueue::RaptorSovaWorkerQueue(const mindex::IndexPtr index,
               const std::shared_ptr<ParamsRaptorSova> params)
    : params_(params),
      index_(index),
      results_(),
      mapping_status_(HifiOverlapperReturnValue::NotRunYet) {}

RaptorSovaWorkerQueue::~RaptorSovaWorkerQueue() {}

void RaptorSovaWorkerQueue::Clear() {
    results_.clear();
    results_.shrink_to_fit();
    mapping_status_ = HifiOverlapperReturnValue::NotRunYet;
}

HifiOverlapperReturnValue RaptorSovaWorkerQueue::Align(const mindex::SequenceFilePtr reads) {
    // Sanity check on parameters.
    if (params_ == nullptr) {
        return HifiOverlapperReturnValue::ParamsNotSet;
    }
    if (params_->num_threads <= 0) {
        return HifiOverlapperReturnValue::InvalidNumThreads;
    }
    if (index_ == nullptr) {
        return HifiOverlapperReturnValue::IndexIsNull;
    }

    #ifdef RAPTOR_TESTING_MODE
            LOG_ALL("Pairwise aligning reads to the reference.\n");
    #endif

    auto pairwise_ret = PairwiseAlign_(reads, results_);

    mapping_status_ = HifiOverlapperReturnValue::OK;

    return HifiOverlapperReturnValue::OK;
}

HifiOverlapperReturnValue RaptorSovaWorkerQueue::PairwiseAlign_(const mindex::SequenceFilePtr reads, std::vector<std::unique_ptr<RaptorResults>>& results) const {

    results.clear();
    results.resize(reads->seqs().size());

    // Sanity check on parameters.
    if (params_ == nullptr) {
        return HifiOverlapperReturnValue::ParamsNotSet;
    }
    if (params_->num_threads <= 0) {
        return HifiOverlapperReturnValue::InvalidNumThreads;
    }
    if (index_ == nullptr) {
        return HifiOverlapperReturnValue::IndexIsNull;
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
            thread_pool->submit_task(MappingWorker_, reads, index_, params_,
                                     std::cref(mapping_jobs[i]), std::ref(results)));
    }

    // Wait for threads to finish.
    for (auto& it : thread_futures) {
        it.wait();
    }

    return HifiOverlapperReturnValue::OK;
}

int RaptorSovaWorkerQueue::MappingWorker_(const mindex::SequenceFilePtr reads, const mindex::IndexPtr index,
                           const std::shared_ptr<ParamsRaptorSova> params,
                           const MappingJob& mapping_job,
                           std::vector<std::unique_ptr<RaptorResults>>& results) {
    // The results vector needs to be large enough to store results for a particular range.
    if (results.size() < mapping_job.read_range.end) {
        return 1;
    }

    // Create the mapper.
    auto mapper = raptor::sova::createSovaMapper(index, params->mapper_params);

    // Process all queries in the batch.
    for (int64_t i = mapping_job.read_range.start; i < mapping_job.read_range.end; i++) {
        #ifdef RAPTOR_TESTING_MODE
                if (params->verbose_level >= 9) {
                    LOG_ALL("Processing: q_id = %ld, q_len = %ld, qname = '%s'\n",
                            i, reads->seqs()[i]->data().size(), reads->seqs()[i]->header().c_str());
                }
        #endif
        if (reads->seqs()[i] == nullptr) {
            LOG_ALL("Warning: Nullptr sequence found: q_id = %ld, q_len = %ld, qname = '%s'. Skipping.\n",
                            i, reads->seqs()[i]->data().size(), reads->seqs()[i]->header().c_str());
            continue;
        }

        // Map the query.
        const auto& qseq = reads->seqs()[i];
        std::shared_ptr<raptor::LinearMappingResult> mapping_result = mapper->Map(qseq);
        if (mapping_result == nullptr) {
            continue;
        }

        // Store the results.
        std::vector<std::shared_ptr<raptor::RegionBase>> regions =
                        mapping_result->CollectRegions(params->one_hit_per_target);
        int32_t mapq = 100;
        int64_t q_id = qseq->abs_id();
        const auto& timings = mapping_result->Timings();
        results[i] = raptor::createRaptorResults(q_id, regions, timings, mapq);
    }

    return 0;
}

} /* namespace sova */
} /* namespace raptor */
