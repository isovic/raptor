/*
 * raptor_sova_worker_queue.h
 *
 *  Created on: Nov 13, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RAPTOR_HOLY_H_
#define SRC_RAPTOR_RAPTOR_HOLY_H_

#include <stdint.h>
#include <memory>
#include <sequences/sequence_file.h>
#include <containers/range.h>
#include <raptor_sova/params_raptor_sova.h>
#include <containers/raptor_results.h>
#include <raptor/yield_index.h>
#include <types/typedefs.h>

namespace raptor {
namespace sova {

class RaptorSovaWorkerQueue;

std::unique_ptr<RaptorSovaWorkerQueue> createRaptorSovaWorkerQueue(const mindex::IndexPtr index,
                                     const std::shared_ptr<ParamsRaptorSova> params);

enum class HifiOverlapperReturnValue {
    OK,
    ParamsNotSet,
    InvalidNumThreads,
    IndexVectorEmpty,
    IndexIsNull,
    NotRunYet,
    Failed
};

typedef struct {
    Range read_range;
} MappingJob;

/* RaptorSovaWorkerQueue API takes an index of the reference sequences,
 * and user-specified parameters.
 * The object then allows mapping of a set of sequences to the index.
 */
class RaptorSovaWorkerQueue {
   public:
    friend std::unique_ptr<RaptorSovaWorkerQueue> createRaptorSovaWorkerQueue(const mindex::IndexPtr index,
                                                const std::shared_ptr<ParamsRaptorSova> params);
    ~RaptorSovaWorkerQueue();

    void Clear();
    HifiOverlapperReturnValue Align(const mindex::SequenceFilePtr reads);

    const std::vector<std::unique_ptr<RaptorResults>>& results() const { return results_; }

   private:
    RaptorSovaWorkerQueue(const mindex::IndexPtr index,
           const std::shared_ptr<ParamsRaptorSova> params);
    RaptorSovaWorkerQueue(const RaptorSovaWorkerQueue&) = delete;
    RaptorSovaWorkerQueue& operator=(const RaptorSovaWorkerQueue&) = delete;

    HifiOverlapperReturnValue PairwiseAlign_(const mindex::SequenceFilePtr reads, std::vector<std::unique_ptr<RaptorResults>>& results) const;

    static int MappingWorker_(const mindex::SequenceFilePtr reads, const mindex::IndexPtr index,
                              const std::shared_ptr<ParamsRaptorSova> params,
                              const MappingJob& mapping_job,
                              std::vector<std::unique_ptr<RaptorResults>>& results
                              );

    const std::shared_ptr<ParamsRaptorSova> params_;
    const mindex::IndexPtr index_;
    std::vector<std::unique_ptr<RaptorResults>> results_;
    HifiOverlapperReturnValue mapping_status_;
};

} /* namespace sova */
} /* namespace raptor */

#endif /* SRC_RAPTOR_RAPTOR_H_ */
