/*
 * raptor.h
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RAPTOR_H_
#define SRC_RAPTOR_RAPTOR_H_

#include <stdint.h>
#include <memory>
#include <sequences/sequence_file.h>
#include <containers/range.h>
#include <params/params_raptor.h>
#include <containers/raptor_results.h>
#include <raptor/yield_index.h>
#include <graph/split_segment_graph.h>

namespace raptor {

class Raptor;

std::unique_ptr<Raptor> createRaptor(const mindex::IndexPtr index, const raptor::GraphPtr graph,
                                     const raptor::SplitSegmentGraphPtr ssg,
                                     const std::shared_ptr<ParamsRaptor> params,
                                     bool cleanup_target_hits);

enum class RaptorReturnValue {
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

/* Raptor API takes an index of the reference sequences,
 * and user-specified parameters.
 * The object then allows mapping of a set of sequences to the index.
 */
class Raptor {
   public:
    friend std::unique_ptr<Raptor> createRaptor(const mindex::IndexPtr index, const raptor::GraphPtr graph,
                                                const raptor::SplitSegmentGraphPtr ssg,
                                                const std::shared_ptr<ParamsRaptor> params,
                                                bool cleanup_target_hits);
    ~Raptor();

    void Clear();
    RaptorReturnValue Align(const mindex::SequenceFilePtr reads);

    const std::vector<RaptorResults>& results() const { return results_; }

   private:
    Raptor(const mindex::IndexPtr index, const raptor::GraphPtr graph,
           const raptor::SplitSegmentGraphPtr ssg,
           const std::shared_ptr<ParamsRaptor> params,
           bool cleanup_target_hits);
    Raptor(const Raptor&) = delete;
    Raptor& operator=(const Raptor&) = delete;

    RaptorReturnValue PairwiseAlign_(const mindex::SequenceFilePtr reads, std::vector<RaptorResults>& results, bool cleanup_target_hits) const;

    static int MappingWorker_(const mindex::SequenceFilePtr reads, const mindex::IndexPtr index,
                              const raptor::GraphPtr graph, const raptor::SplitSegmentGraphPtr ssg,
                              const std::shared_ptr<ParamsRaptor> params,
                              const MappingJob& mapping_job, std::vector<RaptorResults>& results,
                              bool cleanup_target_hits // If true, target_hits will be cleared to reduce memory consumption.
                              );

    const std::shared_ptr<ParamsRaptor> params_;
    bool cleanup_target_hits_;
    const mindex::IndexPtr index_;
    const raptor::GraphPtr graph_;
    const raptor::SplitSegmentGraphPtr ssg_;
    std::vector<RaptorResults> results_;
    RaptorReturnValue mapping_status_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RAPTOR_H_ */
