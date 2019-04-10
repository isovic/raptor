/*
 * anchor_graph_aligner.h
 *
 *  Created on: Feb 01, 2019
 *      Author: isovic
 */

#ifndef SRC_ANCHOR_GRAPH_ALIGNER_GRAPH_ALIGNER_H_
#define SRC_ANCHOR_GRAPH_ALIGNER_GRAPH_ALIGNER_H_

#include <memory>
#include <vector>
#include <aligner/aligner_containers.h>
#include <aligner/pairwise_penalties.h>
#include <aligner/aligner_util.hpp>
#include <graph/split_segment_graph.h>
#include <graph/anchor_graph.h>
#include <index/sequence.h>

namespace raptor {

class AnchorGraphAligner;

std::unique_ptr<AnchorGraphAligner> createAnchorGraphAligner(
                        const raptor::SplitSegmentGraphPtr& _ss_graph,
                        const raptor::AlignmentOptions& _opt);

class AnchorGraphAligner {
   public:
    friend std::unique_ptr<AnchorGraphAligner> createAnchorGraphAligner(
                            const raptor::SplitSegmentGraphPtr& _ss_graph,
                            const raptor::AlignmentOptions& _opt);

    ~AnchorGraphAligner();

    void Align(const raptor::AnchorGraphPtr& _anchor_graph);

    // std::shared_ptr<raptor::AlignmentResult> Global(const char* qseq, int64_t qlen,
    //                                                 const char* tseq,
    //                                                 int64_t tlen);  // Global alignment mode.

    // std::shared_ptr<raptor::AlignmentResult> Local(const char* qseq, int64_t qlen, const char* tseq,
    //                                                int64_t tlen);  // Local alignment mode.

    // std::shared_ptr<raptor::AlignmentResult> Semiglobal(
    //     const char* qseq, int64_t qlen, const char* tseq,
    //     int64_t tlen);  // Semiglobal alignment mode.

   private:
    AnchorGraphAligner(
                const raptor::SplitSegmentGraphPtr& _ss_graph,
                const raptor::AlignmentOptions& _opt);

    AnchorGraphAligner(const AnchorGraphAligner&) = delete;              // No copying.
    AnchorGraphAligner& operator=(const AnchorGraphAligner&) = delete;   // No copying.
    AnchorGraphAligner(AnchorGraphAligner&&) = delete;                   // No move constructor.
    AnchorGraphAligner& operator=(const AnchorGraphAligner&&) = delete;  // No copying.

    // static raptor::SplitSegmentGraphPtr ExtractSubgraphFromSSG_(const raptor::SplitSegmentGraphPtr& ssg,
    //                     const std::unordered_map<int64_t, int8_t>& selected_nodes,
    //                     const std::unordered_map<int64_t, int8_t>& selected_edges);

    const raptor::SplitSegmentGraphPtr& ss_graph_;
    const raptor::AlignmentOptions& opt_;
};

} /* namespace raptor */

#endif /* SRC_ANCHOR_GRAPH_ALIGNER_GRAPH_ALIGNER_H_ */