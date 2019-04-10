/*
 * anchor_graph_aligner.cc
 *
 *  Created on: Feb 1, 2019
 *      Author: Ivan Sovic
 */

#include <graph_aligner/anchor_graph_aligner.h>

namespace raptor {

std::unique_ptr<AnchorGraphAligner> createAnchorGraphAligner(
                        const raptor::SplitSegmentGraphPtr& _ss_graph,
                        const raptor::AlignmentOptions& _opt) {
    return std::unique_ptr<raptor::AnchorGraphAligner>(new raptor::AnchorGraphAligner(_ss_graph, _opt));
}

AnchorGraphAligner::AnchorGraphAligner(
            const raptor::SplitSegmentGraphPtr& _ss_graph,
            const raptor::AlignmentOptions& _opt)
            :   ss_graph_(_ss_graph),
                opt_(_opt)
{
}

void AnchorGraphAligner::Align(const raptor::AnchorGraphPtr& _anchor_graph) {

}

}
