/*
 * graph_mapper.h
 *
 *  Created on: Dec 17, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_GRAPH_MAPPER_H_
#define SRC_RAPTOR_GRAPH_MAPPER_H_

#include <cstdint>
#include <memory>
#include <vector>
#include <unordered_map>
#include <sequences/sequence_file.h>
#include <sequences/sequence.h>
#include <containers/target_hits.hpp>
#include <containers/region/region_mapped.h>
#include <containers/mapping_result/linear_mapping_result.h>
#include <raptor/yield_index.h>
#include <params/params_mapper.h>
#include <index/minimizer_index.h>
#include <containers/mapping_result/graph_mapping_result.h>
#include <types/typedefs.h>
#include <types/typedefs_interval_tree.h>
#include <graph/anchor_graph_edge.h>
#include <algorithm/path.hpp>
#include <graph/anchor_graph.h>
#include <graph/local_path.h>
#include <raptor/backtrack_list_node.h>
#include <graph/split_segment_graph.h>

namespace raptor {

class GraphMapper;

/*
  If an extension aligner is not provided, only heuristic
  filtering will be used.
*/
std::unique_ptr<raptor::GraphMapper> createGraphMapper(
    const mindex::IndexPtr index, const raptor::GraphPtr graph,
    const raptor::SplitSegmentGraphPtr& ssg,
    const std::shared_ptr<raptor::ParamsMapper> params);

/* GraphMapper object is intended for mapping of one query
 * sequence to a sequence graph.
 */
class GraphMapper {
   public:
    friend std::unique_ptr<raptor::GraphMapper> createGraphMapper(
        const mindex::IndexPtr index, const raptor::GraphPtr graph,
        const raptor::SplitSegmentGraphPtr& ssg,
        const std::shared_ptr<raptor::ParamsMapper> params);
    ~GraphMapper();

    std::shared_ptr<raptor::GraphMappingResult> Map(
        const mindex::SequencePtr& qseq, const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result);

    std::shared_ptr<raptor::GraphMappingResult> MapBetter(
        const mindex::SequencePtr& qseq, const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result);

    /*
     * Simply translate each input LinearMappingResult as an individual node, without
     * constructing any additional edges.
     */
    std::shared_ptr<raptor::GraphMappingResult> DummyMap(
        const mindex::SequencePtr& qseq, const std::shared_ptr<raptor::LinearMappingResult> input_mapping_result);


    static std::shared_ptr<raptor::AnchorGraph> ConstructAnchorGraph(
        const raptor::GraphPtr graph,
        const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors,
        int32_t max_allowed_dist,
        int32_t predecessor_lookup_window, int32_t max_path_edges,
        int32_t chain_max_skip,
        bool verbose_debug_qid);

    static std::shared_ptr<raptor::AnchorGraph> GraphDP2(
        std::shared_ptr<raptor::AnchorGraph>& local_graph,
        int32_t chain_max_bandwidth,
        bool verbose_debug_qid);

    // /*
    //  * ReduceMappedAnchorGraph edits the graph in place.
    // */
//    static void ReduceMappedAnchorGraph(
//                                 std::shared_ptr<raptor::AnchorGraph>& graph,
//                                 double allowed_score_diff_frac,
//                                 bool verbose_debug_qid);

    static std::vector<raptor::AnchorGraphPtr> BacktrackReducedMappedAnchorGraph(
                                            const std::shared_ptr<raptor::AnchorGraph>& graph,
                                            const double allowed_score_diff_frac,
                                            std::unordered_map<int64_t, int64_t>& node_to_path_score,
                                            const bool verbose_debug_qid);

    static std::vector<std::tuple<raptor::AnchorGraphPtr, int64_t, double>> CreateAnchorGraphConnectedComponents(
                                                const std::shared_ptr<raptor::AnchorGraph>& graph,
                                                const double allowed_score_diff_frac,
                                                const bool verbose_debug_qid);


    static raptor::SplitSegmentGraphPtr ConstructAlignmentGraph(const raptor::SplitSegmentGraphPtr& ss_graph,
                                        const raptor::AnchorGraphPtr& anchor_graph,
                                        const mindex::SequencePtr& qseq,
                                        bool use_approx_aln,
                                        bool verbose_debug_qid);

   private:
    GraphMapper(const mindex::IndexPtr index, const raptor::GraphPtr graph,
                const raptor::SplitSegmentGraphPtr& ssg,
                const std::shared_ptr<raptor::ParamsMapper> params);

    GraphMapper(const GraphMapper&) = delete;
    GraphMapper& operator=(const GraphMapper&) = delete;

    void AlignAnchors_(const mindex::SequencePtr& qseq,
                        std::vector<raptor::TargetAnchorPtr>& target_anchors) const;

    std::vector<std::shared_ptr<raptor::LocalPath>> GraphToPaths2_(
                                            const std::shared_ptr<raptor::AnchorGraph>& graph,
                                            const std::unordered_map<int64_t, int64_t>& node_to_path_score,
                                            const bool verbose_debug_qid);

    static bool FindLinearSplitSegmentNodesAndEdges_(const raptor::SplitSegmentGraphPtr& ssg,
                        int64_t t_id, bool t_rev, int64_t start, int64_t end,
                        std::vector<raptor::SplitSegmentGraphNameType>& ret_nodes,
                        std::vector<raptor::SplitSegmentGraphNameType>& ret_edges);

    const mindex::IndexPtr index_;
    const raptor::GraphPtr graph_;
    const raptor::SplitSegmentGraphPtr ssg_;
    const std::shared_ptr<raptor::ParamsMapper> params_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_GRAPH_MAPPER_H_ */
