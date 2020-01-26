/*
 * path_aligner.cc
 *
 *  Created on: Dec 29, 2017
 *      Author: Ivan Sovic
 */

#include <raptor/path_aligner.h>
#include <aligner/aligner_util.hpp>
#include <aligner/cigar.h>
#include <algorithm>
#include <log/log_tools.h>
#include <utility/revcmp.hpp>
#include <utility/tictoc.h>

namespace raptor {

std::shared_ptr<raptor::PathAligner> createPathAligner(const mindex::IndexPtr _index,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                                                   std::shared_ptr<raptor::AlignerBase> _aligner_ext) {
    return std::shared_ptr<raptor::PathAligner>(
        new raptor::PathAligner(_index, _aligner, _aligner_gap, _aligner_ext));
}

PathAligner::PathAligner(const mindex::IndexPtr _index, std::shared_ptr<raptor::AlignerBase> _aligner,
                         std::shared_ptr<raptor::AlignerBase> _aligner_gap,
                         std::shared_ptr<raptor::AlignerBase> _aligner_ext)
    : index_(_index), aligner_(_aligner), aligner_gap_(_aligner_gap), aligner_ext_(_aligner_ext) {}

PathAligner::~PathAligner() {}

/*
 * This method expects valid inputs. This means: no nullptrs, and all neighboring segment edges need
 * to be on the same seqs. It composes a sequence between multiple consecutive segment edges.
 */
bool PathAligner::ComposeEdgeSeq_(const mindex::SequencePtr& qseq,
                                  std::shared_ptr<raptor::AnchorGraphEdge> local_edge,
                                  std::string& edge_seq,
                                  std::vector<AlignmentRegion>& edge_seq_regions) const {
    edge_seq = std::string();
    edge_seq_regions.clear();

//    LOG_ALL("local_edge = %s\n", local_edge->Verbose().c_str());
//    LOG_ALL("local_edge.use_count() = %ld\n", local_edge.use_count());
//    fflush(stderr);
//    fflush(stdout);

    for (size_t i = 1; i < local_edge->segment_edges().size(); i++) {
        auto prev_seg_edge = local_edge->segment_edges()[i - 1];
        auto curr_seg_edge = local_edge->segment_edges()[i];

        if (prev_seg_edge == nullptr) {
            LOG_ALL("WARNING: prev_seg_edge == nullptr!\n");
        }
        if (curr_seg_edge == nullptr) {
            LOG_ALL("WARNING: curr_seg_edge == nullptr!\n");
        }

        if (prev_seg_edge == nullptr || curr_seg_edge == nullptr) {
            LOG_ALL("Something funky with the seg_edges in ComposeEdgeSeq_! Skipping.\n");
// #           LOG_ALL("local_edge = %s\n", local_edge->Verbose().c_str());
            return false;
        }

        int32_t rid = curr_seg_edge->source_id();
        int32_t rstart = prev_seg_edge->sink_end();
        int32_t rend = curr_seg_edge->source_start();
        int32_t rstart_orig = rstart;
        int32_t rend_orig = rend;
        bool rrev = curr_seg_edge->source_is_rev();

        if ((rstart - rend) == 0) {
            LOG_ALL("Note: internal edge sequence has 0 length. (Not an error, just an informative remark.)\n");
            continue;
        }

        if (rrev) {
            std::swap(rstart, rend);
            rstart = curr_seg_edge->source_len() - rstart;
            rend = curr_seg_edge->source_len() - rend;
        }

        // printf ("rid = %ld, index_->starts().size() = %ld, index_->lens() = %ld\n", rid,
        // index_->starts().size(), index_->lens().size()); fflush(stdout);

        auto env = raptor::createMappingEnv(rid, 0, index_->len(rid), rrev,
                                        qseq->abs_id(), qseq->data().size(),
                                        false);

        edge_seq_regions.emplace_back(AlignmentRegion(rstart_orig, rend_orig, edge_seq.size(),
                                                      edge_seq.size() + (rend - rstart), env));

        edge_seq += index_->FetchSeqAsString(rid, rstart, rend, rrev);
    }

    return true;
}

/*
 * local_edge can be a nullptr, for example in case when the implicit_last_node is
 * actually the last node in the path, and there are no edges/nodes following it.
 * Otherwise, it points to an edge from which the coordinate will be taken.
 */
bool PathAligner::ComposeImplicitSeq_(
    const mindex::SequencePtr& qseq, std::shared_ptr<raptor::LocalPath> path,
    std::shared_ptr<raptor::AnchorGraphEdge> local_edge,  // The non-implicit edge connecting nodes.
    int32_t implicit_first_node,  // First node of a streak of consecutive implicit nodes.
    int32_t implicit_last_node,   // Last node of a streak of consecutive implicit nodes. Should be
                                  // the same as the source of local_edge.
    std::string& node_seq,        // Return target sequence for the selected region.
    AlignmentRegion& node_seq_region) const {  // Return region of the selected sequence.
    /*
     * Handle nodes linked by implicit edges.
     * In case there are no implicit edges, this chunk will output
     * the sequence covered by the source node of the edge.
     */
    auto& first_node = path->nodes()[implicit_first_node]->data();
    auto& last_node = path->nodes()[implicit_last_node]->data();

    // If the first node of an implicit block actually had a proper non-implicit in-edge,
    // we need to use that edge's start location to compensate for potential missing sequence
    // when mapping.
    int32_t rstart = (implicit_first_node == 0)
                         ? (first_node->TargetStart())
                         : (path->edges()[implicit_first_node - 1]->data()->sink_start());
    int32_t rend =
        (local_edge == nullptr)
            ? last_node->TargetEnd()
            : local_edge->source_start();  // Edge can start within the node, or after the node.
    int32_t rstart_orig = rstart;
    int32_t rend_orig = rend;

    // In this case, the implicit sequences are fully covered by the edge.
    // This concrete region will be skipped, and the sequence covered by the
    // edge will be taken from the other side of the edge.
    if (rend < rstart) {
        return false;
    }

    // Handle reverse complements.
    if (first_node->TargetRev()) {
        std::swap(rstart, rend);
        rstart = first_node->TargetLen() - rstart;
        rend = first_node->TargetLen() - rend;
    }

    node_seq =
        index_->FetchSeqAsString(first_node->TargetID(), rstart, rend, first_node->TargetRev());

    auto env =
        raptor::createMappingEnv(first_node->TargetID(), 0,
                             index_->len(first_node->TargetID()), first_node->TargetRev(),
                             qseq->abs_id(), qseq->data().size(), false);

    // node_seq_region = AlignmentRegion(0, (rend - rstart), env);

    node_seq_region = AlignmentRegion(rstart_orig, rend_orig, 0, (rend - rstart), env);

    return true;
}

bool PathAligner::ComposeTargetSequence_(const mindex::SequencePtr& qseq,
                                         std::shared_ptr<raptor::LocalPath> path, std::string& ret_seq,
                                         std::vector<AlignmentRegion>& regions) const {
    ret_seq = std::string();
    regions.clear();

    if (path->nodes().size() == 0) {
        return true;
    }

    int32_t implicit_first_node = 0, implicit_last_node = 0;
    for (size_t i = 0; i < path->edges().size(); i++) {
        auto& edge = path->edges()[i]->data();
        implicit_last_node = i;

        // Implicit edges will simply be joined into a single
        // start/end block.
        if (edge->IsImplicit()) {
            continue;
        }

        // Get the implicit sequence.
        std::string node_seq;
        AlignmentRegion node_seq_region(0, 0, 0, 0, nullptr);
        auto rv_implicit = ComposeImplicitSeq_(qseq, path, edge, implicit_first_node,
                                               implicit_last_node, node_seq, node_seq_region);
        if (rv_implicit) {
            node_seq_region.spl_start_pos += ret_seq.size();
            node_seq_region.spl_end_pos += ret_seq.size();
            regions.emplace_back(node_seq_region);
            ret_seq += node_seq;
        }

        // Get the sequence spelled by consecutive edges.
        std::string edge_seq;
        std::vector<AlignmentRegion> edge_seq_regions;
        auto rv_edge = ComposeEdgeSeq_(qseq, edge, edge_seq, edge_seq_regions);
        if (rv_edge) {
            for (size_t reg_id = 0; reg_id < edge_seq_regions.size(); reg_id++) {
                edge_seq_regions[reg_id].spl_start_pos += ret_seq.size();
                edge_seq_regions[reg_id].spl_end_pos += ret_seq.size();
            }
            regions.insert(regions.end(), edge_seq_regions.begin(), edge_seq_regions.end());
            ret_seq += edge_seq;
        }

        // Update.
        implicit_first_node = i + 1;
    }

    {  // Add the last implicit chain of nodes.
        implicit_last_node = path->nodes().size() - 1;
        std::string node_seq;
        AlignmentRegion node_seq_region(0, 0, 0, 0, nullptr);
        auto rv_implicit = ComposeImplicitSeq_(qseq, path, nullptr, implicit_first_node,
                                               implicit_last_node, node_seq, node_seq_region);
        node_seq_region.spl_start_pos += ret_seq.size();
        node_seq_region.spl_end_pos += ret_seq.size();
        regions.emplace_back(node_seq_region);
        ret_seq += node_seq;
    }

    // std::cout << "AlignmentRegions:" << std::endl;
    // for (size_t i = 0; i < regions.size(); i++) {
    //     std::cout << "[i = " << i   << "] ref_start_pos = " << regions[i].ref_start_pos
    //                                 << ", ref_end_pos = " << regions[i].ref_end_pos
    //                                 << ", spl_start_pos = " << regions[i].spl_start_pos
    //                                 << ", spl_end_pos = " << regions[i].spl_end_pos
    //                                 << ", " << regions[i].env->Verbose().c_str() << std::endl;
    //     fflush(stdout);
    // }

    return true;
}

std::vector<std::shared_ptr<raptor::RegionAligned>> PathAligner::SplitAlignment_(
    const mindex::SequencePtr& qseq, const std::shared_ptr<raptor::LocalPath> path,
    const std::string& ref_seq, const std::vector<AlignmentRegion>& ref_seq_regions,
    const std::shared_ptr<raptor::AlignmentResult> entire_alignment,
    double min_idt) {
    std::vector<std::shared_ptr<raptor::RegionAligned>> results;

    // Query should be non-empty.
    if (qseq->data().size() == 0) {
        return results;
    }

    if (ref_seq.size() == 0) {
        return results;
    }

    // There should be at least one region for a valid input.
    if (ref_seq_regions.size() == 0) {
        return results;
    }

    // Path should have at least one node.
    if (path->nodes().size() == 0) {
        return results;
    }

    int64_t qstart = path->nodes().front()->data()->QueryStart();
    int64_t qend = path->nodes().back()->data()->QueryEnd();

    int64_t qpos = qstart;
    int64_t rpos = 0;

    auto aln_array = CigarToAlignmentArray(entire_alignment->cigar());
    size_t aln_id = 0;
    int32_t reg_id = 0;  // Simulate a stack.
    int64_t region_op_start = 0, region_op_end = 0;
    int64_t region_qstart = 0, region_qend = 0;

    // A temporary storage for: region_op_start, region_op_end, region_qstart, region_qend;
    std::vector<std::tuple<int64_t, int64_t, int64_t, int64_t>> region_op_pos;
    while (aln_id <= aln_array.size() && reg_id < ref_seq_regions.size()) {
        if (rpos == ref_seq_regions[reg_id].spl_start_pos) {
            region_op_start = aln_id;
            region_qstart = qpos;
        }

        if (rpos == ref_seq_regions[reg_id].spl_end_pos) {
            region_op_end = aln_id;
            region_qend = qpos;
            region_op_pos.emplace_back(
                std::make_tuple(region_op_start, region_op_end, region_qstart, region_qend));
            region_op_start = region_op_end = region_qstart = region_qend = 01;
            reg_id += 1;
            // Skip the increasing of rpos and qpos because the end position
            // of the current region is the same as the start position of
            // the next region.
            continue;
        }

        if (aln_id < aln_array.size() && IsCigarRef(aln_array[aln_id])) {
            rpos += 1;
        }
        if (aln_id < aln_array.size() && IsCigarQuery(aln_array[aln_id])) {
            qpos += 1;
        }

        aln_id += 1;
    }

    // At this point, there should be the same number of region_op_pos tuples as regions.
    if (region_op_pos.size() != ref_seq_regions.size()) {
        LOG_ALL("Failure to properly split regions for path:\n%s\n\n", path->Verbose().c_str());
        for (size_t i = 0; i < ref_seq_regions.size(); i++) {
            LOG_NOHEADER(
                "[ref_seq_regions %ld] ref_start_pos = %ld, ref_end_pos = %ld, spl_start_pos = "
                "%ld, spl_end_pos = %ld\n",
                i, ref_seq_regions[i].ref_start_pos, ref_seq_regions[i].ref_end_pos,
                ref_seq_regions[i].spl_start_pos, ref_seq_regions[i].spl_end_pos);
        }

        LOG_ALL("Wrong number of elements. region_op_pos.size() = %lu, ref_seq_regions.size() "
                "= %lu, qseq_id = %ld\n",
                region_op_pos.size(), ref_seq_regions.size(), qseq->abs_id());

        return results;

    }

    const auto aln_opt = aligner_->GetAlignmentOptions();

    for (size_t i = 0; i < ref_seq_regions.size(); i++) {
        auto aln = raptor::createAlignmentResult();
        auto result = createAlignedRegion(ref_seq_regions[i].env, aln, -1, -1, -1, -1);

        int64_t region_op_start = std::get<0>(region_op_pos[i]);
        int64_t region_op_end = std::get<1>(region_op_pos[i]);
        int64_t region_qstart = std::get<2>(region_op_pos[i]);
        int64_t region_qend = std::get<3>(region_op_pos[i]);

        auto new_cigar = AlignmentArrayToCigar((unsigned char*)(&aln_array[0] + region_op_start),
                                           (region_op_end - region_op_start));
        int64_t new_score = ScoreCigarAlignment(new_cigar, aln_opt.p.match, aln_opt.p.mismatch, aln_opt.p.w[0].open, aln_opt.p.w[0].ext);
        int64_t new_edit_dist = EditDistFromExtCIGAR(new_cigar);

        aln->score(new_score);
        aln->edit_dist(new_edit_dist);
        // Internally store the original target coordinates (on the strand the original mapping was
        // made to).
        aln->position(raptor::AlignmentPosition(region_qstart, region_qend, ref_seq_regions[i].ref_start_pos,
                                  ref_seq_regions[i].ref_end_pos));
        aln->final_band(-1);
        aln->status(raptor::AlignmentReturnValue::OK);

        // Add the soft clippings at front and back of the alignment.
        if (region_qstart > 0) {
            new_cigar.insert(new_cigar.begin(), raptor::CigarOp('S', region_qstart));
        }
        if ((qseq->data().size() - region_qend) > 0) {
            new_cigar.insert(new_cigar.end(),
                              raptor::CigarOp('S', (qseq->data().size() - region_qend)));
        }

        // Reverse the CIGAR if needed. This makes it ready for output.
        // TODO: Think about whether this needs to be reversed here, or is it
        // better to reverse it right before output?
        if (ref_seq_regions[i].env->t_rev) {
            std::reverse(new_cigar.begin(), new_cigar.end());
        }

        aln->cigar(new_cigar);

        if (aln->op_counts().identity_min >= min_idt) {
            // fprintf (stderr, "(Split) %s; aln->op_counts().eq = %d, aln->op_counts().aligned_qlen = %d\n", aln->Verbose().c_str(), aln->op_counts().eq, aln->op_counts().aligned_qlen);
            results.emplace_back(result);
        }
    }

    // // ther is ref_seq_regions.
    // printf ("region_op_pos.size() = %ld, ref_seq_regions.size() = %ld\n", region_op_pos.size(),
    // ref_seq_regions.size()); for (size_t i = 0; i < region_op_pos.size(); i++) {
    //     printf ("[%ld] region_op_start = %ld, region_op_end = %ld, region_qstart = %ld,
    //     region_qend = %ld\n",
    //             i, std::get<0>(region_op_pos[i]), std::get<1>(region_op_pos[i]),
    //             std::get<2>(region_op_pos[i]), std::get<3>(region_op_pos[i]));
    // }

    // printf ("\n");
    // printf ("Results:\n");
    // for (size_t i = 0; i < results.size(); i++) {
    //     printf ("[%3ld] Env: %s\n      Aln: %s\n", i, results[i]->env()->Verbose().c_str(),
    //     results[i]->aln()->Verbose().c_str());
    // }

    // for (int32_t cig_id = 0; cig_id < entire_alignment->cigar.size(); cig_id++) {
    //     if (reg_id >=ref_seq_regions.size()) {
    //         break;
    //     }

    // }

    // int32_t t_id = 0, index_t_start = 0, t_len = 0;
    // bool t_rev = false;
    // int32_t q_id = 0, q_len = 0;

    // auto env = raptor::createMappingEnv(t_id, index_t_start,
    //                                     t_len, t_rev,
    //                                     q_id, q_len, false);
    // auto aln = raptor::createAlignmentResult();
    // auto result = createAlignedRegion(env, aln);
    // results.emplace_back(result);

    // aln->score = entire_alignment->score;
    // aln->edit_dist = EditDistFromExtCIGAR(aln->cigar);
    // // Internally store the original target coordinates (on the strand the original mapping was
    // made to). aln->position = raptor::AlignmentPosition(qstart, qend,
    // path->nodes().front()->data->TargetStart(), path->nodes().back()->data->TargetEnd()); aln->k
    // = -1; aln->status = raptor::AlignmentReturnValue::OK;
    // // Set the CIGAR.
    // aln->cigar = entire_alignment->cigar;
    // // Add the soft clippings at front and back.
    // if (qstart > 0) {
    //     aln->cigar.insert(aln->cigar.begin(), raptor::CigarOp('S', qstart));
    // }
    // if ((qseq.get_sequence_length() - qend) > 0) {
    //     aln->cigar.insert(aln->cigar.end(), raptor::CigarOp('S', (qseq.get_sequence_length() -
    //     qend)));
    // }
    // // Reverse the CIGAR if needed.
    // if (path->nodes().front()->data->TargetRev()) {
    //     std::reverse(aln->cigar.begin(), aln->cigar.end());
    // }

    return results;
}

std::shared_ptr<raptor::RegionMapped> PathAligner::ExtendAlignmentFront_(
                                        const mindex::IndexPtr index,
                                        const mindex::SequencePtr& qseq,
                                        const std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                                        int32_t max_flank_len,  // If < 0, unrestricted.
                                        const std::shared_ptr<raptor::RegionMapped> anchor) {

    auto new_anchor = raptor::createRegionMapped(anchor->id(), anchor->target_hits_id(), anchor->env(),
                                anchor->QueryStart(), anchor->QueryEnd(),
                                anchor->TargetStart(), anchor->TargetEnd(),
                                anchor->CoveredBasesQuery(), anchor->CoveredBasesTarget(),
                                anchor->NumSeeds(), anchor->EditDistance(),
                                anchor->Score(),
                                anchor->PathId(), anchor->PathsNum(),
                                anchor->SegmentId(), anchor->SegmentsNum());

    const char* query = (const char*)&(qseq->data()[0]);
    int32_t qstart = 0;
    int32_t qend = anchor->QueryStart();
    int32_t qspan = qend - qstart;

    int64_t rstart = 0;
    int64_t rend = anchor->TargetStart();
    int32_t rspan = rend - rstart;

    int32_t min_span = std::min(qspan, rspan);

    // rstart = std::max((int64_t) 0, (int64_t) (anchor->TargetStart() - qspan * 2));

    if (anchor->TargetRev()) {
        // Get the full-length flank.
        rstart = anchor->TargetLen() - (anchor->TargetStart());
        rend = anchor->TargetLen();
        rspan = rend - rstart;

        // Find the minimum size flank (either query or target).
        min_span = std::min(qspan, rspan);

        // Update the coordinates.
        rend = std::min((int64_t) anchor->TargetLen(), rstart + min_span * 2);
        qstart = std::max((int64_t) 0, (int64_t) (qend - min_span * 2));
        // rstart = anchor->TargetLen() - (anchor->TargetStart());
        // rend = std::min((int64_t) anchor->TargetLen(), rstart + qspan * 2);
    } else {
        qstart = std::max((int64_t) 0, (int64_t) (qend - min_span * 2));
        rstart = std::max((int64_t) 0, (int64_t) (rend - min_span * 2));
    }

    qspan = qend - qstart;
    rspan = rend - rstart;

    // If there is a restriction on the maximum flank length, and the flank is longer,
    // then do not extend but just return.
    if (max_flank_len >= 0 && std::min(qspan, rspan) > max_flank_len) {
        return new_anchor;
    }

    auto rev_query_front = raptor::ReverseComplement((const char*) &(qseq->data()[qstart]), qspan);
    std::string ref_seq = index->FetchSeqAsString(anchor->TargetID(), rstart, rend, !anchor->TargetRev());

    auto aln = aligner_ext->Extend(rev_query_front.c_str(), rev_query_front.size(), ref_seq.c_str(), ref_seq.size());

    // std::cerr << "qstart = " << qstart << std::endl;
    // std::cerr << "rstart = " << rstart << std::endl;
    // std::cerr << "rend = " << rend << std::endl;
    // std::cerr << anchor->Verbose() << std::endl;
    // std::cerr << "Query:\n" << rev_query_front << "\n";
    // std::cerr << "Target:\n" << ref_seq << "\n";
    // std::cerr << "Extend alignment result:\n" << aln->Verbose() << std::endl;

    if (aln->position().qend >= 0 && aln->position().tend >= 0) {
        new_anchor->qstart(anchor->qstart() - (aln->position().qend));
        new_anchor->tstart(anchor->tstart() - (aln->position().tend));
    }

    return new_anchor;
}

std::shared_ptr<raptor::RegionMapped> PathAligner::ExtendAlignmentBack_(
                                        const mindex::IndexPtr index,
                                        const mindex::SequencePtr& qseq,
                                        const std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                                        int32_t max_flank_len,  // If < 0, unrestricted.
                                        const std::shared_ptr<raptor::RegionMapped> anchor) {
    auto new_anchor = raptor::createRegionMapped(anchor->id(), anchor->target_hits_id(), anchor->env(),
                                anchor->QueryStart(), anchor->QueryEnd(),
                                anchor->TargetStart(), anchor->TargetEnd(),
                                anchor->CoveredBasesQuery(), anchor->CoveredBasesTarget(),
                                anchor->NumSeeds(), anchor->EditDistance(), anchor->Score(),
                                anchor->PathId(), anchor->PathsNum(),
                                anchor->SegmentId(), anchor->SegmentsNum());

    // auto ret = raptor::NodeItem<int64_t, raptor::RegionMapped>::createNodeItem(node->internal_id(), node->name(), new_anchor);

    const char* query = (const char*) &(qseq->data()[0]);
    int32_t qstart = anchor->QueryEnd();
    int32_t qend = anchor->QueryLen(); // - qstart;
    int32_t qspan = qend - qstart;

    int64_t rstart = anchor->TargetEnd();
    int64_t rend = anchor->TargetLen();
    int32_t rspan = rend - rstart;
    // int64_t rend = rstart + std::min(qspan * 2, anchor->TargetLen() - anchor->TargetEnd());

    int32_t min_span = std::min(qspan, rspan);

    if (anchor->TargetRev()) {
        rend = anchor->TargetLen() - anchor->TargetEnd();
        rstart = 0;
        // rstart = std::max((int64_t) 0, rend - min_span * 2);
        rspan = rend - rstart;

        min_span = std::min(qspan, rspan);

        rstart = std::max((int64_t) 0, rend - min_span * 2);
        qend = qstart + std::min(min_span * 2, anchor->QueryLen() - anchor->QueryEnd());

        // rend = anchor->TargetLen() - anchor->TargetEnd();
        // rstart = std::max((int64_t) 0, rend - qspan * 2);
    } else {
        qend = qstart + std::min(min_span * 2, anchor->QueryLen() - anchor->QueryEnd());
        rend = rstart + std::min(min_span * 2, anchor->TargetLen() - anchor->TargetEnd());
    }

    qspan = qend - qstart;
    rspan = rend - rstart;

    // If there is a restriction on the maximum flank length, and the flank is longer,
    // then do not extend but just return.
    if (max_flank_len >= 0 && std::min(qspan, rspan) > max_flank_len) {
        return new_anchor;
    }

    std::string ref_seq = index->FetchSeqAsString(anchor->TargetID(), rstart, rend, anchor->TargetRev());

    auto aln = aligner_ext->Extend(query + qstart, qend - qstart, ref_seq.c_str(), ref_seq.size());

    // std::cerr << "qstart = " << qstart << std::endl;
    // std::cerr << "rstart = " << rstart << std::endl;
    // std::cerr << "rend = " << rend << std::endl;
    // std::cerr << anchor->Verbose() << std::endl;
    // std::cerr << "Query:\n" << std::string((query + qstart), qspan) << "\n";
    // std::cerr << "Target:\n" << ref_seq << "\n";
    // std::cerr << "Extend alignment result:\n" << aln->Verbose() << std::endl;

    if (aln->position().qend >= 0 && aln->position().tend >= 0) {
        new_anchor->qend(anchor->qend() + (aln->position().qend));
        new_anchor->tend(anchor->tend() + (aln->position().tend));
    }

    return new_anchor;
}

// std::shared_ptr<raptor::LocalPath> GraphMapper::NodesToPath_(
//                                             const std::shared_ptr<raptor::AnchorGraph> local_graph,
//                                             const std::vector<int64_t>& path_nodes,
//                                             int64_t path_score) {
//     if (path_nodes.size() == 0) {
//         return nullptr;
//     }
//     auto first_node = local_graph->GetNode(path_nodes[0]);
//     auto path = raptor::LocalPath::createPath(path_nodes[0], first_node);
//     for (size_t node_id = 1; node_id < path_nodes.size(); node_id++) {
//         // Although the generic Graph allows multi-edges, the AnchorGraph should not
//         // contain any. We will use only the first found edge between two nodes.
//         auto local_edges = local_graph->GetEdges(path_nodes[node_id - 1], path_nodes[node_id]);
//         auto node_data = local_graph->GetNode(path_nodes[node_id]);
//         auto edge_data = local_edges[0]->data();
//         path->Add(path_nodes[node_id], node_data, edge_data);
//     }
//     path->score(path_score);
//     return path;
// }

std::shared_ptr<raptor::LocalPath> PathAligner::FlankExtend(const mindex::IndexPtr& index,
                                                            const mindex::SequencePtr& qseq,
                                                            std::shared_ptr<raptor::AlignerBase>& aligner_ext,
                                                            int32_t flank_ext_len,
                                                            const std::shared_ptr<raptor::LocalPath>& path) {
    /*
     * This function shares similarities to PathAligner::Align in the section which
     * extends the front and back. It might even serve as a replacement for it in the future.
    */

    std::shared_ptr<raptor::LocalPath> ext_path = nullptr;

    if (path->nodes().size() == 0) {
        return ext_path;
    }

    //     auto& first_node = path->nodes().front();
    //     ext_path = raptor::LocalPath::createPath(first_node->name(), first_node->data());
    //     int64_t num_nodes = (int64_t) path->nodes().size();
    //     for (int64_t node_id = 1; node_id < (num_nodes); node_id++) {
    //         auto& node = path->nodes()[node_id];
    //         auto& edge_data = path->edges()[node_id-1]->data();
    //         ext_path->Add(node->name(), node->data(), edge_data);
    //     }
    //     ext_path->score(path->score());
    //     return ext_path;

    // std::cerr << "Extending path:\n";
    // std::cerr << path->Verbose() << "\n";

    // std::cerr << "qseq->header() = " << qseq->header();
    // There definitely is at least one node (due to the check above).
    if (path->nodes().size() == 1) {
        auto& node = path->nodes().back();
            // std::cerr << "Front (1):\n";
        auto new_anchor_1 = ExtendAlignmentFront_(index, qseq, aligner_ext, flank_ext_len, node->data());
            // std::cerr << "Back (1):\n";
        auto new_anchor_2 = ExtendAlignmentBack_(index, qseq, aligner_ext, flank_ext_len, new_anchor_1);

        ext_path = raptor::LocalPath::createPath(node->name(), new_anchor_2);

    } else {
        auto& first_node = path->nodes().front();
            // std::cerr << "Front (2):\n";
        auto new_first_anchor = ExtendAlignmentFront_(index, qseq, aligner_ext, flank_ext_len, first_node->data());
        ext_path = raptor::LocalPath::createPath(first_node->name(), new_first_anchor);

        int64_t num_nodes = (int64_t) path->nodes().size();
        for (int64_t node_id = 1; node_id < (num_nodes - 1); node_id++) {
            auto& node = path->nodes()[node_id];
            auto& edge_data = path->edges()[node_id-1]->data();
            ext_path->Add(node->name(), node->data(), edge_data);
        }

        auto& last_node = path->nodes().back();
            // std::cerr << "Back (2):\n";
        auto new_last_anchor = ExtendAlignmentBack_(index, qseq, aligner_ext, flank_ext_len, last_node->data());
        auto& last_edge_data = path->edges().back()->data();
        ext_path->Add(last_node->name(), new_last_anchor, last_edge_data);
    }

    ext_path->score(path->score());

    return ext_path;
}

std::shared_ptr<raptor::PathAlignment> PathAligner::Align(const mindex::SequencePtr& qseq,
                                                      const std::shared_ptr<raptor::LocalPath> path,
                                                      int32_t path_id, int32_t num_paths,
                                                      bool use_extend_alignment,
                                                      const std::shared_ptr<raptor::ParamsAligner> params) {
    TicToc total_time;
    total_time.start();

    std::shared_ptr<raptor::PathAlignment> result = raptor::createPathAlignment(path);

    if (path->nodes().size() == 0) {
        return result;
    }

    std::shared_ptr<raptor::LocalPath> ext_path = nullptr;

    if (use_extend_alignment) {
        // There definitely is at least one node (due to the check above).
        if (path->nodes().size() == 1) {
            auto& node = path->nodes().back();
            auto new_anchor_1 = ExtendAlignmentFront_(index_, qseq, aligner_ext_, -1, node->data());
            auto new_anchor_2 = ExtendAlignmentBack_(index_, qseq, aligner_ext_, -1, new_anchor_1);

            ext_path = raptor::LocalPath::createPath(node->name(), new_anchor_2);

        } else {
            auto& first_node = path->nodes().front();
            auto new_first_anchor = ExtendAlignmentFront_(index_, qseq, aligner_ext_, -1, first_node->data());
            ext_path = raptor::LocalPath::createPath(first_node->name(), new_first_anchor);

            int64_t num_nodes = (int64_t) path->nodes().size();
            for (int64_t node_id = 1; node_id < (num_nodes - 1); node_id++) {
                auto& node = path->nodes()[node_id];
                auto& edge_data = path->edges()[node_id-1]->data();
                ext_path->Add(node->name(), node->data(), edge_data);
            }

            auto& last_node = path->nodes().back();
            auto new_last_anchor = ExtendAlignmentBack_(index_, qseq, aligner_ext_, -1, last_node->data());
            auto& last_edge_data = path->edges().back()->data();
            ext_path->Add(last_node->name(), new_last_anchor, last_edge_data);
        }
    } else {
        ext_path = path;
    }




    const char* query = (const char*) &(qseq->data()[0]);
    std::string ref_seq;
    std::vector<AlignmentRegion> ref_seq_regions;
    bool rv_compose = ComposeTargetSequence_(qseq, ext_path, ref_seq, ref_seq_regions);

    if (rv_compose == false) {
        LOG_ALL("Problem composing the target sequence when processing qid: %d (%s). Skipping.\n", qseq->abs_id(), qseq->header().c_str());
        return result;
    }

    // Sanity check.
    if (ref_seq.size() == 0) {
        return result;
    }

    // Perform the actual alignment.
    int32_t qstart = ext_path->nodes().front()->data()->QueryStart();
    int32_t qend = ext_path->nodes().back()->data()->QueryEnd();

    #ifdef RAPTOR_TESTING_MODE
        DEBUG_QSEQ(params, qseq,
            LOG_ALL("Running aligner_->Global on qname = '%s', id = %d\n", qseq->header().c_str(), qseq->abs_id())
        );
        // fflush(stderr);
        // fprintf (stderr, "Query:\n'%s'\n", query);
        // fprintf (stderr, "Target:\n'%s'\n", ref_seq.c_str());
        DEBUG_QSEQ(params, qseq,
            LOG_ALL("qstart = %ld, qend = %ld, qseq->data().size() = %ld, ref_seq.size() = %ld\n", qstart, qend, qseq->data().size(), ref_seq.size())
        );
    #endif

    auto entire_alignment =
        aligner_->Global(query + qstart, qend - qstart, ref_seq.c_str(), ref_seq.size());

    #ifdef RAPTOR_TESTING_MODE
        if (params->debug_qid == qseq->abs_id() ||
            params->debug_qname == std::string(qseq->header())) {
            if (entire_alignment->status() == raptor::AlignmentReturnValue::Suboptimal) {
                LOG_ALL("Suboptimal");
            }

            LOG_ALL("Identity: %.2f\n", entire_alignment->op_counts().identity_min);
            LOG_ALL("Error rate: %.2f\n", entire_alignment->op_counts().error_rate * 100.0f);
            LOG_ALL("Error rate unique: %.2f\n", entire_alignment->op_counts().error_rate_u * 100.0f);
            LOG_ALL("Running SplitAlignment_.\n");
        }
    #endif

    if (entire_alignment->status() != raptor::AlignmentReturnValue::OK) {
        return result;
    }


    auto region_aln = SplitAlignment_(qseq, ext_path, ref_seq, ref_seq_regions, entire_alignment, params->min_identity);

    #ifdef RAPTOR_TESTING_MODE
        DEBUG_QSEQ(params, qseq, LOG_ALL("region_aln.size() = %ld\n", region_aln.size()));
    #endif


    if (region_aln.size() == 0) {
        return result;
    }

    // Sum the individual alignment scores.
    int64_t total_score = 0;
    // fprintf (stderr, "Entire identity: %.2f\n", entire_alignment->op_counts().identity_min * 100.0f);
    for (auto& aligned_region : region_aln) {
        total_score += aligned_region->aln()->score();

        #ifdef RAPTOR_TESTING_MODE
            if (params->debug_qid == qseq->abs_id() ||
                params->debug_qname == std::string(qseq->header())) {
                LOG_ALL("Identity: %.2f\n", aligned_region->aln()->op_counts().identity_min * 100.0f);
                LOG_ALL("Error rate: %.2f\n", aligned_region->aln()->op_counts().error_rate * 100.0f);
                LOG_ALL("Error rate unique: %.2f\n", aligned_region->aln()->op_counts().error_rate_u * 100.0f);
                LOG_ALL("\n");
            }
        #endif

    }

    // // Extend align the first and the last region.
    // // These should be in order, thanks to SplitAlignment_.
    // if (region_aln.size() > 0) {
    //     ExtendAlignmentFront_(region_aln.front());
    //     ExtendAlignmentBack_(region_aln.back());
    // }

    total_time.stop();

    // fprintf (stderr, "aln_total_time = %lf\n", total_time.get_microsecs());

    result->add_timing("aln_total_time", total_time.get_microsecs());
    result->entire_alignment(entire_alignment);
    result->path_score(total_score);
    result->alns(region_aln);

    // printf ("Aligned:\n");
    // printf ("  CIGAR: %s\n", CigarToString(entire_alignment->cigar).c_str());
    // printf ("  Edit dist: %ld\n", entire_alignment->edit_dist);
    // printf ("  Query len: %ld\n", qseq.get_sequence_length());
    // printf ("q:\n%s\n", qseq.GetSequenceAsString(entire_alignment->position.qstart,
    // entire_alignment->position.qend).c_str()); printf ("r:\n%s\n", ref_seq.c_str());
    // // printf ("chain->env()->t_id = %d, qstart = %d, qend = %d, rstart = %d, rend = %d,
    // chain->env()->t_rev = %d\n", chain->env()->t_id, qstart, qend, rstart, rend,
    // chain->env()->t_rev); fflush(stdout);

    return result;
}

}  // namespace raptor
