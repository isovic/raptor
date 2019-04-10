// /*
//  * anchor_aligner.cc
//  *
//  *  Created on: Aug 23, 2017
//  *      Author: isovic
//  */

// #include "aligner/anchor_aligner.h"
// #include "aligner/aligner_util.hpp"
// #include <algorithm>

// namespace raptor {

// std::shared_ptr<AnchorAligner> createAnchorAligner(std::shared_ptr<raptor::AlignerBase> aligner)
// {
//   return std::shared_ptr<AnchorAligner>(new AnchorAligner(aligner));
// }

// AnchorAligner::AnchorAligner(std::shared_ptr<raptor::AlignerBase> aligner) : aligner_(aligner) {

// }

// AnchorAligner::~AnchorAligner() {

// }

// std::shared_ptr<AlignmentResult> AnchorAligner::GlobalEndToEnd(const char *query, int64_t qlen,
// const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
//   auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);
//   std::vector<AlignmentAnchor> final_anchors;
//   if (anchors.size() > 0) {
//     final_anchors.emplace_back(AlignmentAnchor(anchors.front().qstart, anchors.back().qend, qlen,
//     anchors.front().rstart, anchors.back().rend, rlen, false));
//   }
//   return GlobalAnchored(query, qlen, ref, rlen, final_anchors);
// }

// std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchored(const char *query, int64_t qlen,
// const char *ref, int64_t rlen, const std::vector<AlignmentAnchor>& anchors) {
//   auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

//   if (anchors.size() == 0) {
//     return result;
//   }

//   result->cigar.clear();
//   result->score = 0;

//   // Align between anchors.
//   for (int64_t i = 0; i < (anchors.size() - 1); i++) {
//     auto aln_result = aligner_->Global(query + anchors[i].qstart, anchors[i+1].qend -
//     anchors[i].qstart,
//                       ref + anchors[i].rstart, anchors[i+1].rend - anchors[i].rstart);

//     auto left_part = ExtractCigarBetweenQueryCoords(aln_result->cigar,
//                                                     0,
//                                                     anchors[i+1].qstart - anchors[i].qstart); //
//                                                     Leave next anchor for the next alignment.

//     result->cigar.insert(result->cigar.end(), left_part.begin(), left_part.end());

//     // TODO: This is wrong, because it takes the score for the next anchor twice. Need to rework
//     this. result->score += aln_result->score;
//   }

//   // Align the last anchor.
//   auto aln_result = aligner_->Global(query + anchors.back().qstart, anchors.back().qend -
//   anchors.back().qstart,
//                     ref + anchors.back().rstart, anchors.back().rend - anchors.back().rstart);

//   result->cigar.insert(result->cigar.end(), aln_result->cigar.begin(), aln_result->cigar.end());

//   // Add the soft clippings at front and back.
//   if ( anchors.front().qstart > 0) {
//     result->cigar.insert(result->cigar.begin(), raptor::CigarOp('S', anchors.front().qstart));
//   }
//   if ((qlen - anchors.front().qend) > 0) {
//     result->cigar.insert(result->cigar.end(), raptor::CigarOp('S', (qlen -
//     anchors.back().qend)));
//   }

//   result->score += aln_result->score;

//   // Fill the other alignment info.
//   result->edit_dist = EditDistFromExtCIGAR(result->cigar);

//   result->position = raptor::AlignmentPosition(anchors.front().qstart, anchors.back().qend,
//   anchors.front().rstart, anchors.back().rend); result->k = -1; result->rv =
//   raptor::AlignmentReturnValue::OK;

//   return result;
// }

// std::shared_ptr<AlignmentResult> AnchorAligner::GlobalAnchoredWithExtend(const char *query,
// int64_t qlen,
//                                                                          const char *ref, int64_t
//                                                                          rlen, const
//                                                                          std::vector<AlignmentAnchor>&
//                                                                          anchors, int32_t
//                                                                          bandwidth, int32_t
//                                                                          zdrop) {
//   auto result = std::shared_ptr<AlignmentResult>(new AlignmentResult);

//   if (anchors.size() == 0) {
//     return result;
//   }

//   std::vector<AlignmentAnchor> updated_anchors = anchors;

//   // Extend align last anchor forward to find max position.
//   // Align the last anchor.
//   auto ext_back = aligner_->Extend(query + anchors.back().qstart, qlen - anchors.back().qstart,
//                    ref + anchors.back().rstart, std::min(rlen - anchors.back().rstart, (qlen -
//                    anchors.back().qstart) * 2), bandwidth, zdrop);
//   // Get the correct coordinates from the alignment.
//   int64_t max_q_pos_back = ext_back->max_q_pos + anchors.back().qstart + 1; // The "+1" because
//   end coordinate is non-inclusive in Raptor. int64_t max_t_pos_back = ext_back->max_t_pos +
//   anchors.back().rstart + 1; // The max position is inclusive on the other hand.
//   // If extend did not pan out (e.g. band is too narrow), do not extend.
//   if (ext_back->max_q_pos >= 0 && ext_back->max_t_pos >= 0) {
//     updated_anchors.back().qend = max_q_pos_back;
//     updated_anchors.back().rend = max_t_pos_back;
//   }

//   // Extension of the front part.
//   // Create a reverse copy (but not complemented) of the query and target.
//   std::string rev_q_front, rev_t_front;
//   std::reverse_copy(query , query + anchors.front().qend, std::back_inserter(rev_q_front));
//   std::reverse_copy(ref + std::max((int64_t) 0, anchors.front().rend - 2 * anchors.front().qend),
//                     ref + anchors.front().rend, std::back_inserter(rev_t_front));
//   auto ext_front = aligner_->Extend(rev_q_front.c_str(), rev_q_front.size(),
//                    rev_t_front.c_str(), rev_t_front.size(),
//                    bandwidth, zdrop);
//   // Get the correct coordinates from the alignment.
//   int64_t max_q_pos_front = anchors.front().qend - (ext_front->max_q_pos + 1); // The "+1"
//   because end coordinate is non-inclusive in Raptor. int64_t max_t_pos_front =
//   anchors.front().rend - (ext_front->max_t_pos + 1); // The max position is inclusive on the
//   other hand.
//   // If extend did not pan out (e.g. band is too narrow), do not extend.
//   if (ext_front->max_q_pos >= 0 && ext_front->max_t_pos >= 0) {
//     updated_anchors.front().qstart = max_q_pos_front;
//     updated_anchors.front().rstart = max_t_pos_front;
//   }

//   // printf ("####################################\n");
//   // // printf ("Front rev Q:\n%s\n\nFront rev T:\n%s\n\n", rev_q_front.c_str(),
//   rev_t_front.c_str());
//   // // printf ("Back fwd Q:\n");
//   // // for (int64_t i=anchors.back().qend; i<qlen; i++) {
//   // //   printf ("%c", query[i]);
//   // // }
//   // // printf ("\n\n");
//   // // printf ("Back fwd T:\n");
//   // // for (int64_t i=anchors.back().rend; i<rlen; i++) {
//   // //   printf ("%c", ref[i]);
//   // // }
//   // // printf ("\n\n");

//   // for (int64_t i=0; i<anchors.size(); i++) {
//   //   printf ("  anchor[%ld]: qstart = %ld, qend = %ld, qlen = %ld, rstart = %ld, rend = %ld,
//   rlen = %ld\n",
//   //           i, anchors[i].qstart, anchors[i].qend, qlen, anchors[i].rstart, anchors[i].rend);
//   // }
//   // printf ("\n");
//   // printf ("max_q_pos_back = %ld\n", max_q_pos_back);
//   // printf ("max_t_pos_back = %ld\n", max_t_pos_back);
//   // printf ("ext_back->max_q_pos = %ld\n", ext_back->max_q_pos);
//   // printf ("ext_back->max_t_pos = %ld\n", ext_back->max_t_pos);
//   // printf ("####################################\n");
//   // printf ("\n");
//   // printf ("max_q_pos_front = %ld\n", max_q_pos_front);
//   // printf ("max_t_pos_front = %ld\n", max_t_pos_front);
//   // printf ("ext_front->max_q_pos = %ld\n", ext_front->max_q_pos);
//   // printf ("ext_front->max_t_pos = %ld\n", ext_front->max_t_pos);
//   // printf ("####################################\n");
//   // fflush(stdout);

//   // Align the updated coordinates.
//   return GlobalAnchored(query, qlen, ref, rlen, updated_anchors);
// }

// }
