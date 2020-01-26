/*
 * graph_mapper_tools.cc
 *
 *  Created on: Jan 26, 2019
 *      Author: Ivan Sovic
 */

#include <raptor/mapper_tools.h>

#include <algorithm>
#include <cmath>
#include <utility/range_tools.hpp>
#include <lib/edlib.h>
#include <utility/revcmp.hpp>

namespace raptor {
namespace mapper {

void CalcHitCoverage(const std::vector<mindex::SeedHitPacked>& hits, int32_t seed_len,
                                  int32_t hits_begin, int32_t hits_end, int32_t& cov_bases_q,
                                  int32_t& cov_bases_r) {
    /*
      Expects the seed hits to be sorted!
    */
    cov_bases_q = cov_bases_r = 0;

    if (hits.size() == 0 || hits_begin >= hits.size() || (hits_end - hits_begin) <= 0) {
        return;
    }

    // Add the left part of the tile to the covered bases.
    hits_end = std::min(hits_end, (int32_t) hits.size());
    for (size_t i = (hits_begin + 1); i < hits_end; i++) {
        cov_bases_q += std::min(seed_len, (hits[i].QueryPos() - hits[i - 1].QueryPos()));
        cov_bases_r += std::min(seed_len, (hits[i].TargetPos() - hits[i - 1].TargetPos()));
    }

    // The last seed needs to be covered fully.
    cov_bases_q += seed_len;
    cov_bases_r += seed_len;
}

/*
 * The MakeAnchors_ method takes a vector of target hits, where each vector element contains
 * colinear hits for a single target/query pair. Each vector elements contains data for different
 * colinear hits (e.g. on different diagonal). Several vector elements can have hits to the
 * same target.
 * This method just takes the start and end query/target coordinates for each colinear hit range,
 * and groups them by target ID (keeping in mind the strand; fwd and rev strands are maintained
 * separately). One target ID will have a TargetHits<raptor::RegionMapped> object with (most likey) multiple
 * anchors in the `.hits` vector. These anchors lie on different diagonals.
 * A single TargetHits<raptor::RegionMapped> ".hits" entry represents one colinear range of seed hits.
 */
std::vector<std::shared_ptr<raptor::TargetAnchorType>> MakeAnchors(
    const std::vector<raptor::ChainPtr>& target_hits) {
    // Each target ID is one vector element. Each raptor::TargetHits has a .hits
    // vector with anchors.
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> target_anchors;

    // Just a lookup to identify where to place an anchor.
    std::unordered_map<int32_t, int32_t> t_id_map;

    if (target_hits.size() == 0) {
        return target_anchors;
    }

    int32_t num_created_anchors = 0;

    for (size_t target_hits_id = 0; target_hits_id < target_hits.size(); target_hits_id++) {
        assert(target_hits[target_hits_id]->hits().size() > 0 && "This shouldn't be zero, filtering should have removed zero-length target_hits.");

        auto& th = target_hits[target_hits_id];

        // The t_id of a chain could have clashes with rev cmp. Encoding
        // it again with the rev cmp flag will avoid this.
        int32_t key = th->env()->t_id << 1;
        key |= (th->env()->t_rev) ? 1 : 0;

        auto it = t_id_map.find(key);
        if (it == t_id_map.end()) {
            t_id_map[key] = target_anchors.size();

            auto anchor = std::shared_ptr<raptor::TargetAnchorType>(
                new raptor::TargetAnchorType(th->env()));

            target_anchors.emplace_back(anchor);
            it = t_id_map.find(key);
        }

        // Just fetch the info for an anchor.
        int32_t anchor_id = num_created_anchors;  // Storing the anchor ID will allow tracing the
                                                  // exact seed hits that make an alignment.
        ind_t qstart =
            th->hits().front().QueryPos();  // TODO: Removed this on 13.12.2017. when I
                                                        // added homopolymer suppression and
                                                        // accurate handing of hits! " -
                                                        // (index_->params()->k - 1);"  // Seed pos is
                                                        // the last base of the seed.
        ind_t qend = th->hits().back().QueryPos() + 1;  // Again, the last base of the seed. Add 1 to make the end position non-inclusive.
        ind_t rstart = th->hits().front().TargetPos();  // TODO: Removed this on 13.12.2017. when I
                                                         // added homopolymer suppression and
                                                         // accurate handing of hits! "-
                                                         // (index_->params()->k - 1);"  // Seed pos is
                                                         // the last base of the seed.
        ind_t rend = th->hits().back().TargetPos() + 1;  // Again, the last base of the seed. Add 1 to make the end position non-inclusive.
        int32_t cov_bases_q = th->cov_bases_q();
        int32_t cov_bases_t = th->cov_bases_t();
        int32_t num_seeds = th->hits().size();

        target_anchors[it->second]->hits().emplace_back(
            raptor::createRegionMapped(anchor_id, target_hits_id, th->env(), qstart, qend, rstart, rend,
                             cov_bases_q, cov_bases_t, num_seeds, -1, -1, -1, -1, -1, -1, 0, false));

        cov_bases_q += target_anchors[it->second]->cov_bases_q();
        cov_bases_t += target_anchors[it->second]->cov_bases_t();
        num_seeds += target_anchors[it->second]->num_seeds();

        target_anchors[it->second]->cov_bases_q(cov_bases_q);
        target_anchors[it->second]->cov_bases_t(cov_bases_t);
        target_anchors[it->second]->num_seeds(num_seeds);

        num_created_anchors += 1;
    }

    return target_anchors;
}

std::vector<raptor::ChainPtr> GroupTargetSeedHits(
            std::vector<mindex::SeedHitPacked> seed_hits,  // Copy.
            int32_t k,
            int32_t qid,
            int32_t qlen) {

    std::vector<raptor::ChainPtr> all_target_hits;

    // There is an "operator<" defined in the SeedHitPacked.
    std::sort(seed_hits.begin(), seed_hits.end());

    std::vector<std::pair<size_t, size_t>> ranges = istl::FindRanges<mindex::SeedHitPacked>(seed_hits,
                        [](const mindex::SeedHitPacked& a, const mindex::SeedHitPacked& b) {
                                return a.TargetId() == b.TargetId(); });

    for (const auto& range_pair: ranges) {
        size_t range_start = std::get<0>(range_pair);
        size_t range_end = std::get<1>(range_pair);

        if (range_end <= range_start) {
            continue;
        }

        int32_t t_id = seed_hits[range_start].TargetId();
        int32_t t_len = 0;
        bool t_rev = seed_hits[range_start].TargetRev();

        std::shared_ptr<raptor::MappingEnv> new_env = raptor::createMappingEnv(
                t_id, 0, t_len, t_rev, qid, qlen, false);

        auto single_target_hits = raptor::ChainPtr(
                            new raptor::TargetHits<mindex::SeedHitPacked>(new_env));

        for (size_t seed_id = range_start; seed_id < range_end; ++seed_id) {
            single_target_hits->hits().emplace_back(seed_hits[seed_id]);
        }
        single_target_hits->score(0);
        int32_t cov_bases_q = 0;
        int32_t cov_bases_t = 0;
        raptor::mapper::CalcHitCoverage(single_target_hits->hits(), k, 0, single_target_hits->hits().size(), cov_bases_q, cov_bases_t);
        single_target_hits->cov_bases_q(cov_bases_q);
        single_target_hits->cov_bases_t(cov_bases_t);

        all_target_hits.emplace_back(single_target_hits);
    }

    return all_target_hits;
}

/*
 * Aligns the anchor using Edlib to obtain match rates needed for the DP chaining.
 */
int32_t CalcMatchRate(const mindex::SequencePtr& qseq, mindex::IndexPtr index,
                                const std::shared_ptr<raptor::RegionMapped>& anchor) {
    int32_t diffs = 0, matches = 0;

    LOG_ALL("CalcMatchRate is currently deactivated while refactoring MinimizerIndex.\n");

    // int32_t tseq_start =
    //     anchor->TargetIndexStart() +
    //     ((!anchor->TargetRev()) ? (anchor->TargetStart()) : (anchor->TargetLen() - anchor->TargetEnd()));

    // size_t tseq_len = anchor->TargetEnd() - anchor->TargetStart();

    // size_t qseq_len = anchor->QueryEnd() - anchor->QueryStart();

    // int32_t result = 0;

    // EdlibAlignTask task = EDLIB_TASK_DISTANCE;
    // EdlibAlignMode edlib_mode = EDLIB_MODE_NW;

    // // If not reverse complement, we can speed-up the computation by not copying.
    // if (!anchor->TargetRev()) {
    //     EdlibAlignResult result = edlibAlign((const char*)&(qseq->data()[anchor->QueryStart()]),
    //                                          qseq_len, (const char*)&index->data()[tseq_start],
    //                                          tseq_len, edlibNewAlignConfig(-1, edlib_mode, task));
    //     diffs = result.editDistance;
    //     matches = qseq_len - diffs;  // An approximation. Pessimistic.
    //     edlibFreeAlignResult(result);

    // } else {
    //     std::string target = raptor::ReverseComplement((const char*)&index->data()[tseq_start], tseq_len);

    //     EdlibAlignResult result =
    //         edlibAlign((const char*)&(qseq->data()[anchor->QueryStart()]), qseq_len, target.c_str(),
    //                    target.size(), edlibNewAlignConfig(-1, edlib_mode, task));
    //     diffs = result.editDistance;
    //     matches = qseq_len - diffs;  // An approximation. Pessimistic.
    //     edlibFreeAlignResult(result);
    // }

    // // In case something went wrong with the alignment, the anchor
    // // just won't be used.
    // if (result < 0) {
    //     matches = 0;
    // }

    return matches;
}

}
}
