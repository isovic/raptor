#include <gtest/gtest.h>
#include <tests/unit/test_graph_mapper.h>

#include <raptor/graph_mapper.h>
#include <graph/segment_graph_parser.h>

#include <cstdint>
#include <algorithm>
#include <set>
#include <vector>
#include <log/log_tools.h>
#include <containers/mapping_result/linear_mapping_result.h>
#include <utility/range_tools.hpp>
#include <raptor/mapper_tools.h>
#include <containers/mapping_env.h>
#include <graph/split_segment_graph.h>
#include <raptor/graph_mapper_tools.h>

namespace raptor {
namespace unit {

bool VERBOSE_DEBUG_QID_GM_TOOLS = false;

std::vector<raptor::ChainPtr> WrapGroupTargetSeedHits(
            std::vector<std::pair<int32_t, mindex::SeedHitPacked>> seed_hits,  // Copy.
            int32_t k,
            int32_t qid,
            int32_t qlen) {

    std::vector<raptor::ChainPtr> all_target_hits;

    // There is an "operator<" defined in the SeedHitPacked.
    std::sort(seed_hits.begin(), seed_hits.end());

    std::vector<std::pair<size_t, size_t>> ranges = istl::FindRanges<std::pair<int32_t, mindex::SeedHitPacked>>(seed_hits,
                        [](const std::pair<int32_t, mindex::SeedHitPacked>& a, const std::pair<int32_t, mindex::SeedHitPacked>& b) {
                                return std::get<0>(a) == std::get<0>(b); });

    for (const auto& range_pair: ranges) {
        size_t range_start = std::get<0>(range_pair);
        size_t range_end = std::get<1>(range_pair);

        if (range_end <= range_start) {
            continue;
        }

        const auto& seed_hit_first = std::get<1>(seed_hits[range_start]);

        int32_t t_id = seed_hit_first.TargetId();
        bool t_rev = seed_hit_first.TargetRev();
        int32_t t_len = 0;

        std::shared_ptr<raptor::MappingEnv> new_env = raptor::createMappingEnv(
                t_id, 0, t_len, t_rev, qid, qlen, false);

        auto single_target_hits = raptor::ChainPtr(
                            new raptor::TargetHits<mindex::SeedHitPacked>(new_env));

        for (size_t seed_id = range_start; seed_id < range_end; ++seed_id) {
            const auto& seed_hit = std::get<1>(seed_hits[seed_id]);
            // Sanity check that the test is well defined. The the t_id should be the same
            // for every member of the group.
            EXPECT_EQ(seed_hit.TargetId(), t_id);
            single_target_hits->hits().emplace_back(seed_hit);
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
 * Aside from the graphs (SegmentGraph and SplitSegmentGraph), the
 * LinearMappingResult is the only input to GraphMapper.
 * It consists of a vector of anchors, and a vector of seed target hits, plus
 * the environment specification.
 * The only reason that the seed_hits are required for GraphMapper is so that
 * the anchors can be broken on graph in/out edges.
 * The reason why the GraphMapper accepts the SplitSegmentGraph instead of constructing
 * it is the speed, since the same SSG is used by all threads for all reads.
 * The reason why SegmentGraph is required is legacy code. Perhaps it can be
 * refactored in the future to drop it and simplify the interface.
*/
std::shared_ptr<raptor::LinearMappingResult> UtilCreateLinearMappingResult(
                            mindex::IndexPtr index, int32_t k, const std::string& qname,
                            int32_t qlen, int32_t qid,
                            const std::vector<std::pair<int32_t, mindex::SeedHitPacked>>& seed_hits) {

    // Create an empty container.
    std::shared_ptr<raptor::LinearMappingResult> result = raptor::createMappingResult(qid, qlen, qname, index);

    // Take the plain flat list of seeds and group them by target.
    std::vector<raptor::ChainPtr> hits = raptor::unit::WrapGroupTargetSeedHits(seed_hits, k, qid, qlen);

    // Anchors will be formed automatically from hits.
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> anchors = raptor::mapper::MakeAnchors(hits);

    // Set the linear mapping result.
    result->target_anchors(anchors);
    result->target_hits(hits);
    result->return_value(raptor::MapperReturnValueBase::OK);

    return result;
}

TEST(GraphMapperTools, BreakAnchors1) {
    /*
     * Empty input.
    */
    std::string test_name("BreakAnchors1");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_GM_TOOLS);

    // raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    // Ther first element in pair is the target group, and the second is the MinimizerHit.
    std::vector<std::pair<int32_t, mindex::SeedHitPacked>> seed_hits {
    };

    mindex::IndexPtr index = nullptr;
    int32_t k = 15;
    std::string qname("query1");
    int32_t qlen = 1000;
    int32_t qid = 0;

    std::shared_ptr<raptor::LinearMappingResult> linear_mapping = UtilCreateLinearMappingResult(
                                index, k, qname,
                                qlen, qid, seed_hits);

    std::shared_ptr<raptor::LinearMappingResult> result_linear_mapping = raptor::graphmapper::BreakAnchors(
                        seg_graph,
                        linear_mapping, k);

    std::string result_str = result_linear_mapping->WriteAsCSV(',');

    std::string expected_str = "";
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(GraphMapperTools, BreakAnchors2) {
    /*
     * Test with an empty graph and linear mapping. No anchors should be broken
     * and everything should be as is.
    */
    std::string test_name("BreakAnchors2");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_GM_TOOLS);

    // raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    // Ther first element in pair is the target group, and the second is the MinimizerHit.
    std::vector<std::pair<int32_t, mindex::SeedHitPacked>> seed_hits {
            {0, {0,0,2,0,2}},
            {0, {0,0,12,0,12}},
            {0, {0,0,25,0,25}},
            {0, {0,0,40,0,40}},
            {0, {0,0,55,0,55}},
            {0, {0,0,67,0,67}},
            {0, {0,0,80,0,80}},
            {0, {0,0,95,0,95}},
            {0, {0,0,109,0,109}},
            {0, {0,0,125,0,125}},
            {0, {0,0,135,0,135}},
            {0, {0,0,150,0,150}},
            {0, {0,0,163,0,163}},
            {0, {0,0,178,0,178}},
            {0, {0,0,192,0,192}},
            {0, {0,0,212,0,212}},
            {0, {0,0,220,0,220}},
            {0, {0,0,235,0,235}},
            {0, {0,0,247,0,247}},
            {0, {0,0,262,0,262}},
            {0, {0,0,274,0,274}},
            {0, {0,0,292,0,292}},
            {0, {0,0,302,0,302}},
            {0, {0,0,317,0,317}},
            {0, {0,0,329,0,329}},
            {0, {0,0,343,0,343}},
            {0, {0,0,355,0,355}},
            {0, {0,0,367,0,367}},
            {0, {0,0,382,0,382}},
            {0, {0,0,397,0,397}},
            {0, {0,0,411,0,411}},
            {0, {0,0,425,0,425}},
            {0, {0,0,439,0,439}},
            {0, {0,0,450,0,450}},
            {0, {0,0,467,0,472}},
            {0, {0,0,480,0,485}},
            {0, {0,0,495,0,500}},
            {0, {0,0,509,0,514}},
            {0, {0,0,524,0,529}},
            {0, {0,0,536,0,541}},
            {0, {0,0,551,0,556}},
            {0, {0,0,566,0,571}},
            {0, {0,0,580,0,585}},
            {0, {0,0,605,0,600}},
            {0, {0,0,609,0,604}},
            {0, {0,0,623,0,618}},
            {0, {0,0,638,0,633}},
            {0, {0,0,652,0,647}},
            {0, {0,0,664,0,659}},
            {0, {0,0,679,0,674}},
    };

    mindex::IndexPtr index = nullptr;
    int32_t k = 15;
    std::string qname("query1");
    int32_t qlen = 1000;
    int32_t qid = 0;

    std::shared_ptr<raptor::LinearMappingResult> linear_mapping = UtilCreateLinearMappingResult(
                                index, k, qname,
                                qlen, qid, seed_hits);

    std::shared_ptr<raptor::LinearMappingResult> result_linear_mapping = raptor::graphmapper::BreakAnchors(
                        seg_graph,
                        linear_mapping, k);

    std::string result_str = result_linear_mapping->WriteAsCSV(',');

    std::string expected_str =
            // H (for hit), group_id, tid, trev, tpos, qmask, qpos
            "H,0,0,0,2,0,2\n"
            "H,0,0,0,12,0,12\n"
            "H,0,0,0,25,0,25\n"
            "H,0,0,0,40,0,40\n"
            "H,0,0,0,55,0,55\n"
            "H,0,0,0,67,0,67\n"
            "H,0,0,0,80,0,80\n"
            "H,0,0,0,95,0,95\n"
            "H,0,0,0,109,0,109\n"
            "H,0,0,0,125,0,125\n"
            "H,0,0,0,135,0,135\n"
            "H,0,0,0,150,0,150\n"
            "H,0,0,0,163,0,163\n"
            "H,0,0,0,178,0,178\n"
            "H,0,0,0,192,0,192\n"
            "H,0,0,0,212,0,212\n"
            "H,0,0,0,220,0,220\n"
            "H,0,0,0,235,0,235\n"
            "H,0,0,0,247,0,247\n"
            "H,0,0,0,262,0,262\n"
            "H,0,0,0,274,0,274\n"
            "H,0,0,0,292,0,292\n"
            "H,0,0,0,302,0,302\n"
            "H,0,0,0,317,0,317\n"
            "H,0,0,0,329,0,329\n"
            "H,0,0,0,343,0,343\n"
            "H,0,0,0,355,0,355\n"
            "H,0,0,0,367,0,367\n"
            "H,0,0,0,382,0,382\n"
            "H,0,0,0,397,0,397\n"
            "H,0,0,0,411,0,411\n"
            "H,0,0,0,425,0,425\n"
            "H,0,0,0,439,0,439\n"
            "H,0,0,0,450,0,450\n"
            "H,0,0,0,467,0,472\n"
            "H,0,0,0,480,0,485\n"
            "H,0,0,0,495,0,500\n"
            "H,0,0,0,509,0,514\n"
            "H,0,0,0,524,0,529\n"
            "H,0,0,0,536,0,541\n"
            "H,0,0,0,551,0,556\n"
            "H,0,0,0,566,0,571\n"
            "H,0,0,0,580,0,585\n"
            "H,0,0,0,605,0,600\n"
            "H,0,0,0,609,0,604\n"
            "H,0,0,0,623,0,618\n"
            "H,0,0,0,638,0,633\n"
            "H,0,0,0,652,0,647\n"
            "H,0,0,0,664,0,659\n"
            "H,0,0,0,679,0,674\n"
            // A (for anchor), group_id, qid, tid, score, num_seeds, qrev, qstart, qend, qlen, trev, tstart, tend, tlen, cov_bases_q, cov_bases_t, num_seeds
            "A,0,0,0,671,50,0,2,675,1000,0,2,680,0,671,671,50\n"
        ;
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(GraphMapperTools, BreakAnchors3) {
    /*
     * Test on a more complicated graph with indels, SNPs and multiple
     * edges stemming out of the same nodes.
    */
    std::string test_name("BreakAnchors3");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "695", "0"},
        {"snp1", "1", "1"},
        {"snp2", "1", "2"},
        {"snp3", "1", "3"},
        {"sv1", "5", "4"},
        {"sv2", "10", "5"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "E1-in", "ref1+", "snp1+", "124", "124", "0", "0", "*"},
        {"E", "E1-out", "snp1+", "ref1+", "1$", "1$", "125", "125", "*"},
        {"E", "E2-in", "ref1+", "snp2+", "211", "211", "0", "0", "*"},
        {"E", "E2-out", "snp2+", "ref1+", "1$", "1$", "212", "212", "*"},
        {"E", "E3-in", "ref1+", "snp3+", "291", "291", "0", "0", "*"},
        {"E", "E3-out", "snp3+", "ref1+", "1$", "1$", "292", "292", "*"},
        {"E", "E4-sv1-in", "ref1+", "sv1+", "464", "464", "0", "0", "*"},
        {"E", "E4-sv1-out", "sv1+", "ref1+", "5$", "5$", "464", "464", "*"},
        {"E", "E5-sv2", "ref1+", "ref1+", "596", "596", "605", "605", "*"},
        {"E", "E6", "ref1+", "ref1+", "464", "464", "596", "596", "*"}, // An edge to induce multiple forking at coords 464 and 596.
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_GM_TOOLS);

    // raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    // std::cerr << "Testing JSON output, seg_graph:\n" << seg_graph->ToJSON() << "\n";
    // std::cerr << "Testing JSON output, ssg:\n" << ssg->ToJSON() << "\n";

    // Ther first element in pair is the target group, and the second is the MinimizerHit.
    std::vector<std::pair<int32_t, mindex::SeedHitPacked>> seed_hits {
            {0, {0,0,2,0,2}},
            {0, {0,0,12,0,12}},
            {0, {0,0,25,0,25}},
            {0, {0,0,40,0,40}},
            {0, {0,0,55,0,55}},
            {0, {0,0,67,0,67}},
            {0, {0,0,80,0,80}},
            {0, {0,0,95,0,95}},
            {0, {0,0,109,0,109}},
            {0, {0,0,125,0,125}},
            {0, {0,0,135,0,135}},
            {0, {0,0,150,0,150}},
            {0, {0,0,163,0,163}},
            {0, {0,0,178,0,178}},
            {0, {0,0,192,0,192}},
            {0, {0,0,212,0,212}},
            {0, {0,0,220,0,220}},
            {0, {0,0,235,0,235}},
            {0, {0,0,247,0,247}},
            {0, {0,0,262,0,262}},
            {0, {0,0,274,0,274}},
            {0, {0,0,292,0,292}},
            {0, {0,0,302,0,302}},
            {0, {0,0,317,0,317}},
            {0, {0,0,329,0,329}},
            {0, {0,0,343,0,343}},
            {0, {0,0,355,0,355}},
            {0, {0,0,367,0,367}},
            {0, {0,0,382,0,382}},
            {0, {0,0,397,0,397}},
            {0, {0,0,411,0,411}},
            {0, {0,0,425,0,425}},
            {0, {0,0,439,0,439}},
            {0, {0,0,450,0,450}},
            {0, {0,0,467,0,472}},
            {0, {0,0,480,0,485}},
            {0, {0,0,495,0,500}},
            {0, {0,0,509,0,514}},
            {0, {0,0,524,0,529}},
            {0, {0,0,536,0,541}},
            {0, {0,0,551,0,556}},
            {0, {0,0,566,0,571}},
            {0, {0,0,580,0,585}},
            {0, {0,0,605,0,600}},
            {0, {0,0,609,0,604}},
            {0, {0,0,623,0,618}},
            {0, {0,0,638,0,633}},
            {0, {0,0,652,0,647}},
            {0, {0,0,664,0,659}},
            {0, {0,0,679,0,674}},
    };

    mindex::IndexPtr index = nullptr;
    int32_t k = 15;
    std::string qname("query1");
    int32_t qlen = 1000;
    int32_t qid = 0;

    std::shared_ptr<raptor::LinearMappingResult> linear_mapping = UtilCreateLinearMappingResult(
                                index, k, qname,
                                qlen, qid, seed_hits);

    std::shared_ptr<raptor::LinearMappingResult> result_linear_mapping = raptor::graphmapper::BreakAnchors(
                        seg_graph,
                        linear_mapping, k);

    // std::cerr << result_linear_mapping->WriteAsCSV(',') << "\n";

    std::string result_str = result_linear_mapping->WriteAsCSV(',');

    std::string expected_str =
            // H (for hit), group_id, tid, trev, tpos, qmask, qpos
            "H,0,0,0,2,0,2\n"
            "H,0,0,0,12,0,12\n"
            "H,0,0,0,25,0,25\n"
            "H,0,0,0,40,0,40\n"
            "H,0,0,0,55,0,55\n"
            "H,0,0,0,67,0,67\n"
            "H,0,0,0,80,0,80\n"
            "H,0,0,0,95,0,95\n"
            "H,0,0,0,109,0,109\n"
            "H,1,0,0,125,0,125\n"
            "H,1,0,0,135,0,135\n"
            "H,1,0,0,150,0,150\n"
            "H,1,0,0,163,0,163\n"
            "H,1,0,0,178,0,178\n"
            "H,1,0,0,192,0,192\n"
            "H,2,0,0,212,0,212\n"
            "H,2,0,0,220,0,220\n"
            "H,2,0,0,235,0,235\n"
            "H,2,0,0,247,0,247\n"
            "H,2,0,0,262,0,262\n"
            "H,2,0,0,274,0,274\n"
            "H,3,0,0,292,0,292\n"
            "H,3,0,0,302,0,302\n"
            "H,3,0,0,317,0,317\n"
            "H,3,0,0,329,0,329\n"
            "H,3,0,0,343,0,343\n"
            "H,3,0,0,355,0,355\n"
            "H,3,0,0,367,0,367\n"
            "H,3,0,0,382,0,382\n"
            "H,3,0,0,397,0,397\n"
            "H,3,0,0,411,0,411\n"
            "H,3,0,0,425,0,425\n"
            "H,3,0,0,439,0,439\n"
            "H,3,0,0,450,0,450\n"
            "H,4,0,0,467,0,472\n"
            "H,4,0,0,480,0,485\n"
            "H,4,0,0,495,0,500\n"
            "H,4,0,0,509,0,514\n"
            "H,4,0,0,524,0,529\n"
            "H,4,0,0,536,0,541\n"
            "H,4,0,0,551,0,556\n"
            "H,4,0,0,566,0,571\n"
            "H,4,0,0,580,0,585\n"
            "H,5,0,0,605,0,600\n"
            "H,5,0,0,609,0,604\n"
            "H,5,0,0,623,0,618\n"
            "H,5,0,0,638,0,633\n"
            "H,5,0,0,652,0,647\n"
            "H,5,0,0,664,0,659\n"
            "H,5,0,0,679,0,674\n"
            // A (for anchor), group_id, qid, tid, score, num_seeds, qrev, qstart, qend, qlen, trev, tstart, tend, tlen, cov_bases_q, cov_bases_t, num_seeds
            "A,0,0,0,122,9,0,2,110,1000,0,2,110,0,122,122,9\n"
            "A,0,0,0,82,6,0,125,193,1000,0,125,193,0,82,82,6\n"
            "A,0,0,0,77,6,0,212,275,1000,0,212,275,0,77,77,6\n"
            "A,0,0,0,173,13,0,292,451,1000,0,292,451,0,173,173,13\n"
            "A,0,0,0,128,9,0,472,586,1000,0,467,581,0,128,128,9\n"
            "A,0,0,0,89,7,0,600,675,1000,0,605,680,0,89,89,7\n"
        ;
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

}
}
