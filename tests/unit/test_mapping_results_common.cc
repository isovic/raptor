#include <gtest/gtest.h>

#include <containers/mapping_result/mapping_result_common.h>
#include <containers/region/region_mapped.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

std::vector<std::shared_ptr<raptor::RegionBase>> HelperRunRelabelSupplementary(
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions,
    int32_t qlen, int32_t tlen,
    double min_sec_to_prim_ratio, int32_t allowed_overlap_bp) {
    // Dummy mapping env, needed for constructing the regions.
    std::shared_ptr<raptor::MappingEnv> env =
            raptor::createMappingEnv(0, 0, tlen, false,    // t_id, index_t_start, t_len, t_rev
                             0, qlen, false);       // q_id, q_len, q_rev
    // Initialize the alns based on the regions and the mapping env.
    std::vector<std::shared_ptr<raptor::RegionBase>> alns;
    for (int32_t i = 0; i < static_cast<int32_t>(input_regions.size()); ++i) {
        const auto& reg = input_regions[i];
        std::shared_ptr<raptor::RegionBase> new_aln = raptor::createRegionMapped(
                i, -1, env,  // id, target_hits_id, env
                std::get<0>(reg), std::get<1>(reg), std::get<2>(reg), std::get<3>(reg), // qstart, qend, tstart, tend
                -1, -1, -1, -1, std::get<4>(reg),         // cov_bases_q, cov_bases_t, num_seeds, edit_dist, score
                -1, -1, -1, -1,                 // path_id, num_paths, segment_id, num_segments
                std::get<5>(reg), std::get<6>(reg));    // region_priority, region_is_supplementary
        alns.emplace_back(new_aln);
    }

    // Convert the mapped regions into interval trees.
    IntervalVectorInt64 qi_prim;
    std::unordered_map<int64_t, IntervalVectorInt64> ti_prim;
    IntervalTreeInt64 qt_prim;
    std::unordered_map<int64_t, IntervalTreeInt64> tt_prim;
    CreateRegionIntervalTrees(alns, [](int32_t a){ return a == 0; },
                                qi_prim, qt_prim, ti_prim, tt_prim);

    // Label the supplementary, and mark them as primary too.
    // FindSupplementary(qt_prim, qi_prim, tt_prim, ti_prim, alns, min_sec_to_prim_ratio, allowed_overlap_bp);
    RelabelSupplementary(alns, min_sec_to_prim_ratio, allowed_overlap_bp, allowed_overlap_bp);

    return alns;
}

TEST(MappingResultsCommon, RelabelSupplementaryEmptyInput) {
    /*
     * There are zero input regions. Nothing should happen.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryNormal) {
    /*
     * Normal test case where we have 1 primary alignment, 2 supplementary and 1
     * secondary which aligns with one of the supplementary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,47",
        "0,0,-1,-1,0,1000,2000,5000,0,1000,2000,4000,-1,-1,-1,0,1,47",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,1,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryOverlappingSupp) {
    /*
     * One region slightly overlaps the primary region. Since allowed_overlap_bp is 0,
     * this region will be marked as secondary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {990, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,3",
        "0,0,-1,-1,0,990,2000,5000,0,1000,2000,4000,-1,-1,-1,1,0,3",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,2,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryAllowedOverlapInQuery) {
    /*
     * One region slightly overlaps the primary region.
     * The allowed_overlap_bp is 50, which means that this region should be marked
     * as supplementary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 50;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {990, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,49",
        "0,0,-1,-1,0,990,2000,5000,0,1000,2000,4000,-1,-1,-1,0,1,49",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,1,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryAllowedOverlapInQueryButAboveLimit) {
    /*
     * One region slightly overlaps the primary region.
     * The allowed_overlap_bp is 7 and overlap is of 50bp, which means that this
     * region should be marked as secondary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 7;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {990, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,3",
        "0,0,-1,-1,0,990,2000,5000,0,1000,2000,4000,-1,-1,-1,1,0,3",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,2,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryAllowedOverlapInTarget) {
    /*
     * One region slightly overlaps the primary region.
     * The allowed_overlap_bp is 50, which means that this region should be marked
     * as supplementary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 50;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 990, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,48",
        "0,0,-1,-1,0,1000,2000,5000,0,990,2000,4000,-1,-1,-1,0,1,48",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,1,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryAllowedOverlapInTargetButAboveLimit) {
    /*
     * One region slightly overlaps the primary region.
     * The allowed_overlap_bp is 7 and overlap is of 50bp, which means that this
     * region should be marked as secondary.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 7;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 990, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,3",
        "0,0,-1,-1,0,1000,2000,5000,0,990,2000,4000,-1,-1,-1,1,0,3",
        "0,0,-1,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,3",
        "0,0,-1,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,2,0,3"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryNormalMapq2) {
    /*
     * Normal test case. There are multiple secondary regions which overlap one
     * supplementary region.
     * The mapping quality of these regions should be equal to 2.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 6000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
        {5000, 6000, 3000, 4000, score, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,6000,0,0,1000,4000,-1,-1,-1,0,0,32",
        "0,0,-1,-1,0,1000,2000,6000,0,1000,2000,4000,-1,-1,-1,0,1,32",
        "0,0,-1,-1,0,3000,4000,6000,0,3000,4000,4000,-1,-1,-1,0,1,2",
        "0,0,-1,-1,0,4000,5000,6000,0,3000,4000,4000,-1,-1,-1,1,0,2",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,2,0,2"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryNormalMapq1) {
    /*
     * Normal test case. There are multiple secondary regions which overlap one
     * supplementary region.
     * The mapping quality of these regions should be equal to 1.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 6000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
        {5000, 6000, 3000, 4000, score, 4, false},
        {5000, 6000, 3000, 4000, score, 5, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,6000,0,0,1000,4000,-1,-1,-1,0,0,32",
        "0,0,-1,-1,0,1000,2000,6000,0,1000,2000,4000,-1,-1,-1,0,1,32",
        "0,0,-1,-1,0,3000,4000,6000,0,3000,4000,4000,-1,-1,-1,0,1,1",
        "0,0,-1,-1,0,4000,5000,6000,0,3000,4000,4000,-1,-1,-1,1,0,1",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,2,0,1",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,3,0,1"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryNormalMapq0) {
    /*
     * Normal test case. There are multiple secondary regions which overlap one
     * supplementary region.
     * The mapping quality of these regions should be equal to 0.
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 6000;
    int32_t tlen = 4000;
    int32_t score = -1;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, score, 0, false},
        {1000, 2000, 1000, 2000, score, 1, false},
        {3000, 4000, 3000, 4000, score, 2, false},
        {4000, 5000, 3000, 4000, score, 3, false},
        {5000, 6000, 3000, 4000, score, 4, false},
        {5000, 6000, 3000, 4000, score, 5, false},
        {5000, 6000, 3000, 4000, score, 6, false},
        {5000, 6000, 3000, 4000, score, 7, false},
        {5000, 6000, 3000, 4000, score, 8, false},
        {5000, 6000, 3000, 4000, score, 9, false},
        {5000, 6000, 3000, 4000, score, 10, false},
        {5000, 6000, 3000, 4000, score, 11, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,-1,-1,0,0,1000,6000,0,0,1000,4000,-1,-1,-1,0,0,32",
        "0,0,-1,-1,0,1000,2000,6000,0,1000,2000,4000,-1,-1,-1,0,1,32",
        "0,0,-1,-1,0,3000,4000,6000,0,3000,4000,4000,-1,-1,-1,0,1,0",
        "0,0,-1,-1,0,4000,5000,6000,0,3000,4000,4000,-1,-1,-1,1,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,2,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,3,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,4,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,5,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,6,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,7,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,8,0,0",
        "0,0,-1,-1,0,5000,6000,6000,0,3000,4000,4000,-1,-1,-1,9,0,0",
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}

TEST(MappingResultsCommon, RelabelSupplementaryMinSecToPrimRatio) {
    /*
     * There is a secondary mapping, but it has a lower score than the overlapping primary
     * region. It's below the allowed ratio (min_sec_to_prim_ratio), so it should be
     * marked as filtered (region_priority = -1).
    */

    // Input values.
    double min_sec_to_prim_ratio = 0.80;
    int32_t allowed_overlap_bp = 0;
    int32_t qlen = 5000;
    int32_t tlen = 4000;
    // Tuple: qstart, qend, tstart, tend, score, region_priority, region_is_supplementary
    std::vector<std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, bool>> input_regions = {
        {0, 1000, 0, 1000, 1000, 0, false},
        {1000, 2000, 1000, 2000, 1000, 1, false},
        {3000, 4000, 3000, 4000, 1000, 2, false},
        {4000, 5000, 3000, 4000, 700, 3, false},
    };

    // Expected.
    // Columns:
    //      Aid,Bid,score,NumSeeds,Arev,Astart,Aend,Alen,Brev,Bstart,Bend,Blen,CovBasesA,CovBasesB,NumSeeds,RegionPriority,RegionIsSupplementary,Mapq
    std::vector<std::string> expected = {
        "0,0,1000,-1,0,0,1000,5000,0,0,1000,4000,-1,-1,-1,0,0,47",
        "0,0,1000,-1,0,1000,2000,5000,0,1000,2000,4000,-1,-1,-1,0,1,47",
        "0,0,1000,-1,0,3000,4000,5000,0,3000,4000,4000,-1,-1,-1,0,1,47",
        "0,0,700,-1,0,4000,5000,5000,0,3000,4000,4000,-1,-1,-1,-1,0,0"
    };

    // This constructs the RegionMapped objects from input regions, and runs the RelabelSupplementary.
    auto alns = HelperRunRelabelSupplementary(input_regions, qlen, tlen, min_sec_to_prim_ratio, allowed_overlap_bp);

    // Collect results as strings for comparison.
    std::vector<std::string> results;
    for (size_t i = 0; i < alns.size(); ++i) {
        std::ostringstream oss;
        oss << alns[i]->WriteAsCSV(',') << "," << alns[i]->GetRegionPriority() << ","
            << alns[i]->GetRegionIsSupplementary() << "," << alns[i]->MappingQuality();
        results.emplace_back(oss.str());
    }

    EXPECT_EQ(expected, results);
}
