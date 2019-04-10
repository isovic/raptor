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

bool VERBOSE_DEBUG_QID_SG = false;

TEST(SplitSegmentGraph, CreateSplitSegmentGraph1) {
    /*
     * Simple empty case.
    */
    std::string test_name("CreateSplitSegmentGraph1");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "],\n"
        "\"eitems\": [\n"
        "],\n"
        "\"adj\": {\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(SplitSegmentGraph, CreateSplitSegmentGraph2) {
    /*
     * A single node, no edges.
    */
    std::string test_name("CreateSplitSegmentGraph2");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "1000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":0,\"te\":1000}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":1000}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "],\n"
        "\"adj\": {\n"
        "0: [],\n"
        "1: []\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(SplitSegmentGraph, CreateSplitSegmentGraph3) {
    /*
     * Single node with a circular edge.
    */
    std::string test_name("CreateSplitSegmentGraph3");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "1000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "E1", "ref1+", "ref1+", "1000$", "1000$", "0", "0", "*"},
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":0,\"te\":1000}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":1000}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 0, \"w\": 0, \"v_int_id\": 0, \"w_int_id\": 0, \"data\": []}]\n"
        "],\n"
        "\"adj\": {\n"
        "0: [0],\n"
        "1: []\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(SplitSegmentGraph, CreateSplitSegmentGraph4) {
    /*
     * Test the construction of the split segment graph from a regular
     * segment graph.
    */
    std::string test_name("CreateSplitSegmentGraph4");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "1000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "E1", "ref1+", "ref1+", "500", "500", "700", "700", "*"},
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":0,\"te\":500}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":500,\"te\":700}}],\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":700,\"te\":1000}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":1000}}]\n"
        /// The reverse edges are not present in the original SegmentGraph, so the following is not used here.
        // "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":300}}]\n"
        // "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":300,\"te\":500}}]\n"
        // "[\"nitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":500,\"te\":1000}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 0, \"w\": 1, \"v_int_id\": 0, \"w_int_id\": 1, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 1, \"w\": 2, \"v_int_id\": 1, \"w_int_id\": 2, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"v\": 0, \"w\": 2, \"v_int_id\": 0, \"w_int_id\": 2, \"data\": []}]\n"
        // "[\"eitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"v\": 3, \"w\": 4, \"v_int_id\": 3, \"w_int_id\": 4, \"data\": []}]\n"
        // "[\"eitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"v\": 5, \"w\": 4, \"v_int_id\": 3, \"w_int_id\": 5, \"data\": []}]\n"
        // "[\"eitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"v\": 5, \"w\": 4, \"v_int_id\": 4, \"w_int_id\": 5, \"data\": []}]\n"
        "],\n"
        "\"adj\": {\n"
        "0: [1,2],\n"
        "1: [2],\n"
        "2: [],\n"
        "3: []\n"
        // "3: [4,5],\n"
        // "4: [5],\n"
        // "5: []\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(SplitSegmentGraph, CreateSplitSegmentGraph5) {
    /*
     * Test the construction of the split segment graph from a regular
     * segment graph.
    */
    std::string test_name("CreateSplitSegmentGraph5");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "1000", "0"},
    };
    std::vector<std::vector<std::string>> edges_gfa2 = {
        {"E", "E1", "ref1+", "ref1+", "500", "500", "700", "700", "*"},
        {"E", "E2", "ref1-", "ref1-", "700", "700", "500", "500", "*"},
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    // std::cerr << "Testing JSON output, seg_graph:\n" << seg_graph->ToJSON() << "\n";
    // std::cerr << "Testing JSON output, ssg:\n" << ssg->ToJSON() << "\n";

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":0,\"te\":500}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":500,\"te\":700}}],\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":700,\"te\":1000}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":300}}],\n"
        "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":300,\"te\":500}}],\n"
        "[\"nitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":500,\"te\":1000}}]\n"
        "],\n"
        "\"eitems\": [\n"
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 0, \"w\": 1, \"v_int_id\": 0, \"w_int_id\": 1, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 1, \"w\": 2, \"v_int_id\": 1, \"w_int_id\": 2, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"v\": 3, \"w\": 4, \"v_int_id\": 3, \"w_int_id\": 4, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"v\": 4, \"w\": 5, \"v_int_id\": 4, \"w_int_id\": 5, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"v\": 0, \"w\": 2, \"v_int_id\": 0, \"w_int_id\": 2, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"v\": 3, \"w\": 5, \"v_int_id\": 3, \"w_int_id\": 5, \"data\": []}]\n"
        "],\n"
        "\"adj\": {\n"
        "0: [1,2],\n"
        "1: [2],\n"
        "2: [],\n"
        "3: [4,5],\n"
        "4: [5],\n"
        "5: []\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

TEST(SplitSegmentGraph, CreateSplitSegmentGraph6) {
    /*
     * Test the construction of the split segment graph from a regular
     * segment graph.
    */
    std::string test_name("CreateSplitSegmentGraph6");
	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(9);

    // Define the test inputs.
    std::vector<std::vector<std::string>> nodes_str = {
        {"ref1", "695", "0"},
        {"snp1", "1", "1"},
        {"snp2", "1", "2"},
        {"snp3", "1", "3"},
        {"sv1", "5", "4"},
        {"sv3", "10", "5"},
        {"ref2", "100", "6"},
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

        {"E", "E7", "ref2-", "ref2-", "10", "10", "3", "3", "*"},   // Edgs always have fwd coordinates.
        {"E", "E8", "ref2+", "ref2-", "71", "71", "57", "57", "*"},
    };
    raptor::SegmentGraphPtr seg_graph = WrapConstructSegmentGraph(test_name,
                                                nodes_str,
                                                edges_gfa2,
                                                VERBOSE_DEBUG_QID_SG);

    raptor::SplitSegmentGraphPtr ssg = raptor::createSplitSegmentGraph(seg_graph);

    // std::cerr << "Testing JSON output, seg_graph:\n" << seg_graph->ToJSON() << "\n";
    // std::cerr << "Testing JSON output, ssg:\n" << ssg->ToJSON() << "\n";

    std::string result_str = ssg->ToJSON();

    std::string expected_str =
        "[\"graph\", {\n"
        "\"nitems\": [\n"
        // Fwd of ref1.
        "[\"nitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":0,\"te\":124}}],\n"
        "[\"nitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":124,\"te\":125}}],\n"
        "[\"nitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":125,\"te\":211}}],\n"
        "[\"nitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":211,\"te\":212}}],\n"
        "[\"nitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":212,\"te\":291}}],\n"
        "[\"nitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":291,\"te\":292}}],\n"
        "[\"nitem\", {\"name\": 6, \"removed\": 0, \"int_id\": 6, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":292,\"te\":464}}],\n"
        "[\"nitem\", {\"name\": 7, \"removed\": 0, \"int_id\": 7, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":464,\"te\":596}}],\n"
        "[\"nitem\", {\"name\": 8, \"removed\": 0, \"int_id\": 8, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":596,\"te\":605}}],\n"
        "[\"nitem\", {\"name\": 9, \"removed\": 0, \"int_id\": 9, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":0,\"ts\":605,\"te\":695}}],\n"
        // Rev of ref1, there are no edges, so just one node.
        "[\"nitem\", {\"name\": 10, \"removed\": 0, \"int_id\": 10, \"data\": {\"seg_name\":0,\"tid\":0,\"tr\":1,\"ts\":0,\"te\":695}}],\n"

        // Nodes for snp1.
        "[\"nitem\", {\"name\": 11, \"removed\": 0, \"int_id\": 11, \"data\": {\"seg_name\":1,\"tid\":1,\"tr\":0,\"ts\":0,\"te\":1}}],\n"
        "[\"nitem\", {\"name\": 12, \"removed\": 0, \"int_id\": 12, \"data\": {\"seg_name\":1,\"tid\":1,\"tr\":1,\"ts\":0,\"te\":1}}],\n"

        // Nodes for snp2.
        "[\"nitem\", {\"name\": 13, \"removed\": 0, \"int_id\": 13, \"data\": {\"seg_name\":2,\"tid\":2,\"tr\":0,\"ts\":0,\"te\":1}}],\n"
        "[\"nitem\", {\"name\": 14, \"removed\": 0, \"int_id\": 14, \"data\": {\"seg_name\":2,\"tid\":2,\"tr\":1,\"ts\":0,\"te\":1}}],\n"

        // Nodes for snp3.
        "[\"nitem\", {\"name\": 15, \"removed\": 0, \"int_id\": 15, \"data\": {\"seg_name\":3,\"tid\":3,\"tr\":0,\"ts\":0,\"te\":1}}],\n"
        "[\"nitem\", {\"name\": 16, \"removed\": 0, \"int_id\": 16, \"data\": {\"seg_name\":3,\"tid\":3,\"tr\":1,\"ts\":0,\"te\":1}}],\n"

        // Nodes for sv1.
        "[\"nitem\", {\"name\": 17, \"removed\": 0, \"int_id\": 17, \"data\": {\"seg_name\":4,\"tid\":4,\"tr\":0,\"ts\":0,\"te\":5}}],\n"
        "[\"nitem\", {\"name\": 18, \"removed\": 0, \"int_id\": 18, \"data\": {\"seg_name\":4,\"tid\":4,\"tr\":1,\"ts\":0,\"te\":5}}],\n"

        // Nodes for sv3.
        "[\"nitem\", {\"name\": 19, \"removed\": 0, \"int_id\": 19, \"data\": {\"seg_name\":5,\"tid\":5,\"tr\":0,\"ts\":0,\"te\":10}}],\n"
        "[\"nitem\", {\"name\": 20, \"removed\": 0, \"int_id\": 20, \"data\": {\"seg_name\":5,\"tid\":5,\"tr\":1,\"ts\":0,\"te\":10}}],\n"

        // Nodes for ref2.
        "[\"nitem\", {\"name\": 21, \"removed\": 0, \"int_id\": 21, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":0,\"ts\":0,\"te\":71}}],\n"
        "[\"nitem\", {\"name\": 22, \"removed\": 0, \"int_id\": 22, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":0,\"ts\":71,\"te\":100}}],\n"
        "[\"nitem\", {\"name\": 23, \"removed\": 0, \"int_id\": 23, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":1,\"ts\":0,\"te\":43}}],\n"
        "[\"nitem\", {\"name\": 24, \"removed\": 0, \"int_id\": 24, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":1,\"ts\":43,\"te\":90}}],\n"
        "[\"nitem\", {\"name\": 25, \"removed\": 0, \"int_id\": 25, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":1,\"ts\":90,\"te\":97}}],\n"
        "[\"nitem\", {\"name\": 26, \"removed\": 0, \"int_id\": 26, \"data\": {\"seg_name\":6,\"tid\":6,\"tr\":1,\"ts\":97,\"te\":100}}]\n"

        "],\n"
        "\"eitems\": [\n"
        // First, the implicit edges on ref1.
        "[\"eitem\", {\"name\": 0, \"removed\": 0, \"int_id\": 0, \"v\": 0, \"w\": 1, \"v_int_id\": 0, \"w_int_id\": 1, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 1, \"removed\": 0, \"int_id\": 1, \"v\": 1, \"w\": 2, \"v_int_id\": 1, \"w_int_id\": 2, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 2, \"removed\": 0, \"int_id\": 2, \"v\": 2, \"w\": 3, \"v_int_id\": 2, \"w_int_id\": 3, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 3, \"removed\": 0, \"int_id\": 3, \"v\": 3, \"w\": 4, \"v_int_id\": 3, \"w_int_id\": 4, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 4, \"removed\": 0, \"int_id\": 4, \"v\": 4, \"w\": 5, \"v_int_id\": 4, \"w_int_id\": 5, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 5, \"removed\": 0, \"int_id\": 5, \"v\": 5, \"w\": 6, \"v_int_id\": 5, \"w_int_id\": 6, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 6, \"removed\": 0, \"int_id\": 6, \"v\": 6, \"w\": 7, \"v_int_id\": 6, \"w_int_id\": 7, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 7, \"removed\": 0, \"int_id\": 7, \"v\": 7, \"w\": 8, \"v_int_id\": 7, \"w_int_id\": 8, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 8, \"removed\": 0, \"int_id\": 8, \"v\": 8, \"w\": 9, \"v_int_id\": 8, \"w_int_id\": 9, \"data\": []}],\n"
        // Implicit edges on ref2.
        "[\"eitem\", {\"name\": 9, \"removed\": 0, \"int_id\": 9, \"v\": 21, \"w\": 22, \"v_int_id\": 21, \"w_int_id\": 22, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 10, \"removed\": 0, \"int_id\": 10, \"v\": 23, \"w\": 24, \"v_int_id\": 23, \"w_int_id\": 24, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 11, \"removed\": 0, \"int_id\": 11, \"v\": 24, \"w\": 25, \"v_int_id\": 24, \"w_int_id\": 25, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 12, \"removed\": 0, \"int_id\": 12, \"v\": 25, \"w\": 26, \"v_int_id\": 25, \"w_int_id\": 26, \"data\": []}],\n"

        // Next, all of the explicit segment edges.

        // In-and-out into snp1.
        "[\"eitem\", {\"name\": 13, \"removed\": 0, \"int_id\": 13, \"v\": 0, \"w\": 11, \"v_int_id\": 0, \"w_int_id\": 11, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 14, \"removed\": 0, \"int_id\": 14, \"v\": 11, \"w\": 2, \"v_int_id\": 11, \"w_int_id\": 2, \"data\": []}],\n"

        // In-and-out into snp2.
        "[\"eitem\", {\"name\": 15, \"removed\": 0, \"int_id\": 15, \"v\": 2, \"w\": 13, \"v_int_id\": 2, \"w_int_id\": 13, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 16, \"removed\": 0, \"int_id\": 16, \"v\": 13, \"w\": 4, \"v_int_id\": 13, \"w_int_id\": 4, \"data\": []}],\n"

        // In-and-out into snp3.
        "[\"eitem\", {\"name\": 17, \"removed\": 0, \"int_id\": 17, \"v\": 4, \"w\": 15, \"v_int_id\": 4, \"w_int_id\": 15, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 18, \"removed\": 0, \"int_id\": 18, \"v\": 15, \"w\": 6, \"v_int_id\": 15, \"w_int_id\": 6, \"data\": []}],\n"

        // In-and-out into sv1.
        "[\"eitem\", {\"name\": 19, \"removed\": 0, \"int_id\": 19, \"v\": 6, \"w\": 17, \"v_int_id\": 6, \"w_int_id\": 17, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 20, \"removed\": 0, \"int_id\": 20, \"v\": 17, \"w\": 7, \"v_int_id\": 17, \"w_int_id\": 7, \"data\": []}],\n"

        // The deletion edge.
        "[\"eitem\", {\"name\": 21, \"removed\": 0, \"int_id\": 21, \"v\": 7, \"w\": 9, \"v_int_id\": 7, \"w_int_id\": 9, \"data\": []}],\n"

        // The extra edge which causes multiple forks for some nodes.
        "[\"eitem\", {\"name\": 22, \"removed\": 0, \"int_id\": 22, \"v\": 6, \"w\": 8, \"v_int_id\": 6, \"w_int_id\": 8, \"data\": []}],\n"

        // The additional reverse edges between ref2 and ref1-ref2.
        "[\"eitem\", {\"name\": 23, \"removed\": 0, \"int_id\": 23, \"v\": 24, \"w\": 26, \"v_int_id\": 24, \"w_int_id\": 26, \"data\": []}],\n"
        "[\"eitem\", {\"name\": 24, \"removed\": 0, \"int_id\": 24, \"v\": 21, \"w\": 24, \"v_int_id\": 21, \"w_int_id\": 24, \"data\": []}]\n"

        "],\n"
        "\"adj\": {\n"
        "0: [1,11],\n"
        "1: [2],\n"
        "2: [3,13],\n"
        "3: [4],\n"
        "4: [5,15],\n"
        "5: [6],\n"
        "6: [7,17,8],\n"
        "7: [8,9],\n"
        "8: [9],\n"
        "9: [],\n"
        "10: [],\n"
        "11: [2],\n"
        "12: [],\n"
        "13: [4],\n"
        "14: [],\n"
        "15: [6],\n"
        "16: [],\n"
        "17: [7],\n"
        "18: [],\n"
        "19: [],\n"
        "20: [],\n"
        "21: [22,24],\n"
        "22: [],\n"
        "23: [24],\n"
        "24: [25,26],\n"
        "25: [26],\n"
        "26: []\n"
        "}\n"
        "}\n"
        "]\n";

	/// LogSystem::GetInstance().SetProgramVerboseLevelFromInt(0);

    ASSERT_EQ(result_str, expected_str);
}

}
}
