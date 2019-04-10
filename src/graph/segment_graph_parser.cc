#include <graph/segment_graph_parser.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <utility/trim.hpp>
#include <utility/stringutil.h>
#include <aligner/cigar.h>

namespace raptor {

std::shared_ptr<raptor::SegmentGraph> GraphLoader::FromGFA(const std::string& graph_path, mindex::IndexPtr index, bool add_symmetric_arcs) {
    // Find the GFA version.
    std::ifstream ifs(graph_path);
    if (ifs.is_open() == false) {
        std::cerr << "ERROR: Can not open file " << graph_path << " for reading! Exiting." << std::endl;
        exit(1);
    }
    float gfa_version = GraphLoader::FindGFAVersion_(ifs);
    if (gfa_version < 0) {
        std::cerr << "ERROR: Version header not found in graph file '" << graph_path << "'. Exiting." << std::endl;
        exit(1);
    }

    // Setup an empty graph.
    std::shared_ptr<raptor::SegmentGraph> graph = createSegmentGraph();

    // Initialize the nodes from the index.
    GraphLoader::AddGraphNodesFromIndex(index, graph);

    if (gfa_version >= 1.0 && gfa_version < 2.0) {
        LoadEdgesFromGFA1_(ifs, index, graph);

    } else if (gfa_version >= 2.0 && gfa_version < 3.0) {
        LoadEdgesFromGFA2_(ifs, index, graph);

    } else {
        std::cerr << "ERROR: GFA version " << gfa_version << " is not currently supported. Exiting." << std::endl;
        exit(1);
    }

    ifs.close();

    AssertGraphValidity(index, graph);

    if (add_symmetric_arcs) {
        AddSymmetricArcs_(graph);
    }

    graph->BuildSegmentTrees();

    return graph;
}

void GraphLoader::AssertGraphValidity(const mindex::IndexPtr& index, const std::shared_ptr<raptor::SegmentGraph>& graph) {
    for (const auto& edge_item: graph->edges()) {
        const auto& edge = edge_item->data();
        int64_t source_id = edge->source_id();
        int64_t sink_id = edge->sink_id();

        const mindex::SequencePtr& source_seq = index->seqs()->GetSeqByAbsID(source_id);
        if (source_seq == nullptr) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Absolute source sequence ID %ld not found in index!", source_id);
        }

        const mindex::SequencePtr& sink_seq = index->seqs()->GetSeqByAbsID(sink_id);
        if (sink_seq == nullptr) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Absolute sink sequence ID %ld not found in index!", sink_id);
        }

        if (edge->source_start() < 0) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The source_start value of edge '%s' is < 0. source_start = %ld, source_id = %ld.", edge->symbolic_edge_name().c_str(), edge->source_start(), source_id);
        }
        if (edge->source_end() < 0) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The source_end value of edge '%s' is < 0. source_start = %ld, source_id = %ld.", edge->symbolic_edge_name().c_str(), edge->source_end(), source_id);
        }
        if (edge->source_end() < edge->source_start()) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "source_end < source_start for edge '%s'. source_start = %ld, source_end = %ld, source_id = %ld.", edge->symbolic_edge_name().c_str(), edge->source_start(), edge->source_end(), source_id);
        }
        if (edge->source_end() > source_seq->len()) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The source_end value of edge '%s' larger than the indexed sequence size. source_end = %ld, source_seq->len() = %ld, source_id = %ld.", edge->symbolic_edge_name().c_str(), edge->source_end(), source_seq->len(), source_id);
        }

        if (edge->sink_start() < 0) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The sink_start value of edge '%s' is < 0. sink_start = %ld, sink_id = %ld.", edge->symbolic_edge_name().c_str(), edge->sink_start(), sink_id);
        }
        if (edge->sink_end() < 0) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The sink_end value of edge '%s' is < 0. sink_start = %ld, sink_id = %ld.", edge->symbolic_edge_name().c_str(), edge->sink_end(), sink_id);
        }
        if (edge->sink_end() < edge->sink_start()) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "sink_end < sink_start for edge '%s'. sink_start = %ld, sink_end = %ld, sink_id = %ld.", edge->symbolic_edge_name().c_str(), edge->sink_start(), edge->sink_end(), sink_id);
        }
        if (edge->sink_end() > sink_seq->len()) {
            FATAL_REPORT(ERR_UNEXPECTED_VALUE, "The sink_end value of edge '%s' larger than the indexed sequence size. sink_end = %ld, sink_seq->len() = %ld, sink_id = %ld.", edge->symbolic_edge_name().c_str(), edge->sink_end(), sink_seq->len(), sink_id);
        }
    }
}

void GraphLoader::AddGraphNodesFromIndex(mindex::IndexPtr index, std::shared_ptr<raptor::SegmentGraph> graph) {
    // Add sequences from the index as nodes to the graph.
    for (size_t i = 0; i < index->num_seqs(); ++i) {
        const mindex::SequencePtr& seq = index->seqs()->seqs()[i];
        if (seq == nullptr) {
            std::cerr << "[GraphLoader::AddGraphNodesFromIndex] Warning: seq_id does not match a sequence in index! Skipping hit.\n";
            continue;
        }

        auto node = createSegmentNode(seq->GetTrimmedHeader(), seq->len(), false, seq->abs_id());
        node->internal_id(graph->nodes().size());

        graph->AddNode(seq->abs_id(), node);
    }
}

float GraphLoader::FindGFAVersion_(std::ifstream& ifs) {
    // Determine the GFA version.
    float gfa_version = -1.0;
    std::string line;
    while(std::getline(ifs, line)) {
        if (line.size() == 0) { continue; }
        // Tokenize the line.
        std::istringstream iss(line);
        std::vector<std::string> params;
        std::string param;
        while(iss >> param) {
            params.emplace_back(param);
        }
        if (params[0] != std::string("H")) {
            break;
        }
        for (auto& p: params) {
            auto tokens = raptor::Tokenize(p, ':');
            if (tokens.size() != 3) {
                continue;
            }
            if (tokens[0] == std::string("VN")) {
                std::string::size_type sz;
                gfa_version = std::stof(tokens.back(), &sz);
                break;
            }
        }
    }
    return gfa_version;
}

void GraphLoader::AddSymmetricArcs_(std::shared_ptr<raptor::SegmentGraph> graph) {
    std::vector<SegmentEdgePtr> new_edges;

    for (auto& edge_item: graph->edges()) {
        auto& orig_edge = edge_item->data();
        auto edge = createSegmentEdge();

        edge->symbolic_edge_name(edge_item->data()->symbolic_edge_name() + std::string("_rev"));
        edge->id(graph->edges().size());

        edge->source_name(orig_edge->sink_name());
        edge->source_id(orig_edge->sink_id());
        edge->source_is_rev(!orig_edge->sink_is_rev());
        edge->source_start(orig_edge->sink_len() - orig_edge->sink_end());
        edge->source_end(orig_edge->sink_len() - orig_edge->sink_start());
        edge->source_len(orig_edge->sink_len());

        edge->sink_name(orig_edge->source_name());
        edge->sink_id(orig_edge->source_id());
        edge->sink_is_rev(!orig_edge->source_is_rev());
        edge->sink_start(orig_edge->source_len() - orig_edge->source_end());
        edge->sink_end(orig_edge->source_len() - orig_edge->source_start());
        edge->sink_len(orig_edge->source_len());

        edge->alignment(orig_edge->alignment());    // TODO: Check if this needs to be reversed if both sequences were reversed.

        // Check if the same edge already exists. If so, perhaps
        // the input graph is already symmetric.
        auto existing_edges = graph->GetEdges(edge->source_id(), edge->sink_id());
        bool symmetric_edge_exists = false;
        for (auto& existing_edge: existing_edges) {
            if (existing_edge->data()->source_is_rev() == edge->source_is_rev() &&
                    existing_edge->data()->sink_is_rev() == edge->sink_is_rev() &&
                    existing_edge->data()->source_start() == edge->source_start() &&
                    existing_edge->data()->source_end() == edge->source_end() &&
                    existing_edge->data()->sink_start() == edge->sink_start() &&
                    existing_edge->data()->sink_end() == edge->sink_end()) {
                symmetric_edge_exists = true;
                break;
            }
        }

        if (symmetric_edge_exists) {
            continue;
        }

        new_edges.emplace_back(edge);
    }

    for (auto& edge: new_edges) {
        graph->AddEdge(edge->source_id(), edge->sink_id(), edge);
    }
}

void GraphLoader::LoadEdgesFromGFA1_(
                                std::ifstream& ifs,
                                mindex::IndexPtr index,
                                std::shared_ptr<raptor::SegmentGraph> graph) {
    std::string line;
    while(std::getline(ifs, line)) {
        if (line.size() == 0) {
            continue;
        }

        // Tokenize the line.
        std::istringstream iss(line);
        std::vector<std::string> params;
        std::string param;
        while(iss >> param) {
            params.emplace_back(param);
        }

        if (params.size() == 0) {
            continue;
        }

        if (params[0] == "#") {
            // Comment lines.

        } else if (params[0] == "H") {

        } else if (params[0] == "S") {

        } else if (params[0] == "C") {

        } else if (params[0] == "P") {

        } else if (params[0] == "L") {

            raptor::SegmentEdgePtr edge = GraphLoader::ParseGFA1Edge(graph, params);
            if (edge != nullptr) {
                graph->AddEdge(edge->source_id(), edge->sink_id(), edge);
            }
        }
    }
}

void GraphLoader::LoadEdgesFromGFA2_(
                                std::ifstream& ifs,
                                mindex::IndexPtr index,
                                std::shared_ptr<raptor::SegmentGraph> graph) {
    std::string line;
    while(std::getline(ifs, line)) {
        if (line.size() == 0) {
            continue;
        }

        // Tokenize the line.
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string param;
        while(iss >> param) {
            tokens.emplace_back(param);
        }

        if (tokens.size() == 0) {
            continue;
        }

        if (tokens[0] == "H") {

        } else if (tokens[0] == "S") {

        } else if (tokens[0] == "E") {
            raptor::SegmentEdgePtr edge = GraphLoader::ParseGFA2Edge(graph, tokens);
            if (edge != nullptr) {
                graph->AddEdge(edge->source_id(), edge->sink_id(), edge);
            }
        }
    }
}

raptor::SegmentEdgePtr GraphLoader::ParseGFA1Edge(const std::shared_ptr<raptor::SegmentGraph>& graph,
                                                    const std::vector<std::string>& params) {
    auto edge = createSegmentEdge();

    std::string source_name = params[1];
    raptor::SegmentNodePtr source_node = graph->GetNodeByHeader(source_name);
    if (source_node == nullptr) {
        LOG_ALL("Warning: Edge skipped, because the source segment '%s' is not part of the index.\n", source_name.c_str());
        return nullptr;
    }
    uint64_t source_id = source_node->seq_id();
    bool source_is_rev = (params[2] == "+") ? false : true;

    std::string sink_name = params[3];
    raptor::SegmentNodePtr sink_node = graph->GetNodeByHeader(sink_name);
    if (sink_node == nullptr) {
        LOG_ALL("Warning: Edge skipped, because the sink segment '%s' is not part of the index.\n", sink_name.c_str());
        return nullptr;
    }
    uint64_t sink_id = sink_node->seq_id();
    bool sink_is_rev = (params[4] == "+") ? false : true;

    const std::string& aln_str = params[5];

    int64_t source_overlap_len = 0;
    int64_t target_overlap_len = 0;
    if (aln_str != std::string("*")) {
        auto cigar = raptor::CigarStringToVector(aln_str);
        source_overlap_len = raptor::QueryLengthFromCigar(cigar, true);
        target_overlap_len = raptor::ReferenceLengthFromCigar(cigar);
    }

    int64_t len1 = source_node->len();
    int64_t start1 = len1 - source_overlap_len;
    int64_t end1 = len1;

    int64_t start2 = 0;
    int64_t end2 = target_overlap_len;
    int64_t len2 = sink_node->len();

    // // In case this edge is connecting two reverse coordinates on the same sequence,
    // // make sure that the edge is reversed. Otherwise, we'll create cycles.
    // if (source_is_rev && sink_is_rev && source_name == sink_name) {
    //     edge->symbolic_edge_name("");       // TODO: There is an optional custom tag ID.
    //     edge->id(graph->edges().size());

    //     edge->source_name(sink_name);
    //     edge->source_id(sink_id);
    //     edge->source_is_rev(sink_is_rev);

    //     edge->sink_name(source_name);
    //     edge->sink_id(source_id);
    //     edge->sink_is_rev(source_is_rev);

    //     edge->alignment(aln_str);
    //     edge->source_len(len2);
    //     edge->source_start(start2);
    //     edge->source_end(end2);
    //     edge->sink_len(len1);
    //     edge->sink_start(start1);
    //     edge->sink_end(end1);

    // } else {
        edge->symbolic_edge_name("");       // TODO: There is an optional custom tag ID.
        edge->id(graph->edges().size());

        edge->source_name(source_name);
        edge->source_id(source_id);
        edge->source_is_rev(source_is_rev);

        edge->sink_name(sink_name);
        edge->sink_id(sink_id);
        edge->sink_is_rev(sink_is_rev);

        edge->alignment(aln_str);
        edge->source_len(len1);
        edge->source_start(start1);
        edge->source_end(end1);
        edge->sink_len(len2);
        edge->sink_start(start2);
        edge->sink_end(end2);
    // }

    return nullptr;
}

raptor::SegmentEdgePtr GraphLoader::ParseGFA2Edge(const std::shared_ptr<raptor::SegmentGraph>& graph,
                                                    const std::vector<std::string>& params) {
    raptor::SegmentEdgePtr edge = createSegmentEdge();

    std::string source_name = params[2].substr(0, params[2].size() - 1);
    raptor::SegmentNodePtr source_node = graph->GetNodeByHeader(source_name);
    if (source_node == nullptr) {
        LOG_ALL("Warning: Edge skipped, because the source segment '%s' is not part of the index.\n", source_name.c_str());
        return nullptr;
    }
    uint64_t source_id = source_node->seq_id();
    bool source_is_rev = (params[2].back() == '+') ? false : true;

    std::string sink_name = params[3].substr(0, params[3].size() - 1);
    raptor::SegmentNodePtr sink_node = graph->GetNodeByHeader(sink_name);
    if (sink_node == nullptr) {
        LOG_ALL("Warning: Edge skipped, because the sink segment '%s' is not part of the index.\n", sink_name.c_str());
        return nullptr;
    }
    uint64_t sink_id = sink_node->seq_id();
    bool sink_is_rev = (params[3].back() == '+') ? false : true;

    // Start and End positions are always in the fwd strand, according to the standard.
    int64_t start1 = (params[4].back() != '$') ? (atoi(params[4].c_str())) : (atoi(params[4].substr(0, params[4].size() - 1).c_str()));
    int64_t end1 = (params[5].back() != '$') ? (atoi(params[5].c_str())) : (atoi(params[5].substr(0, params[5].size() - 1).c_str()));
    int64_t len1 = source_node->len();

    // Start and End positions are always in the fwd strand, according to the standard.
    int64_t start2 = (params[6].back() != '$') ? (atoi(params[6].c_str())) : (atoi(params[6].substr(0, params[6].size() - 1).c_str()));
    int64_t end2 = (params[7].back() != '$') ? (atoi(params[7].c_str())) : (atoi(params[7].substr(0, params[7].size() - 1).c_str()));
    int64_t len2 = sink_node->len();

    const std::string& aln_str = params[8];

    // // In case this edge is connecting two reverse coordinates on the same sequence,
    // // make sure that the edge is reversed. Otherwise, we'll create cycles.
    // if (source_is_rev && sink_is_rev && source_name == sink_name) {
    //     edge->symbolic_edge_name(params[1]);
    //     edge->id(graph->edges().size());

    //     edge->source_name(sink_name);
    //     edge->source_id(sink_id);
    //     edge->source_is_rev(sink_is_rev);

    //     edge->sink_name(source_name);
    //     edge->sink_is_rev(source_is_rev);
    //     edge->sink_id(source_id);

    //     edge->alignment(aln_str);

    //     edge->source_len(len2);
    //     if (!edge->source_is_rev()) {
    //         edge->source_start(start2);
    //         edge->source_end(end2);
    //     } else {
    //         // Convert the coordinates to the strand of the sequence.
    //         edge->source_start(len2 - end2);
    //         edge->source_end(len2 - start2);
    //     }

    //     edge->sink_len(len1);
    //     if (!edge->sink_is_rev()) {
    //         edge->sink_start(start1);
    //         edge->sink_end(end1);
    //     } else {
    //         // Convert the coordinates to the strand of the sequence.
    //         edge->sink_start(len1 - end1);
    //         edge->sink_end(len1 - start1);
    //     }

    // } else {
        edge->symbolic_edge_name(params[1]);
        edge->id(graph->edges().size());

        edge->source_name(source_name);
        edge->source_id(source_id);
        edge->source_is_rev(source_is_rev);

        edge->sink_name(sink_name);
        edge->sink_is_rev(sink_is_rev);
        edge->sink_id(sink_id);

        edge->alignment(aln_str);

        edge->source_len(len1);
        if (!edge->source_is_rev()) {
            edge->source_start(start1);
            edge->source_end(end1);
        } else {
            // Convert the coordinates to the strand of the sequence.
            edge->source_start(len1 - end1);
            edge->source_end(len1 - start1);
        }

        edge->sink_len(len2);
        if (!edge->sink_is_rev()) {
            edge->sink_start(start2);
            edge->sink_end(end2);
        } else {
            // Convert the coordinates to the strand of the sequence.
            edge->sink_start(len2 - end2);
            edge->sink_end(len2 - start2);
        }
    // }

    return edge;
}

}
