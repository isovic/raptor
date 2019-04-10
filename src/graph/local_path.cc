/*
 * local_path.cc
 *
 *  Created on: Feb 1, 2018
 *      Author: Ivan Sovic
 */

#include <graph/local_path.h>
#include <sstream>

namespace raptor {

std::string LocalPathTools::SummarizeFullPathAsString(const std::shared_ptr<raptor::LocalPath> path) {
    std::ostringstream oss;

    if (path->nodes().size() == 0) {
        return oss.str();
    }

    auto& first_node = path->nodes().front()->data();

    oss << first_node->id();

    for (size_t node_id = 1 ; node_id < path->nodes().size(); node_id++) {
        size_t edge_id = node_id - 1;
        auto& edge = path->edges()[edge_id]->data();
        auto& node = path->nodes()[node_id]->data();

        if (edge->IsImplicit()) {
            oss << "(i)";
        } else {
            oss << "(";
            for (size_t i = 0; i < edge->segment_edges().size(); i++) {
                auto& seg_edge = edge->segment_edges()[i];
                // oss << "(" << seg_edge->id() << ")";
                if (i > 0) {
                    oss << "=";
                }
                oss << seg_edge->id();
            }
            oss << ")";
        }
        oss << node->id();
    }

    return oss.str();
}

std::string LocalPathTools::SummarizeJumpPathAsString(const std::shared_ptr<raptor::LocalPath> path) {
    std::ostringstream oss;

    if (path->nodes().size() == 0) {
        return oss.str();
    }

    auto& first_node = path->nodes().front()->data();
    oss << first_node->id();

    for (size_t node_id = 1 ; node_id < path->nodes().size(); node_id++) {
        size_t edge_id = node_id - 1;
        auto& edge = path->edges()[edge_id]->data();
        if (edge == nullptr) {
            LOG_ALL("Edge ptr is null! In SummarizeJumpPathAsString. Skipping.\n");
            break;
        }

        if (edge->IsImplicit() == false) {
            for (size_t i = 0; i < edge->segment_edges().size(); i++) {
                auto& seg_edge = edge->segment_edges()[i];
                if (seg_edge == nullptr) {
                    LOG_ALL("Seg edge ptr is null! In SummarizeJumpPathAsString. Skipping.\n");
                    break;
                }
                oss << "(" << seg_edge->id() << ")";
            }
        }
    }

    auto& last_node = path->nodes().back()->data();
    oss << last_node->id();

    return oss.str();
}

std::shared_ptr<raptor::LocalPath> LocalPathTools::MergeImplicitEdges(const std::shared_ptr<raptor::LocalPath> path) {

    if (path->nodes().size() == 0) {
        return nullptr;
    }

    int32_t num_chains = 0;

    auto& first_node = path->nodes().front()->data();
    // This creates a new anchor which wraps the entire span of
    // an implicit path, and sets the scores to match.
    // std::shared_ptr<raptor::RegionMapped> new_anchor(new raptor::RegionMapped(*(first_node->hits().front())));
    std::shared_ptr<raptor::RegionMapped> new_anchor(new raptor::RegionMapped(*(first_node)));
    new_anchor->qend(first_node->QueryEnd());
    new_anchor->tend(first_node->TargetEnd());
    new_anchor->cov_bases_q(first_node->cov_bases_q());
    new_anchor->cov_bases_t(first_node->cov_bases_t());
    new_anchor->num_seeds(first_node->num_seeds());
    new_anchor->score(first_node->score());

    auto merged_path = raptor::LocalPath::createPath(num_chains, new_anchor);
    num_chains += 1;

    for (size_t node_id = 1 ; node_id < path->nodes().size(); node_id++) {
        size_t edge_id = node_id - 1;
        auto& edge = path->edges()[edge_id]->data();
        auto& node = path->nodes()[node_id]->data();

        if (edge->IsImplicit()) {
            new_anchor->qend(node->QueryEnd());
            new_anchor->tend(node->TargetEnd());
            new_anchor->cov_bases_q(new_anchor->cov_bases_q() + node->cov_bases_q());
            new_anchor->cov_bases_t(new_anchor->cov_bases_t() + node->cov_bases_t());
            new_anchor->num_seeds(new_anchor->num_seeds() + node->num_seeds());
            new_anchor->score(new_anchor->score() + node->score());

        } else {
            new_anchor = std::shared_ptr<raptor::RegionMapped>(new raptor::RegionMapped(*(node)));
            new_anchor->qend(node->QueryEnd());
            new_anchor->tend(node->TargetEnd());
            new_anchor->cov_bases_q(node->cov_bases_q());
            new_anchor->cov_bases_t(node->cov_bases_t());
            new_anchor->num_seeds(node->num_seeds());
            new_anchor->score(node->score());

            merged_path->Add(num_chains, new_anchor, edge);

            num_chains += 1;
        }
    }

    return merged_path;
}

}
