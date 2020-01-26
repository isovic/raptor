#include <containers/mapping_result/linear_mapping_result.h>
#include <containers/mapping_result/mapping_result_common.h>
#include <tuple>

namespace raptor {

std::shared_ptr<raptor::LinearMappingResult> createMappingResult(int64_t qseq_id,
                                                       int64_t qseq_len,
                                                       std::string qseq_header,
                                                       mindex::IndexPtr index) {
    return std::shared_ptr<raptor::LinearMappingResult>(new LinearMappingResult(qseq_id, qseq_len, qseq_header, index));
}

LinearMappingResult::LinearMappingResult(int64_t _qseq_id,
                             int64_t _qseq_len,
                             std::string _qseq_header,
                             mindex::IndexPtr _index) :
                                      qseq_id_(_qseq_id),
                                      qseq_len_(_qseq_len),
                                      qseq_header_(_qseq_header),
                                      index_(_index) {
}

std::vector<std::shared_ptr<raptor::RegionBase>> LinearMappingResult::CollectRegions(bool one_hit_per_target, bool do_relabel_sec_supp) const {
    int32_t num_paths = static_cast<int32_t>(target_anchors_.size());
    std::vector<std::pair<int64_t, size_t>> path_scores;
    for (size_t path_id = 0; path_id < target_anchors_.size(); ++path_id) {
        const auto& path_aln = target_anchors_[path_id];
        path_scores.emplace_back(std::make_pair(path_aln->score(), path_id));
    }
    std::sort(path_scores.begin(), path_scores.end());
    std::reverse(path_scores.begin(), path_scores.end());

    std::vector<std::shared_ptr<raptor::RegionBase>> ret;

    // Used for filtering multiple hits to the same target.
    // Also, label the alignment priorities and path IDs.
    std::unordered_map<std::string, int32_t> query_target_pairs;
    for (int64_t i = 0; i < static_cast<int64_t>(path_scores.size()); ++i) {
        auto path_score = std::get<0>(path_scores[i]);
        auto path_id = std::get<1>(path_scores[i]);
        auto& tanchors = target_anchors_[path_id];
        for (auto& aln: tanchors->hits()) {
            aln->SetRegionPriority(i);
            aln->SetRegionIsSupplementary(0);
            if (one_hit_per_target) {
                std::string pair_name = std::to_string(aln->QueryID()) + std::string("->") + std::to_string(aln->TargetID());
                if (query_target_pairs.find(pair_name) != query_target_pairs.end()) {
                  continue;
                }
                query_target_pairs[pair_name] = 1;
            }
            aln->path_id(static_cast<int32_t>(path_id));
            aln->num_paths(num_paths);
            aln->segment_id(static_cast<int32_t>(0));
            aln->num_segments(1);
            ret.emplace_back(aln);
        }
    }

    if (do_relabel_sec_supp) {
        raptor::RelabelSupplementary(ret, 0.0);
    }

    return ret;
}

std::string LinearMappingResult::Verbose() const {
    std::ostringstream oss;

    oss << "[LinearMappingResult::Verbose()] qseq_id = " << qseq_id_ << ", qseq_len = "
        << qseq_len_ << ", qseq_header = '" << qseq_header_
        << "'\n";

    oss << "target_anchors:\n";

    for (size_t i = 0; i < target_anchors_.size(); ++i) {
        const auto& single_target_anchors = target_anchors_[i];
        oss << "[" << i << "] TargetAnchors for target_id = " << single_target_anchors->env()->t_id << ":\n";
        oss << single_target_anchors->VerbosePointers() << "\n";
    }

    oss << "target_hits:\n";

    for (size_t i = 0; i < target_hits_.size(); ++i) {
        const auto& single_target_hits = target_hits_[i];
        oss << "[" << i << "] TargetHits for target_id = " << single_target_hits->env()->t_id << ":\n";
        oss << single_target_hits->Verbose() << "\n";
    }

    return oss.str();
}

std::string LinearMappingResult::WriteAsCSV(const char separator) const {
    std::ostringstream oss;

    // oss << "[LinearMappingResult::Verbose()] qseq_id = " << qseq_id_ << ", qseq_len = "
    //     << qseq_len_ << ", qseq_header = '" << qseq_header_
    //     << "'\n";

    // oss << "target_hits:\n";
    for (size_t i = 0; i < target_hits_.size(); ++i) {
        const auto& single_target_hits = target_hits_[i];
        for (const auto& hit: single_target_hits->hits()) {
            oss << "H" << separator << i << separator << hit.WriteAsCSV(separator) << "\n";
        }
    }

    // oss << "target_anchors:\n";
    for (size_t i = 0; i < target_anchors_.size(); ++i) {
        const auto& single_target_anchors = target_anchors_[i];
        for (const auto& hit: single_target_anchors->hits()) {
            oss << "A" << separator << i << separator << hit->WriteAsCSV(separator) << "\n";
        }
    }

    return oss.str();
}

int64_t LinearMappingResult::QueryId() const {
    return qseq_id_;
}

int64_t LinearMappingResult::QueryLen() const {
    return qseq_len_;
}

std::string LinearMappingResult::QueryHeader() const {
    return qseq_header_;
}

const mindex::IndexPtr LinearMappingResult::Index() const {
    return index_;
}

MapperReturnValueBase LinearMappingResult::ReturnValue() const {
    return return_value_;
}

const std::unordered_map<std::string, double>& LinearMappingResult::Timings() const {
    return timings_;
}

}
