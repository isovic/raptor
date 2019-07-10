#include <algorithm>
#include <raptor_fetch/overlaps/overlap_file.h>
#include <log/log_tools.h>
#include <utility/files.hpp>
#include <iostream>
#include <memory>

namespace raptor {

OverlapFilePtr createOverlapFile(const std::string& path) {
    auto ret = std::unique_ptr<raptor::OverlapFile>(new raptor::OverlapFile());
    ret->OpenFile(path);
    // ret->HashAllQnames();
    return ret;
}

void OverlapFile::OpenFile(const std::string& path) {
    fp_handlers_ = raptor::createOverlapFileHandlers(path);
}

// bool OverlapFile::HashAllQnames() {
//     if (fp_handlers_ == nullptr) {
//         WARNING_REPORT(ERR_OPENING_FILE, "No overlap file was opened.");
//         return false;
//     }

//     qname_map_.clear();

//     auto curr_pos = gztell(fp_handlers_->fp_gzip);
//     gzseek(fp_handlers_->fp_gzip, 0, SEEK_SET);

//     std::string line;
//     while (raptor::ReadGZLine(fp_handlers_->fp_gzip, line)) {
//         // std::unique_ptr<raptor::Overlap> overlap = raptor::ParseOverlap(line);
//         // std::cout << overlap->a_name() << " " << overlap->b_name() << "\n";
//     }

//     gzseek(fp_handlers_->fp_gzip, curr_pos, SEEK_SET);

//     return true;
// }

std::unique_ptr<raptor::Overlap> OverlapFile::YieldOverlap_(const std::string& line) const {

    std::unique_ptr<raptor::Overlap> overlap = raptor::ParseOverlap(line);

    if (overlap == nullptr) {
        return overlap;
    }

    int32_t a_id = -1, b_id = -1;

    auto it_a = qname_map_.find(overlap->a_name());
    if (it_a != qname_map_.end()) {
        a_id = it_a->second;
    }

    auto it_b = qname_map_.find(overlap->b_name());
    if (it_b != qname_map_.end()) {
        b_id = it_b->second;
    }

    overlap->a_id(a_id);
    overlap->b_id(b_id);

    return overlap;
}

// bool OverlapFile::HashAllQnames() {
//     if (fp_handlers_ == nullptr) {
//         WARNING_REPORT(ERR_OPENING_FILE, "No overlap file was opened.");
//         return false;
//     }

//     qname_map_.clear();
//     query_data_.clear();

//     auto curr_pos = ftell(fp_handlers_->fp);
//     fseek(fp_handlers_->fp, 0, SEEK_SET);

//     std::string line;

//     while (fp_handlers_->ReadLine(line)) {
//         auto ovl_parsed = YieldOverlap_(line);

//         int32_t a_id = -1, b_id = -1;

//         if (ovl_parsed->a_id() < 0) {
//             a_id = query_data_.size();
//             query_data_.emplace_back(raptor::QueryData(ovl_parsed->a_name(), ovl_parsed->a_len(), ovl_parsed->clip_a_start(), ovl_parsed->clip_a_end(), 0));
//             qname_map_[ovl_parsed->a_name()] = a_id;
//             ovl_parsed->a_id(a_id);
//         }

//         if (ovl_parsed->b_id() < 0) {
//             b_id = (ovl_parsed->b_name() == ovl_parsed->a_name()) ? (a_id) : query_data_.size();
//             query_data_.emplace_back(raptor::QueryData(ovl_parsed->b_name(), ovl_parsed->b_len(), ovl_parsed->clip_b_start(), ovl_parsed->clip_b_end(), 0));
//             qname_map_[ovl_parsed->b_name()] = b_id;
//             ovl_parsed->b_id(b_id);
//         }

//         // std::cout << ovl_parsed->a_name() << " " << ovl_parsed->b_name() << "\n";
//         // std::cout << "    " << query_data_[ovl_compact->a_id()].Verbose() << "\n";
//         // std::cout << "    " << query_data_[ovl_compact->b_id()].Verbose() << "\n";
//     }

//     fseek(fp_handlers_->fp, curr_pos, SEEK_SET);

//     return true;
// }

bool OverlapFile::LoadNextBatch(size_t batch_size, bool single_query_only, int32_t min_len, double min_score, double min_idt) {
    if (fp_handlers_ == nullptr) {
        WARNING_REPORT(ERR_OPENING_FILE, "No overlap file was opened.");
        return false;
    }

    batch_start_id_ += static_cast<int32_t>(batch_.size());
    batch_.clear();
    int32_t prev_a_id = -1;
    int32_t num_loaded = 0;

    std::string line;
    while (fp_handlers_->ReadLine(line)) {
        auto fp_pos = ftell(fp_handlers_->fp);
        auto ovl_parsed = YieldOverlap_(line);

        if (ovl_parsed == nullptr) {
            continue;
        }

        ++num_loaded;

        int32_t a_id = ovl_parsed->a_id(), b_id = ovl_parsed->b_id();

        if (a_id < 0) {
            a_id = static_cast<int32_t>(query_data_.size());
            auto new_query_data = std::make_unique<raptor::QueryData>(ovl_parsed->a_name(), ovl_parsed->a_len(), ovl_parsed->clip_a_start(), ovl_parsed->clip_a_end(), 0);
            query_data_.emplace_back(std::move(new_query_data));
            qname_map_[ovl_parsed->a_name()] = a_id;
            ovl_parsed->a_id(a_id);
        }

        if (b_id < 0) {
            if (ovl_parsed->b_name() == ovl_parsed->a_name()) {
                b_id = a_id;
            } else {
                b_id = static_cast<int32_t>(query_data_.size());
                auto new_query_data = std::make_unique<raptor::QueryData>(ovl_parsed->b_name(), ovl_parsed->b_len(), ovl_parsed->clip_b_start(), ovl_parsed->clip_b_end(), 0);
                query_data_.emplace_back(std::move(new_query_data));
            }
            // query_data_.emplace_back(raptor::QueryData(ovl_parsed->b_name(), ovl_parsed->b_len(), ovl_parsed->clip_b_start(), ovl_parsed->clip_b_end(), 0));
            qname_map_[ovl_parsed->b_name()] = b_id;
            ovl_parsed->b_id(b_id);
        }

        if (single_query_only && num_loaded > 1 && ovl_parsed->a_id() != prev_a_id) {
            // Rewind to before the current query.
            fseek(fp_handlers_->fp, fp_pos, SEEK_SET);
            break;
        }
        prev_a_id = a_id;

        std::unique_ptr<raptor::OverlapCompact> ovl_compact = raptor::createOverlapCompact(
            ovl_parsed->a_id(), ovl_parsed->b_id(), ovl_parsed->score(), ovl_parsed->identity(),
            ovl_parsed->a_rev(), ovl_parsed->a_start(), ovl_parsed->a_end(),
            ovl_parsed->b_rev(), ovl_parsed->b_start(), ovl_parsed->b_end()
            );

        if (ovl_parsed->a_len() < min_len || ovl_parsed->b_len() < min_len) {
            ovl_compact->SetFlagLowLen(true);
        }
        if (ovl_parsed->a_span() < min_len || ovl_parsed->b_span() < min_len) {
            ovl_compact->SetFlagLowSpan(true);
        }
        if (abs(ovl_parsed->score()) < min_score) {
            ovl_compact->SetFlagLowScore(true);
        }
        if (ovl_parsed->identity() < min_idt) {
            ovl_compact->SetFlagLowIdentity(true);
        }

        batch_.emplace_back(std::move(ovl_compact));

        if (batch_size > 0 && batch_.size() >= batch_size) {
            break;
        }
    }

    if (batch_.size() == 0) {
        return false;
    }

    return true;
}

int64_t FindQuerySpan(const std::vector<raptor::OverlapCompactPtr>& overlaps, int64_t start) {
	int64_t num_overlaps = static_cast<int64_t>(overlaps.size());
	int64_t span = 0;
	if (start < 0 || start >= num_overlaps) {
		return span;
	}
	for (int64_t i = start; i < num_overlaps; ++i) {
		if (overlaps[i]->a_id() != overlaps[start]->a_id()) {
			break;
		}
		++span;
	}
	return span;
}

}
