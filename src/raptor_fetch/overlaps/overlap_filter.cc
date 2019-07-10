#include <raptor_fetch/overlaps/overlap_filter.h>
#include <algorithm>
#include <deque>

namespace raptor {

int64_t FindQuerySpan(const std::vector<raptor::OverlapCompactPtr>& overlaps, int64_t start) {
	int64_t num_overlaps = static_cast<int64_t>(overlaps.size());
	int64_t span = 0;
	if (start < 0 || start >= num_overlaps) {
		return span;
	}
	for (int64_t i = start; i < num_overlaps; ++i) {
		++span;
		if (overlaps[i]->a_id() != overlaps[start]->a_id()) {
			break;
		}
	}
	return span;
}

bool FindClippingRegions(raptor::OverlapFilePtr& ovl_file, size_t batch_start, size_t batch_end,
                         int32_t fwd_dist,
						 int32_t min_user_cov_threshold, bool clip_ends_on_coverage,
						 std::vector<QueryRegion>& ret_regions) {

	ret_regions.clear();

	if (ovl_file->empty()) {
		return false;
	}

	// std::cerr << "ovl_file->batch()[0]->a_id() = " << ovl_file->batch()[0]->a_id() << "\n";

	int32_t a_id = ovl_file->batch()[batch_start]->a_id();
	auto& a_data = ovl_file->GetQueryDataById(a_id);
	if (a_data == nullptr) {
		std::cerr << "Warning: a_data == nullptr in FindClippingRegions!\n";
		return false;
	}
	int32_t a_len = a_data->len;

	std::vector<CoverageEvent> events;
	for (size_t i = batch_start; i < batch_end; ++i) {
        if (ovl_file->batch()[i]->IsFiltered()) {
            continue;
        }
		events.emplace_back(CoverageEvent(ovl_file->batch()[i]->a_start(), true, i, 0, 0));
		events.emplace_back(CoverageEvent(ovl_file->batch()[i]->a_end(), false, i, 0, 0));
	}
	std::sort(events.begin(), events.end(), [](const CoverageEvent& a, const CoverageEvent& b) { return (a.pos < b.pos || (a.pos == b.pos && a.is_start < b.is_start)); });

	std::vector<int32_t> covs(a_len + 1, 0);
	std::vector<int32_t> fwd_dists(a_len + 1, 0);

	// Find the coverage at each of the events.
	int32_t curr_cov = 0;
	for (size_t i = 0; i < events.size(); ++i) {
		if (events[i].is_start) {
			++curr_cov;
		} else {
			--curr_cov;
		}
		covs[events[i].pos] = curr_cov;
	}
	int32_t prev_cov = 0;
	double avg_cov = 0.0;
	for (size_t i = 0; i < events.size(); ++i) {
		events[i].cov = covs[events[i].pos];
		if (i > 0) {
			avg_cov += prev_cov * (events[i].pos - events[i-1].pos);
		}
		prev_cov = events[i].cov;
	}
	if (events.size() > 0) {
		avg_cov += (a_len - events.back().pos) * prev_cov;
	}
	avg_cov = (a_len > 0) ? (avg_cov / static_cast<double>(a_len)) : 0.0;

	// Find out how many closings there are in a window in front of each event.
	std::deque<size_t> window;
	int32_t num_events = static_cast<int32_t>(events.size());
	for (int32_t i = (num_events - 1); i >= 0; --i) {
		if (events[i].is_start == false) {
			window.push_back(i);
		}
		while (window.size() > 0 && (events[window.front()].pos - events[i].pos) > fwd_dist) {
			window.pop_front();
		}
		fwd_dists[events[i].pos] = window.size();
	}
	for (size_t i = 0; i < events.size(); ++i) {
		events[i].fwd_closing = fwd_dists[events[i].pos];
	}

	// std::cerr << "avg_cov = " << avg_cov << ", num_overlaps = " << ovl_file->batch().size() << ", a_len = " << a_len << "\n";
	// for (size_t i = 0; i < events.size(); ++i) {
	// 	std::cerr << "  [" << i << "] pos = " << events[i].pos << ", type = " << ((events[i].is_start) ? "start" : "end") <<
	// 		", source_id = " << events[i].source_id << ", cov = " << events[i].cov << ", fwd_closing = " << events[i].fwd_closing << "\n";
	// }
	// std::cerr << "\n";

    /*
    One good region:
    <------------------->
    |       GOOD        |

    Large insertion:
    <------>---<---------->
    | GOOD |BAD|   GOOD   |

    Chimera:
    <------><---------->
    | GOOD ||   GOOD   |

    Coverage drops at the end. This should be extended all the way.
    -<---------->-
    |   GOOD   |

    If clip_ends_on_coverage == false, we only want to clip on chimera, so low coverage on both ends is still good and should be:
    -<---------->-
    |    GOOD    |
    */

   	// Find the start and stop of every region.
	// const int32_t min_cov_threshold = std::max(std::max(static_cast<int32_t>(1), min_user_cov_threshold), static_cast<int32_t>(avg_cov / 2));
	const int32_t min_cov_threshold = (min_user_cov_threshold >= 0) ? min_user_cov_threshold : std::max(static_cast<int32_t>(1), static_cast<int32_t>(avg_cov / 2));
	// const int32_t max_num_closings = static_cast<int32_t>(avg_cov / 2);
	std::vector<CoverageEvent> region_events;
	bool curr_state = (!clip_ends_on_coverage);
	bool next_state = curr_state;
	int32_t prev_event_cov = 0;
	for (int32_t i = 0; i < num_events; ++i) {
		const auto& event = events[i];
		if (curr_state == false) {
			next_state = false;
			if (event.cov >= min_cov_threshold) {
				if (region_events.empty() && event.pos > 0) {
					region_events.emplace_back(CoverageEvent(0, curr_state, -1, 0, 0));
				}
				region_events.emplace_back(CoverageEvent(event.pos, true, event.source_id, event.cov, event.fwd_closing));
				next_state = true;
			}
		} else {
			next_state = true;
			const int32_t max_num_closings = prev_event_cov / 2;
			if (clip_ends_on_coverage && event.cov < min_cov_threshold) {
				if (region_events.empty() && event.pos > 0) {
					region_events.emplace_back(CoverageEvent(0, curr_state, -1, 0, 0));
				}
				region_events.emplace_back(CoverageEvent(event.pos, false, event.source_id, event.cov, event.fwd_closing));
				next_state = false;
			} else if (event.fwd_closing > max_num_closings) {
				if (region_events.empty() && event.pos > 0) {
					region_events.emplace_back(CoverageEvent(0, curr_state, -1, 0, 0));
				}
				region_events.emplace_back(CoverageEvent(event.pos, false, event.source_id, event.cov, event.fwd_closing));
				next_state = false;
			}
		}
		// Skip all events at the same position.
		if (curr_state != next_state) {
			int32_t new_i = i + 1;
			while(new_i < num_events) {
				if (events[new_i].pos != events[i].pos) {
					break;
				}
				++new_i;
			}
			i = new_i - 1;
		}
		curr_state = next_state;
		prev_event_cov = event.cov;
	}

	if (clip_ends_on_coverage) {
		if (curr_state == true) {
            // If there was a coverage drop such as happens near ends
            // of the queries, then artificially extend the good region to the end,
            // so that the last few bases get picked up.
            // In case of chimeras, there should be a jump from GOOD->BAD->GOOD, so in
            // that case we wouldn't have a BAD state in any case other than the lack of
            // coverage at the back of the query.
			region_events.emplace_back(CoverageEvent(a_len, false, -1, 0, 0));
		}
	} else {
		if (region_events.size() > 0 && curr_state == false) {
			region_events.back().pos = a_len;
		} else {
			// Close the last good interval if required.
			region_events.emplace_back(CoverageEvent(a_len, false, -1, 0, 0));
		}
	}

	ret_regions.clear();
	for (size_t i = 1; i < region_events.size(); ++i) {
		const auto& event = region_events[i];
		if (event.is_start) {
			continue;
		}
		ret_regions.emplace_back(QueryRegion(a_id, region_events[i-1].pos, region_events[i].pos, a_len));
	}


	// std::cerr << "QueryRegion events:\n";
	// for (size_t i = 0; i < region_events.size(); ++i) {
	// 	const auto& event = region_events[i];
	// 	std::cerr << "  [" << i << "] pos = " << event.pos << ", type = " << ((event.is_start) ? "start" : "end") <<
	// 		", source_id = " << event.source_id << ", cov = " << event.cov << ", fwd_closing = " << event.fwd_closing << ", a_len = " << a_len << "\n";
	// }
	// std::cerr << "QueryRegions:\n";
	// for (size_t i = 0; i < ret_regions.size(); ++i) {
	// 	std::cerr << "  [" << i << "] a_id = " << ret_regions[i].id << ", start = " << ret_regions[i].start << ", end = " << ret_regions[i].end << ", len = " << ret_regions[i].len << "\n";
	// }
	// std::cerr << "\n";

	// for (size_t i = 0; i < ret_regions.size(); ++i) {
	// 	std::cout << "a_id = " << ret_regions[i].id << ", start = " << ret_regions[i].start << ", end = " << ret_regions[i].end << ", len = " << ret_regions[i].len << "\n";
	// }

	return true;
}

int32_t FilterRegions(const std::vector<QueryRegion>& in_regions, int32_t min_span, int32_t bestn, std::vector<QueryRegion>& out_regions) {
	/*
	Takes a set of input regions, skips the ones of span < min_span, and returns the bestn of those regions (best being defined by
	span length).
	Returns the number of total number of regions above min_span, even if the bestn is lower. This is useful because:
		- 0 - No valid regions.
		- 1 - One good valid region. Non-chimeric read.
		- >= 2 - 2 or more valid regions. Most likely a chimeric read.
	*/

	int32_t num_valid = 0;

	out_regions.clear();

	std::vector<int32_t> sorted_reg_ids(in_regions.size(), 0);
	for (size_t i = 0; i < in_regions.size(); ++i) {
		sorted_reg_ids[i] = i;
	}
	std::sort(sorted_reg_ids.begin(), sorted_reg_ids.end(),
				[&](int32_t a, int32_t b) {
					return (in_regions[a].end - in_regions[a].start) > (in_regions[b].end - in_regions[b].start);
				} );

	for (size_t i = 0; i < sorted_reg_ids.size(); ++i) {
		int32_t reg_id = sorted_reg_ids[i];
		const auto& region = in_regions[reg_id];
		int32_t span = region.end - region.start;
		if (span < min_span) {
			continue;
		}
		if (out_regions.size() < static_cast<size_t>(bestn)) {
			out_regions.emplace_back(region);
		}
		++num_valid;
	}

	return num_valid;
}

raptor::OverlapCompactPtr GetDovetailOverlap(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 int32_t max_allowed_clip, int32_t min_edge_len) {

	raptor::OverlapCompactPtr ret = raptor::createOverlapCompact(ovl);

	int32_t left_a = ovl->a_start();
	int32_t right_a = ovl->a_end();
	int32_t left_b = ovl->b_start();
	int32_t right_b = ovl->b_end();

	if (ovl->b_rev()) {
		// Calculation is easier in the strand of the overlap.
		std::swap(left_b, right_b);
		left_b = qdb->len - left_b;
		right_b = qdb->len - right_b;
	}

	int32_t offset_left = 0;
	int32_t offset_right = 0;
	int32_t type = 0;

	if (left_a < max_allowed_clip && (qdb->len - right_b) < max_allowed_clip) {
		/*
			Dovetail 5'.
                             left_a            right_a
		                     >|//|<          >|///////|<
			A:	              o--=============-------->
			B:     o-------------============->
                  >|/////////////|<        >|/|<
                        left_b            right_b
		*/
		offset_left = left_a;
		offset_right = qdb->len - right_b;
		type = 5;
	}
	else if ((qda->len - right_a) < max_allowed_clip && (left_b) < max_allowed_clip) {
		/*
			Dovetail 5'.
                        left_a            right_a
                  >|/////////////|<        >|/|<
			A:     o-------------============->
			B:	              o--=============-------->
		                     >|//|<          >|///////|<
                             left_b            right_b
		*/
		offset_left = left_b;
		offset_right = qda->len - right_a;
		type = 3;
	}

	if (type == 3 || type == 5) {
		ret->a_start(std::max(0, left_a - offset_left));
		ret->a_end(std::min(qda->len, right_a + offset_right));

		int32_t b_start = std::max(0, left_b - offset_left);
		int32_t b_end = std::min(qdb->len, right_b + offset_right);

		if (ret->b_rev()) {
			// Make the coords in the fwd orientation again.
			std::swap(b_start, b_end);
			b_start = qdb->len - b_start;
			b_end = qdb->len - b_end;
		}

		ret->b_start(b_start);
		ret->b_end(b_end);

		if (type == 5) {
			ret->SetType5Prime(true);
		} else if (type == 3) {
			ret->SetType3Prime(true);
		}
	} else {
		ret = nullptr;
	}

	if (left_a == 0 && right_a == 0) {
		ret = nullptr;
	}
	if (left_b == 0 && right_b == 0) {
		ret = nullptr;
	}

	if (ret == nullptr) {
		return ret;
	}

	int32_t new_left_hang_a = ret->a_start();
	int32_t new_right_hang_a = qda->len - ret->a_end();
	int32_t new_left_hang_b = ret->b_start();
	int32_t new_right_hang_b = qdb->len - ret->b_end();

    if (ret->b_rev()) {
        std::swap(new_left_hang_b, new_right_hang_b);
    }

    if (type == 3 && (new_left_hang_a < min_edge_len || new_right_hang_b < min_edge_len)) {
		ret = nullptr;
	} else if (type == 5 && (new_right_hang_a < min_edge_len || new_left_hang_b < min_edge_len)) {
		ret = nullptr;
	}

	return ret;
}

void FindQueryClips(int32_t fwd_dist, int32_t min_cov, int32_t min_len, raptor::OverlapFilePtr& overlap_file) {
	if (overlap_file->batch().empty()) {
		return;
	}

	auto& batch = overlap_file->batch();

	size_t start_i = 0;
	int32_t a_id = batch[start_i]->a_id();

	for (size_t end_i = 1; end_i <= batch.size(); ++end_i) {
		int32_t new_a_id = (end_i < batch.size()) ? batch[end_i]->a_id() : -1;

		if (new_a_id != a_id) {
			auto& qd = overlap_file->GetQueryDataById(a_id);

			if (qd == nullptr) {
				continue;
			}

			qd->info += "n_ovl=" + std::to_string((end_i - start_i)) + ";";

			// Find clipping regions. This will detect possible chimera too.
			std::vector<raptor::QueryRegion> all_regions;
			bool rv_find_clip = raptor::FindClippingRegions(overlap_file, start_i, end_i, fwd_dist, min_cov, true, all_regions);
			if (rv_find_clip == false) {
				continue;
			}

			// Filter the very short regions, such as the front and back of each sequence.
			std::vector<raptor::QueryRegion> regions;
			int32_t rv_filter = FilterRegions(all_regions, min_len, 1, regions);
			// bool is_chimera = (rv_filter > 1);	// Just for show, how to detect chimera.
			qd->info += "rv_filter=" + std::to_string(rv_filter) + ";regions=";
			for (size_t region_id = 0; region_id < regions.size(); ++region_id) {
				if (region_id > 0) {
					qd->info += ",";
				}
				qd->info += std::to_string(regions[region_id].start) + "-" + std::to_string(regions[region_id].end);
			}
			qd->info += ";";

			if (rv_filter > 0) {
				qd->clip_start = regions[0].start;
				qd->clip_end = regions[0].end;
			}
			if (rv_filter == 0) {
				qd->SetFlagFiltered(true);
			} else if (rv_filter > 1) {
				qd->SetFlagChimeric(true);
			}

			start_i = end_i;
		}

		a_id = new_a_id;
	}
}

void ApplyClips(raptor::OverlapFilePtr& overlap_file) {
	if (overlap_file->batch().empty()) {
		return;
	}

	auto& batch = overlap_file->batch();

	for (size_t oid = 0; oid < batch.size(); ++oid) {

		auto& ovl = batch[oid];
		auto& qd_a = overlap_file->GetQueryDataById(ovl->a_id());
		auto& qd_b = overlap_file->GetQueryDataById(ovl->b_id());
		if (qd_a == nullptr || qd_b == nullptr) {
			continue;
		}
		// First update, then checl;
		if (qd_a->IsFiltered()) {
			ovl->SetFlagFilteredA(true);
		}
		if (qd_b->IsFiltered()) {
			ovl->SetFlagFilteredB(true);
		}
		if (ovl->IsFiltered()) {
			continue;
		}

		int32_t new_a_start = ovl->a_start();
		int32_t new_a_end = ovl->a_end();
		int32_t new_b_start = ovl->b_start();
		int32_t new_b_end = ovl->b_end();

		if (ovl->b_rev()) {
			// Target.
			new_a_start = (ovl->b_end() < qd_b->clip_end) ? ovl->a_start() : (ovl->a_start() + (ovl->b_end() - qd_b->clip_end));
			new_a_end = (ovl->b_start() > qd_b->clip_start) ? ovl->a_end() : (ovl->a_end() - (qd_b->clip_start - ovl->b_start()));
			// Query
			new_b_start = (ovl->a_end() < qd_a->clip_end) ? ovl->b_start() : (ovl->b_start() + (ovl->a_end() - qd_a->clip_end));
			new_b_end = (ovl->a_start() > qd_a->clip_start) ? ovl->b_end() : (ovl->b_end() - (qd_a->clip_start - ovl->a_start()));

		} else {
			// Target.
			new_a_start = (ovl->b_start() > qd_b->clip_start) ? ovl->a_start() : (ovl->a_start() + (qd_b->clip_start - ovl->b_start()));
			new_a_end = (ovl->b_end() < qd_b->clip_end) ? ovl->a_end() : (ovl->a_end() - (ovl->b_end() - qd_b->clip_end));
			// Query.
			new_b_start = (ovl->a_start() > qd_a->clip_start) ? ovl->b_start() : (ovl->b_start() + (qd_a->clip_start - ovl->a_start()));
			new_b_end = (ovl->a_end() < qd_a->clip_end) ? ovl->b_end() : (ovl->b_end() - (ovl->a_end() - qd_a->clip_end));
		}

		new_a_start = (new_a_start > qd_a->clip_start) ? new_a_start : qd_a->clip_start;
		new_a_end = (new_a_end < qd_a->clip_end) ? new_a_end : qd_a->clip_end;
		new_b_start = (new_b_start > qd_b->clip_start) ? new_b_start : qd_b->clip_start;
		new_b_end = (new_b_end < qd_b->clip_end) ? new_b_end : qd_b->clip_end;

		ovl->a_start(new_a_start);
		ovl->a_end(new_a_end);
		ovl->b_start(new_b_start);
		ovl->b_end(new_b_end);
	}
}

std::string FormatOverlapToCage11(const raptor::OverlapCompactPtr& ovl, const raptor::QueryDataPtr& qda, const raptor::QueryDataPtr& qdb) {
	std::ostringstream oss;
	oss << qda->name << " " << qdb->name << " " << ovl->score() << " " << ovl->identity() << " "
		<< (ovl->a_rev() ? 1 : 0) << " " << ovl->a_start() << " " << ovl->a_end() << " " << qda->len << " "
		<< (ovl->b_rev() ? 1 : 0) << " " << ovl->b_start() << " " << ovl->b_end() << " " << qdb->len << " "
		<< qda->clip_start << " " << qda->clip_end << " "
		<< qdb->clip_start << " " << qdb->clip_end << " "
		<< ovl->flag();
	return oss.str();
}

std::string FormatOverlapForFalcon(const raptor::OverlapCompactPtr& ovl, const raptor::QueryDataPtr& qda, const raptor::QueryDataPtr& qdb) {
	std::ostringstream oss;
	oss << qda->name << " " << qdb->name << " " << ovl->score() << " " << ovl->identity() << " "
		<< (ovl->a_rev() ? 1 : 0) << " " << ovl->a_start() << " " << ovl->a_end() << " " << qda->len << " "
		<< (ovl->b_rev() ? 1 : 0) << " " << ovl->b_start() << " " << ovl->b_end() << " " << qdb->len;
	if (ovl->GetType5Prime()) {
		oss << " 5";
	}
	if (ovl->GetType3Prime()) {
		oss << " 3";
	}
	return oss.str();
}

// 000000001 000010111 1948 0.149452 0 4651 17692 18319 1 48 13403 13505 299 18244 267 13330 0 contained_A



// raptor::OverlapCompactPtr ClipOverlapCoords(const raptor::OverlapCompactPtr& ovl, const raptor::QueryDataPtr& qda, const raptor::QueryDataPtr& qdb) {
// 	raptor::OverlapCompactPtr ret = raptor::createOverlapCompact(ovl);

// 	ret->SetFlagClippedToRegion(true);
// 	ret->a_start(ret->a_start() - qda->clip_start);
// 	ret->a_end(ret->a_end() - qda->clip_start);
// 	ret->b_start(ret->b_start() - qdb->clip_start);
// 	ret->b_end(ret->b_end() - qdb->clip_start);

// 	return ret;
// }

// bool CheckOverlapContained(const raptor::OverlapCompactPtr& ovl) {
// }

raptor::OverlapPtr CalculateClippedOverlap(const raptor::OverlapCompactPtr& ovl,
								 			const raptor::QueryDataPtr& qda,
								 			const raptor::QueryDataPtr& qdb) {
	// Calculate the clipped coordinates.
	int32_t a_start = ovl->a_start() - qda->clip_start;
	int32_t a_end = ovl->a_end() - qda->clip_start;
	int32_t a_len = qda->clip_end - qda->clip_start;
	int32_t a_span = a_end - a_start;
	int32_t b_start = ovl->b_start() - qdb->clip_start;
	int32_t b_end = ovl->b_end() - qdb->clip_start;
	int32_t b_len = qdb->clip_end - qdb->clip_start;
	int32_t b_span = b_end - b_start;
	// // Get the coordinates in the strand of the overlap, to make things clearer.
	// if (ovl->b_rev()) {
	// 	std::swap(b_start, b_end);
	// 	b_start = b_len - b_start;
	// 	b_end = b_len - b_end;
	// }

	std::string a_name = qda->name + ":" + std::to_string(qda->clip_start) + "-" + std::to_string(qda->clip_end);
	std::string b_name = qdb->name + ":" + std::to_string(qdb->clip_start) + "-" + std::to_string(qdb->clip_end);
	raptor::OverlapPtr ret = raptor::createOverlap(
			a_name, b_name, ovl->score(), ovl->identity(),
			ovl->a_rev(), a_start, a_end, a_len,
			ovl->b_rev(), b_start, b_end, b_len
	);
	ret->a_id(ovl->a_id());
	ret->b_id(ovl->b_id());

	return ret;
}

OverlapType DetermineOverlapType(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 int32_t max_valid_aln_clip,							// Used for calling contained reads, should be smaller tan max_internal_hang, but not very high to not call SV regions falsely.
								 int32_t contained_min_dist_from_edge,						// For containment. If e.g. A is end-to-end contained within B, it still has to be this far from one of the edges of B, just in case.
								 int32_t internal_min_hang, float min_frac_cov		// Used for checking internal overlaps.
								 ) {

	raptor::OverlapPtr ovlc = CalculateClippedOverlap(ovl, qda, qdb);
	if (ovlc->b_rev()) {
		ovlc->FlipBCoords();
	}

	int32_t alhang = ovlc->a_start();
	int32_t arhang = ovlc->a_len() - ovlc->a_end();
	int32_t blhang = ovlc->b_start();
	int32_t brhang = ovlc->b_len() - ovlc->b_end();

	int32_t min_lhang = (alhang < blhang) ? alhang : blhang;
	int32_t max_lhang = (alhang > blhang) ? alhang : blhang;
	int32_t min_rhang = (arhang < brhang) ? arhang : brhang;
	int32_t max_rhang = (arhang > brhang) ? arhang : brhang;

	// Check if the overlap is an internal hit.
	if (min_lhang >= internal_min_hang || min_rhang >= internal_min_hang) {
		std::cerr << "-> Internal (1) because: min_lhang = " << min_lhang << ", min_rhang = " << min_rhang << ", internal_min_hang = " << internal_min_hang << "\n";
		return OverlapType::Internal;
	}
	// if ((alhang - blhang) >= internal_min_hang && ( arhang - brhang) >= internal_min_hang) {
	// 	return OverlapType::Internal;
	// }

	// if ((a_span + min_lhang + min_rhang) < a_len * min_frac_cov && (b_span + min_lhang + min_rhang) < b_len * min_frac_cov) {
	// 	// Both have to be true, otherwise it can be a containment.
	// 	std::cerr << "-> Internal (1) because: min_lhang = " << min_lhang << ", min_rhang = " << min_rhang << ", internal_min_hang = " << internal_min_hang << "\n";
	// 	return OverlapType::Internal;
	// }

	// If either side has an overhang larger than max_allowed_hang, we best not apply the
	// containment calling, as we might be in an SV region.
	if (min_lhang > max_valid_aln_clip || min_rhang > max_valid_aln_clip) {
		return OverlapType::Sketchy;
	}

	alhang -= min_lhang;
	arhang -= min_rhang;
	blhang -= min_lhang;
	brhang -= min_rhang;

	if (alhang == 0 && arhang == 0 && (blhang >= contained_min_dist_from_edge || brhang >= contained_min_dist_from_edge)) {
		return OverlapType::ContainedA;
	}
	if (blhang == 0 && brhang == 0 && (alhang >= contained_min_dist_from_edge || arhang >= contained_min_dist_from_edge)) {
		return OverlapType::ContainedB;
	}

	return OverlapType::Overlap;
}

OverlapType DetermineOverlapTypeNew(const raptor::OverlapCompactPtr& ovl,
								 const raptor::QueryDataPtr& qda,
								 const raptor::QueryDataPtr& qdb,
								 // If min_lhang or min_hang is > max_valid_aln_clip, it's called an internal overlap.
								 int32_t internal_max_valid_aln_clip,
								 // Used for checking internal overlaps.
								 float internal_min_frac_cov,
								 // Perform coordinate offsetting only if the distance is <= offset_max_dist. Offsetting is performed before calling contained and dovetail overlaps.
								 int32_t offset_max_dist,
                                 // Similar to internal_max_valid_aln_clip, but should be smaller. If unaligned portion on the shorter end is <= internal_max_valid_aln_clip,
                                 // this value will be used to offset the coordinates and create fake dovetail overlaps.
                                 // The fake dovetail overlap will be used to call contained reads or dovetail overlaps.
								 int32_t dovetail_max_valid_aln_clip,
 								 // For containment. If e.g. A is end-to-end contained within B, it still has to be this far from one of the edges of B, just in case.
								 int32_t contained_min_dist_from_edge
								 ) {

	raptor::OverlapPtr ovlc = CalculateClippedOverlap(ovl, qda, qdb);
	if (ovl->b_rev()) {
		ovlc->FlipBCoords();
        // std::cerr << "ovlc->b_start = " << ovlc->b_start() << ", ovlc->b_end() = " << ovlc->b_end() << "\n";
	}

	int32_t alhang = ovlc->a_start();
	int32_t arhang = ovlc->a_len() - ovlc->a_end();
	int32_t blhang = ovlc->b_start();
	int32_t brhang = ovlc->b_len() - ovlc->b_end();

	int32_t min_lhang = (alhang < blhang) ? alhang : blhang;
	int32_t max_lhang = (alhang > blhang) ? alhang : blhang;
	int32_t min_rhang = (arhang < brhang) ? arhang : brhang;
	int32_t max_rhang = (arhang > brhang) ? arhang : brhang;
	int32_t min_ahang = std::min(alhang, arhang);
	int32_t min_bhang = std::min(blhang, brhang);

#ifdef DEBUG_OVERLAP_TYPE
	std::cerr << "alhang = " << alhang << ", arhang = " << arhang << ", blhang = " << blhang << ", brhang = " << brhang
			<< ", min_lhang = " << min_lhang << ", max_lhang = " << max_lhang << ", min_rhang = " << min_rhang << ", max_rhang = " << max_rhang << "\n";
#endif

	// If either side is longer than a specific threshold, call it internal.
	if (min_ahang > internal_max_valid_aln_clip && min_bhang > internal_max_valid_aln_clip) {
#ifdef DEBUG_OVERLAP_TYPE
		std::cerr << "-> Internal (1) because: minimum ahang is longer than allowed fixed length (min_ahang = " << min_ahang << ") > (internal_max_valid_aln_clip = " << internal_max_valid_aln_clip << ")\n";
		std::cerr << "-> or because: minimum bhang is longer than allowed fixed length (min_bhang = " << min_ahang << ") > (internal_max_valid_aln_clip = " << internal_max_valid_aln_clip << ")\n";
#endif
		return OverlapType::Internal;
	}

	// // Check if the overlap is an internal hit.
	// if (min_lhang > max_valid_aln_clip || min_rhang > max_valid_aln_clip) {
	// 	std::cerr << "-> Internal (5) because: min_lhang = " << min_lhang << ", min_rhang = " << min_rhang << ", max_valid_aln_clip = " << max_valid_aln_clip << "\n";
	// 	return OverlapType::Internal;
	// }

	// if (ovlc->a_span() < (ovlc->a_span() + min_lhang + min_rhang) * min_frac_cov)

	// if (min_lhang >=

	// Make an internal dovetail overlap.
    int32_t offset_left = std::min(alhang, blhang);
    int32_t offset_right = std::min(arhang, brhang);
    if (offset_left > internal_max_valid_aln_clip) {
        offset_left = 0;
    }
    if (offset_right > internal_max_valid_aln_clip) {
        offset_right = 0;
    }

    alhang -= offset_left;
    arhang -= offset_right;
    blhang -= offset_left;
    brhang -= offset_right;

    // int32_t type = 0;

	// if (alhang < dovetail_max_valid_aln_clip && brhang < dovetail_max_valid_aln_clip) {
	// 	offset_left = alhang;
	// 	offset_right = brhang;
    //     type = 5;
	// }
	// else if (arhang < dovetail_max_valid_aln_clip && blhang < dovetail_max_valid_aln_clip) {
	// 	offset_left = blhang;
	// 	offset_right = arhang;
	// 	type = 3;
	// }

	// alhang = std::max(0, alhang - offset_left);
	// arhang = std::max(0, arhang - offset_right);
	// blhang = std::max(0, blhang - offset_left);
	// brhang = std::max(0, brhang - offset_right);
#ifdef DEBUG_OVERLAP_TYPE
    std::cerr << "After updating hangs:\n";
	std::cerr << "alhang = " << alhang << ", arhang = " << arhang << ", blhang = " << blhang << ", brhang = " << brhang
			<< ", min_lhang = " << min_lhang << ", max_lhang = " << max_lhang << ", min_rhang = " << min_rhang << ", max_rhang = " << max_rhang << "\n";
#endif

	if (alhang <= dovetail_max_valid_aln_clip && arhang <= dovetail_max_valid_aln_clip && (blhang >= contained_min_dist_from_edge || brhang >= contained_min_dist_from_edge)) {
		return OverlapType::ContainedA;
	}
	if (blhang <= dovetail_max_valid_aln_clip && brhang <= dovetail_max_valid_aln_clip && (alhang >= contained_min_dist_from_edge || arhang >= contained_min_dist_from_edge)) {
		return OverlapType::ContainedB;
	}

    if (alhang <= dovetail_max_valid_aln_clip && brhang <= dovetail_max_valid_aln_clip && arhang > dovetail_max_valid_aln_clip && blhang > dovetail_max_valid_aln_clip) {
        return OverlapType::Dovetail5;
    }
    if (arhang <= dovetail_max_valid_aln_clip && blhang <= dovetail_max_valid_aln_clip && alhang > dovetail_max_valid_aln_clip && brhang > dovetail_max_valid_aln_clip) {
        return OverlapType::Dovetail3;
    }

	/*
	This section (with internal_min_frac_cov) needs to go below containment, otherwise overlaps like this:
        {"A B -2156 19.4866 0 0 11468 11468 1 3528 14426 24275", "contained_A"},
	would be called "interna" instead of "contained_A".
	*/
	if (min_ahang > ovlc->a_len() * internal_min_frac_cov) {
#ifdef DEBUG_OVERLAP_TYPE
		std::cerr << "-> Internal (2) because: minimum ahang is longer than allowed fraction of length. (min_ahang = " << min_ahang << ") > (ovlc->a_len() * internal_min_frac_cov = " << ovlc->a_len() * internal_min_frac_cov << ")\n";
#endif
		return OverlapType::Internal;
	}
	// if (min_bhang > internal_max_valid_aln_clip) {
	// 	std::cerr << "-> Internal (3) because: minimum ahang is longer than allowed fixed length (min_bhang = " << min_ahang << ") > (internal_max_valid_aln_clip = " << internal_max_valid_aln_clip << ")\n";
	// 	return OverlapType::Internal;
	// }
	if (min_bhang > ovlc->b_len() * internal_min_frac_cov) {
#ifdef DEBUG_OVERLAP_TYPE
		std::cerr << "-> Internal (4) because: minimum bhang is longer than allowed fraction of length. (min_bhang = " << min_bhang << ") > (ovlc->a_len() * internal_min_frac_cov = " << ovlc->a_len() * internal_min_frac_cov << ")\n";
#endif
		return OverlapType::Internal;
	}


    // if (type == 5) {
    //     return OverlapType::Dovetail5;
    // }
    // if (type == 3) {
    //     return OverlapType::Dovetail3;
    // }

    return OverlapType::Sketchy;

	// return OverlapType::Overlap;
}

}
