/*
 * output_formatter.cc
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#include <writer/output_formatter.h>
#include <sstream>

namespace raptor {

double OutputFormatter::CalcIdentity_(
            int32_t q_start, int32_t q_end,
            int32_t t_start, int32_t t_end,
            int32_t cov_bases_q, int32_t cov_bases_t,
            int32_t edit_dist, int32_t match_bases) {
    double qspan = q_end - q_start;
    double tspan = t_end - t_start;
    double identity_q = (qspan != 0.0) ? (static_cast<double>(cov_bases_q)) / (qspan) : 0.0;
    double identity_t = (tspan != 0.0) ? (static_cast<double>(cov_bases_t)) / (tspan) : 0.0;
    double identity_min = std::min(identity_q, identity_t);

    if (edit_dist >= 0) {
        double match_bases_q = (match_bases > 0) ? match_bases : std::max(0.0, qspan  - edit_dist);
        double match_bases_t = (match_bases > 0) ? match_bases : std::max(0.0, tspan  - edit_dist);
        identity_q = (qspan != 0) ? (match_bases_q / qspan) : 0.0;
        identity_t = (tspan != 0) ? (match_bases_t / tspan) : 0.0;
        identity_min = std::min(identity_q, identity_t);
    }

    return identity_min;
}

std::string OutputFormatter::TimingMapToString(const std::unordered_map<std::string, double>& timings) {
    std::ostringstream oss;
    oss << std::fixed;
    bool is_first = false;
    for (auto& it: timings) {
        if (is_first == false) {
            oss << "/";
        }
        oss << it.first << "=" << it.second << "us";
    }
    return oss.str();
}

std::string OutputFormatter::UnmappedSAM(const mindex::SequencePtr& qseq, bool write_custom_tags) {
    std::stringstream ss;

    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    uint32_t flag = 4;
    std::string seq = qseq->GetSequenceAsString();
    std::string qual = (qseq->qual().size() > 0) ? qseq->GetQualityAsString() : std::string("*");

    ss  << q_name << "\t"
        << flag << "\t"
        << "*" << "\t"
        << 0 << "\t"
        << 255 << "\t"
        << "*" << "\t"
        << "*" << "\t"
        << "0" << "\t"
        << "0" << "\t"
        << seq << "\t"
        << qual;

    // if (write_custom_tags) {
        // Extra tags provided in the alignment.
        for (const auto& vals: qseq->tags()) {
            ss << "\t" << vals.FormatAsSAM();
        }
    // }

    ss << "\n";

    return ss.str();
}

std::string OutputFormatter::ToSAM(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping,
                            bool write_custom_tags,
                            const std::string& timings) {
    std::stringstream ss;

    // int32_t q_id = mapping->QueryID();
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    int32_t score = mapping->Score();
    int32_t mapq = mapping->MappingQuality();

    int64_t path_id = mapping->PathId();
    int64_t num_paths = mapping->PathsNum();
    int64_t segment_in_path = mapping->SegmentId();
    int64_t num_segments_in_path = mapping->SegmentsNum();

    uint32_t flag = (t_is_rev) ? 16 : 0;
    if (mapping->IsSupplementary()) {
        flag |= 2048;
    }
    if (mapping->IsSecondary()) {
        flag |= 256;
    }

    std::string seq = qseq->GetSequenceAsString();
    std::string qual = (qseq->qual().size() > 0) ? qseq->GetQualityAsString() : std::string("*");

    std::string cigar = CigarToString(mapping->Cigar(), false);
    if (cigar.size() == 0) {
        cigar = "*";
    }

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;

        seq = ReverseComplement(seq);
        std::reverse(qual.begin(), qual.end());
    }

    ss  << q_name << "\t"
        << flag << "\t"
        << t_name << "\t"
        << t_start + 1 << "\t"
        << mapq << "\t"
        << cigar << "\t"
        << "*" << "\t"
        << "0" << "\t"
        << "0" << "\t"
        << seq << "\t"
        << qual;

    // Extra tags provided in the alignment.
    for (const auto& vals: qseq->tags()) {
        ss << "\t" << vals.FormatAsSAM();
    }

    if (write_custom_tags) {
        // Specific tags to Raptor.
        ss << "\t"
            << "NM:i:" << edit_dist << "\t"
            << "AS:i:" << score << "\t"
            << "QS:i:" << q_start << "\t"
            << "QE:i:" << q_end << "\t"
            << "QL:i:" << q_len << "\t"
            << "TS:i:" << t_start << "\t"
            << "TE:i:" << t_end << "\t"
            << "TL:i:" << t_len << "\t"
            << "pi:i:" << path_id << "\t"
            << "pj:i:" << segment_in_path << "\t"
            << "pn:i:" << num_segments_in_path << "\t"
            << "ps:i:" << ((num_segments_in_path > 1) ? 1 : 0);

#ifdef RAPTOR_DEBUG_TIMINGS
    ss << "\t" << "tt:Z:" << timings;
#endif
    }

    ss << "\n";

    return ss.str();
}

std::string OutputFormatter::ToPAF(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping,
                            bool write_custom_tags,
                            const std::string& timings) {
    std::stringstream ss;

    // int32_t q_id = mapping->QueryID();
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    int32_t score = mapping->Score();
    int32_t cov_bases_q = mapping->CoveredBasesQuery();
    int32_t num_seeds = mapping->NumSeeds();
    int32_t mapq = mapping->MappingQuality();

    int64_t path_id = mapping->PathId();
    int64_t num_paths = mapping->PathsNum();
    int64_t segment_in_path = mapping->SegmentId();
    int64_t num_segments_in_path = mapping->SegmentsNum();

    // This is how a flag is calculated for SAM/BAM alignments. It's not a standard in PAF
    // but will  be output as an optional tag named "fg".
    uint32_t flag = (t_is_rev) ? 16 : 0;
    if (mapping->IsSupplementary()) {
        flag |= 2048;
    }
    if (mapping->IsSecondary()) {
        flag |= 256;
    }

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

    std::string cigar = CigarToString(mapping->Cigar(), true);
    if (cigar.size() == 0) {
        cigar = "*";
    }

    ss <<   q_name << "\t" <<
            q_len << "\t" <<
            q_start << "\t" <<
            q_end << "\t" <<
            ((t_is_rev) ? "-" : "+") << "\t" <<
            t_name << "\t" <<
            t_len << "\t" <<
            t_start << "\t" <<
            t_end << "\t" <<
            cov_bases_q << "\t" <<
            (t_end - t_start) << "\t" <<
            mapq;

    if (write_custom_tags) {
        ss << "\t"
            // << "pr:i:" << mapping->GetRegionPriority() << "\t"
            // << "su:i:" << (mapping->GetRegionIsSupplementary() ? 1 : 0) << "\t"
            << "cm:i:" << num_seeds << "\t"
            << "NM:i:" << edit_dist << "\t"
            << "AS:i:" << score << "\t"
            << "fg:i:" << flag << "\t"
            << "pi:i:" << path_id << "\t"
            << "pj:i:" << segment_in_path << "\t"
            << "pn:i:" << num_segments_in_path << "\t"
            << "ps:i:" << ((num_segments_in_path > 1) ? 1 : 0) << "\t"
            << "cg:Z:" << cigar;

        #ifdef RAPTOR_DEBUG_TIMINGS
            ss << "\t" << "tt:Z:" << timings;
        #endif
    }

    // Extra tags provided in the alignment.
    for (const auto& vals: qseq->tags()) {
        ss << "\t" << vals.FormatAsSAM();
    }

    ss << "\n";

    return ss.str();
}

std::string OutputFormatter::ToMHAP(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                              const std::shared_ptr<raptor::RegionBase> mapping) {
    std::stringstream ss;

    int32_t q_id = mapping->QueryID();
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    // std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    // std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    // int32_t score = mapping->Score();
    int32_t cov_bases_q = mapping->CoveredBasesQuery();
    int32_t cov_bases_t = mapping->CoveredBasesTarget();
    int32_t match_bases = mapping->MatchBases();
    int32_t num_seeds = mapping->NumSeeds();

    double identity = CalcIdentity_(
                    q_start, q_end, t_start, t_end,
                    cov_bases_q, cov_bases_t, edit_dist, match_bases);

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

    ss  << q_id << " "
        << t_id << " "
        << identity << " "
        << num_seeds << " "
        << 0 << " "
        << q_start << " "
        << q_end << " "
        << q_len << " "
        << ((t_is_rev) ? "1" : "0") << " "
        << t_start << " "
        << t_end << " "
        << t_len;

    ss << std::endl;

    return ss.str();
}

std::string OutputFormatter::ToM4(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping) {
    std::stringstream ss;

    // int32_t q_id = mapping->QueryID();
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    int32_t score = mapping->Score();
    int32_t cov_bases_q = mapping->CoveredBasesQuery();
    int32_t cov_bases_t = mapping->CoveredBasesTarget();
    int32_t match_bases = mapping->MatchBases();
    // int32_t num_seeds = mapping->NumSeeds();

    double identity = CalcIdentity_(
                    q_start, q_end, t_start, t_end,
                    cov_bases_q, cov_bases_t, edit_dist, match_bases);

    if (t_is_rev) {
        // Output in the fwd strand always.
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

    ss  << q_name << " "
        << t_name << " "
        << -score << " "
        << 100.0 * identity << " "
        << 0 << " "
        << q_start << " "
        << q_end << " "
        << q_len << " "
        << ((t_is_rev) ? "1" : "0") << " "
        << t_start << " "
        << t_end << " "
        << t_len;
        // << " "
        // << mapq;

    ss << std::endl;

    return ss.str();
}

std::string OutputFormatter::ToGFA2Edge(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping,
                            bool write_custom_tags,
                            const std::string& timings) {
    std::stringstream ss;

    int32_t q_id = mapping->QueryID();
    int32_t q_start = mapping->QueryStart();
    int32_t q_end = mapping->QueryEnd();
    int32_t q_len = qseq->data().size();
    int32_t t_id = mapping->TargetID();
    int32_t t_start = mapping->TargetStart();
    int32_t t_end = mapping->TargetEnd();
    int32_t t_len = mapping->TargetLen();
    bool t_is_rev = mapping->TargetRev();
    std::string q_name = raptor::TrimToFirstSpace(qseq->header());
    std::string t_name = raptor::TrimToFirstSpace(index->header(t_id));
    int32_t edit_dist = mapping->EditDistance();
    int32_t score = mapping->Score();
    int32_t cov_bases_q = mapping->CoveredBasesQuery();
    int32_t num_seeds = mapping->NumSeeds();

    int64_t path_id = mapping->PathId();
    int64_t num_paths = mapping->PathsNum();
    int64_t segment_in_path = mapping->SegmentId();
    int64_t num_segments_in_path = mapping->SegmentsNum();

    // This is how a flag is calculated for SAM/BAM alignments. It's not a standard in PAF
    // but will  be output as an optional tag named "fg".
    uint32_t flag = (t_is_rev) ? 16 : 0;
    if (mapping->IsSupplementary()) {
        flag |= 2048;
    }
    if (mapping->IsSecondary()) {
        flag |= 256;
    }

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
    }

    std::string cigar = CigarToString(mapping->Cigar(), true);
    if (cigar.size() == 0) {
        cigar = "*";
    }

    std::ostringstream edge_name;
    edge_name << "aln-" << q_id << "-" << t_id << ":" << path_id << "-" << segment_in_path;

    std::string q_end_symbol = (q_end == q_len) ? "$" : "";
    std::string t_end_symbol = (t_end == t_len) ? "$" : "";

    ss << "E" << "\t"
        << edge_name.str() << "\t"
        << q_name << "+" << "\t"
        << t_name << (t_is_rev ? "-" : "+") << "\t"
        << q_start << "\t"
        << q_end << q_end_symbol << "\t"
        << t_start << "\t"
        << t_end << t_end_symbol << "\t"
        << cigar;

    if (write_custom_tags) {
        ss << "\t"
            << "cm:i:" << num_seeds << "\t"
            << "NM:i:" << edit_dist << "\t"
            << "AS:i:" << score << "\t"
            << "fg:i:" << flag << "\t"
            << "pi:i:" << path_id << "\t"
            << "pj:i:" << segment_in_path << "\t"
            << "pn:i:" << num_segments_in_path << "\t"
            << "ps:i:" << ((num_segments_in_path > 1) ? 1 : 0)
            << "\t" << "cg:Z:" << cigar;

        #ifdef RAPTOR_DEBUG_TIMINGS
            ss << "\t" << "tt:Z:" << timings;
        #endif
    }

    ss << std::endl;

    return ss.str();
}



}
