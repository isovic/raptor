#include <writer/raptor_results_writer_bam.h>
#include <writer/output_formatter.h>
#include <fstream>
#include <log/log_tools.h>
#include <pbbam/BamRecord.h>
#include <pbbam/Cigar.h>

namespace raptor {

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterBAM(const std::string& out_fn, const mindex::IndexPtr index) {
    return std::unique_ptr<raptor::RaptorResultsWriterBase>(new raptor::RaptorResultsWriterBAM(out_fn, index));
}

RaptorResultsWriterBAM::RaptorResultsWriterBAM(const std::string& out_fn,
                                        const mindex::IndexPtr index)
                                        : out_fn_(out_fn),
                                            bam_writer_(nullptr),
                                            index_(index) {

}

RaptorResultsWriterBAM::~RaptorResultsWriterBAM() {

}

std::unique_ptr<PacBio::BAM::BamWriter> RaptorResultsWriterBAM::OpenBAMWriter_(const std::string out_fn, PacBio::BAM::BamHeader& header) const {
    return std::make_unique<PacBio::BAM::BamWriter>(out_fn, header);
}

void RaptorResultsWriterBAM::WriteHeader(const mindex::HeaderGroupType header_groups) {
    std::ostringstream sam_header;
    sam_header << "@HD\tVN:1.5" << std::endl;
    for (const auto& it_field: header_groups) {
        sam_header << "@" << it_field.first;
        for (const auto& it_ids: it_field.second) {
            for (const auto& it_tags: it_ids.second) {
                sam_header << "\t" << it_tags.name << ":" << it_tags.val;
            }
        }
        sam_header << "\n";
    }
    for (size_t i = 0; i < index_->seqs()->size(); ++i) {
        std::string qname = TrimToFirstSpace(index_->header(i));
        sam_header << "@SQ\tSN:" << qname << "\tLN:" << index_->len(i) << std::endl;
    }

    PacBio::BAM::BamHeader header(sam_header.str());
    bam_writer_ = std::move(OpenBAMWriter_(out_fn_, header));

}

void RaptorResultsWriterBAM::WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<RaptorResults>& results, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (bam_writer_ == nullptr) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Output BAM file not opened. Ensure that WriteHeader is called first.");
    }

    for (auto& result: results) {
        WriteSingleResult(seqs, result, is_alignment_applied, write_custom_tags, one_hit_per_target);
    }
}

void RaptorResultsWriterBAM::WriteSingleResult(const mindex::SequenceFilePtr seqs, const RaptorResults& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (bam_writer_ == nullptr) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Output BAM file not opened. Ensure that WriteHeader is called first.");
    }

    int32_t mapq = 0;
    std::string timings_all;
    bool do_output = false;
    const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result.regions;

    // Collect the results to write.
    // Branching is because the output can either be aligned or unaligned.
    if (is_alignment_applied) {
        do_output = (result.aln_result != nullptr && result.aln_result->path_alignments().empty() == false);
        if (do_output) {
            mapq = result.aln_result->CalcMapq();
            std::string map_timings = OutputFormatter::TimingMapToString(result.mapping_result->timings());
            std::string graphmap_timings = OutputFormatter::TimingMapToString(result.graph_mapping_result->timings());
            timings_all = map_timings + "///" + graphmap_timings;
        }

    } else {
        do_output = (result.graph_mapping_result != nullptr && result.graph_mapping_result->paths().empty() == false);
        if (do_output) {
            mapq = result.graph_mapping_result->CalcMapq();
            std::string map_timings = OutputFormatter::TimingMapToString(result.mapping_result->timings());
            std::string graphmap_timings = OutputFormatter::TimingMapToString(result.graph_mapping_result->timings());
            timings_all = map_timings + "///" + graphmap_timings;
        }
    }

    // The writing code is generic.
    if (do_output && !regions_to_write.empty()) {
        for (size_t i = 0; i < regions_to_write.size(); i++) {
            auto aln = regions_to_write[i];
            bool is_secondary = (aln->PathId() > 0);
            bool is_supplementary = (aln->SegmentId() > 0);
            auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());

            auto record = ToBAM(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            bam_writer_->Write(record);

            // std::unique_ptr<PacBio::BAM::BamRecord> record = std::move(ToBAM(index_, qseq, aln, mapq, write_custom_tags, timings_all));
            // bam_writer_->Write(*record);
        }
    } else {
        // If the q_id_in_batch < 0, it means that the result was initialized, but never
        // updated.
        // This can happen in processing a certain range of reads, and all the other ones
        // outside this range will not be processed, but could have still been allocated.
        if (result.q_id_in_batch >= 0) {
            const auto& qseq = seqs->GetSeqByID(result.q_id_in_batch);

            auto record = ToUnmappedBAM(qseq);
            bam_writer_->Write(record);

        }
    }
}

PacBio::BAM::BamRecord RaptorResultsWriterBAM::ToUnmappedBAM(const mindex::SequencePtr& qseq) {
    PacBio::BAM::BamRecord record;

    if (qseq->apriori_bam()) {
        record = PacBio::BAM::BamRecord(*qseq->apriori_bam());

    } else {
        std::string q_name = raptor::TrimToFirstSpace(qseq->header());
        int32_t q_len = qseq->data().size();
        std::string seq = qseq->GetSequenceAsString();
        std::string qual = (qseq->qual().size() > 0) ? qseq->GetQualityAsString() : std::string(q_len, '!');

        PacBio::BAM::BamRecordImpl record_impl;
        record_impl.Name(q_name);
        record_impl.SetSequenceAndQualities(seq, qual);
        record = PacBio::BAM::BamRecord(record_impl);
    }

    return record;
}

PacBio::BAM::BamRecord RaptorResultsWriterBAM::ToBAM(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping,
                            int32_t mapq, bool write_custom_tags,
                            const std::string& timings) {

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

    // The mapq is generated from the outside because we need the number of alternative
    // alignment positions to generate it.

    std::string seq = qseq->GetSequenceAsString();
    std::string qual = (qseq->qual().size() > 0) ? qseq->GetQualityAsString() : std::string(q_len, '!');
    std::string cigar = CigarToString(mapping->Cigar(), false);

    if (t_is_rev) {
        std::swap(t_start, t_end);
        t_start = t_len - t_start;
        t_end = t_len - t_end;
        // seq = ReverseComplement(seq);
        // std::reverse(qual.begin(), qual.end());
    }

    PacBio::BAM::BamRecord record;

    if (qseq->apriori_bam()) {
        record = qseq->apriori_bam()->Mapped(
                static_cast<int32_t>(t_id),
                static_cast<int32_t>(t_start),
                (t_is_rev ? PacBio::BAM::Strand::REVERSE : PacBio::BAM::Strand::FORWARD),
                PacBio::BAM::Cigar(cigar),
                static_cast<uint8_t>(mapq)
                );

    } else {
        PacBio::BAM::BamRecordImpl record_impl;
        record_impl.Name(q_name);
        record_impl.SetSequenceAndQualities(seq, qual);
        record = PacBio::BAM::BamRecord(record_impl);
        record.Map(
                    static_cast<int32_t>(t_id),
                    static_cast<int32_t>(t_start),
                    (t_is_rev ? PacBio::BAM::Strand::REVERSE : PacBio::BAM::Strand::FORWARD),
                    PacBio::BAM::Cigar(cigar),
                    static_cast<uint8_t>(mapq)
                );
    }

    return record;




//     ss  << q_name << "\t"
//         << flag << "\t"
//         << t_name << "\t"
//         << t_start + 1 << "\t"
//         << mapq << "\t"
//         << cigar << "\t"
//         << "*" << "\t"
//         << "0" << "\t"
//         << "0" << "\t"
//         << seq << "\t"
//         << qual;

//     // Extra tags provided in the alignment.
//     for (const auto& vals: qseq->tags()) {
//         ss << "\t" << vals.FormatAsSAM();
//     }

//     if (write_custom_tags) {
//         // Specific tags to Raptor.
//         ss << "\t"
//             << "NM:i:" << edit_dist << "\t"
//             << "AS:i:" << score << "\t"
//             << "QS:i:" << q_start << "\t"
//             << "QE:i:" << q_end << "\t"
//             << "QL:i:" << q_len << "\t"
//             << "TS:i:" << t_start << "\t"
//             << "TE:i:" << t_end << "\t"
//             << "TL:i:" << t_len << "\t"
//             << "pi:i:" << path_id << "\t"
//             << "pj:i:" << segment_in_path << "\t"
//             << "pn:i:" << num_segments_in_path << "\t"
//             << "ps:i:" << ((num_segments_in_path > 1) ? 1 : 0);

// #ifdef RAPTOR_DEBUG_TIMINGS
//     ss << "\t" << "tt:Z:" << timings;
// #endif
//     }

//     ss << "\n";

//     return ss.str();
}

}