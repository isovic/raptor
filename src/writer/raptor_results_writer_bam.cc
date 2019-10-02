
#ifdef RAPTOR_COMPILED_WITH_PBBAM

#include <writer/raptor_results_writer_bam.h>
#include <writer/output_formatter.h>
#include <fstream>
#include <log/log_tools.h>
#include <pbbam/BamRecord.h>
#include <pbbam/BamRecordImpl.h>
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

    // Pbbam expects "-" to write to stdout instead of a logical "".
    if (out_fn_.empty()) {
        out_fn_ = "-";
    }

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
        for (const auto& it_ids: it_field.second) {
            sam_header << "@" << it_field.first;
            for (const auto& it_tags: it_ids.second) {
                sam_header << "\t" << it_tags.name << ":" << it_tags.val;
            }
            sam_header << "\n";
        }
    }
    for (size_t i = 0; i < index_->seqs()->size(); ++i) {
        std::string qname = TrimToFirstSpace(index_->header(i));
        sam_header << "@SQ\tSN:" << qname << "\tLN:" << index_->len(i) << std::endl;
    }

    PacBio::BAM::BamHeader header(sam_header.str());
    bam_writer_ = std::move(OpenBAMWriter_(out_fn_, header));

}

void RaptorResultsWriterBAM::WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<std::unique_ptr<raptor::RaptorResults>>& results, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (bam_writer_ == nullptr) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Output BAM file not opened. Ensure that WriteHeader is called first.");
    }

    for (auto& result: results) {
        WriteSingleResult(seqs, result, is_alignment_applied, write_custom_tags, one_hit_per_target);
    }
}

void RaptorResultsWriterBAM::WriteSingleResult(const mindex::SequenceFilePtr seqs, const std::unique_ptr<raptor::RaptorResults>& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (bam_writer_ == nullptr) {
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Output BAM file not opened. Ensure that WriteHeader is called first.");
    }

    if (result == nullptr) {
        return;
    }

    bool do_output = !result->regions().empty();
    std::string timings_all = OutputFormatter::TimingMapToString(result->timings());
    int32_t mapq = result->mapq();
    const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result->regions();
    int64_t q_id = result->q_id();

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
        // If the q_id < 0, it means that the result was initialized, but never
        // updated.
        // This can happen in processing a certain range of reads, and all the other ones
        // outside this range will not be processed, but could have still been allocated.
        if (q_id >= 0) {
            const auto& qseq = seqs->GetSeqByAbsID(q_id);

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

    int32_t path_id = mapping->PathId();
    int32_t num_paths = mapping->PathsNum();
    int32_t segment_in_path = mapping->SegmentId();
    int32_t num_segments_in_path = mapping->SegmentsNum();
    int32_t is_split_mapped = num_segments_in_path > 1;

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
        PacBio::BAM::BamRecordImpl& impl = record.Impl();
        impl.Flag(flag);

        if (write_custom_tags) {
            impl.AddTag("NM", PacBio::BAM::Tag(edit_dist));
            impl.AddTag("AS", PacBio::BAM::Tag(score));
            impl.AddTag("pi", PacBio::BAM::Tag(path_id));
            impl.AddTag("pj", PacBio::BAM::Tag(segment_in_path));
            impl.AddTag("pn", PacBio::BAM::Tag(num_segments_in_path));
            impl.AddTag("ps", PacBio::BAM::Tag(is_split_mapped));
            #ifdef RAPTOR_DEBUG_TIMINGS
                impl.AddTag("tt", PacBio::BAM::Tag(timings));
            #endif
        }
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
        PacBio::BAM::BamRecordImpl& impl = record.Impl();
        impl.Flag(flag);

        if (write_custom_tags) {
            impl.AddTag("NM", PacBio::BAM::Tag(edit_dist));
            impl.AddTag("AS", PacBio::BAM::Tag(score));
            impl.AddTag("pi", PacBio::BAM::Tag(path_id));
            impl.AddTag("pj", PacBio::BAM::Tag(segment_in_path));
            impl.AddTag("pn", PacBio::BAM::Tag(num_segments_in_path));
            impl.AddTag("ps", PacBio::BAM::Tag(is_split_mapped));
            #ifdef RAPTOR_DEBUG_TIMINGS
                impl.AddTag("tt", PacBio::BAM::Tag(timings));
            #endif
        }
    }

    return record;
}

}

#endif
