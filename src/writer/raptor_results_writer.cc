#include <writer/raptor_results_writer.h>
#include <writer/output_formatter.h>

namespace raptor {

std::unique_ptr<raptor::RaptorResultsWriter> createRaptorResultsWriter(std::ostream& oss, const mindex::IndexPtr index, OutputFormat outfmt) {
    return std::unique_ptr<raptor::RaptorResultsWriter>(new raptor::RaptorResultsWriter(oss, index, outfmt));
}

RaptorResultsWriter::RaptorResultsWriter(std::ostream& oss,
                                        const mindex::IndexPtr index,
                                        raptor::OutputFormat outfmt)
                                        : oss_(oss), index_(index),
                                        outfmt_(outfmt) {

}

RaptorResultsWriter::~RaptorResultsWriter() {

}

void RaptorResultsWriter::WriteHeader(const mindex::HeaderGroupType header_groups) {
    if (outfmt_ == raptor::OutputFormat::SAM) {
        oss_ << "@HD\tVN:1.5" << std::endl;

        for (const auto& it_field: header_groups) {
            oss_ << "@" << it_field.first;
            for (const auto& it_ids: it_field.second) {
                for (const auto& it_tags: it_ids.second) {
                    oss_ << "\t" << it_tags.first << ":" << it_tags.second;
                }
            }
            oss_ << "\n";
        }

        for (size_t i = 0; i < index_->seqs()->size(); ++i) {
            std::string header = TrimToFirstSpace(index_->header(i));
            oss_ << "@SQ\tSN:" << header << "\tLN:" << index_->len(i) << std::endl;
        }

    } else if (outfmt_ == raptor::OutputFormat::GFA2) {
        oss_ << "H\tVN:Z:2.0" << std::endl;

        for (size_t i = 0; i < index_->seqs()->size(); ++i) {
            std::string header = TrimToFirstSpace(index_->header(i));
            oss_ << "S\t" << header << "\t" << index_->len(i) << "\t" << "*" << std::endl;
        }
    }
}

void RaptorResultsWriter::Write(const mindex::SequenceFilePtr seqs, const std::vector<RaptorResults>& results, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    for (auto& result: results) {
        WriteSingleResult(seqs, result, is_alignment_applied, write_custom_tags, one_hit_per_target);
    }
    oss_.flush();
}

void RaptorResultsWriter::WriteSingleResult(const mindex::SequenceFilePtr seqs, const RaptorResults& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
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

    // Write the header for the current sequences.
    if (outfmt_ == raptor::OutputFormat::GFA2) {
        for (const auto& seq: seqs->seqs()) {
            std::string header = TrimToFirstSpace(seq->header());
            oss_ << "S\t" << header << "\t" << seq->len() << "\t" << "*" << std::endl;
        }
    }

    // The writing code is generic.
    if (do_output && !regions_to_write.empty()) {
        for (size_t i = 0; i < regions_to_write.size(); i++) {
            auto aln = regions_to_write[i];
            bool is_secondary = (aln->PathId() > 0);
            bool is_supplementary = (aln->SegmentId() > 0);
            auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());

            if (outfmt_ == raptor::OutputFormat::SAM) {
                oss_ << OutputFormatter::ToSAM(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::PAF) {
                oss_ << OutputFormatter::ToPAF(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::GFA2) {
                oss_ << OutputFormatter::ToGFA2Edge(index_, qseq, aln, mapq, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::MHAP) {
                oss_ << OutputFormatter::ToMHAP(index_, qseq, aln, mapq);
            } else if (outfmt_ == raptor::OutputFormat::M4) {
                oss_ << OutputFormatter::ToM4(index_, qseq, aln, mapq);
            }
        }
    } else {
        // If the q_id_in_batch < 0, it means that the result was initialized, but never
        // updated.
        // This can happen in processing a certain range of reads, and all the other ones
        // outside this range will not be processed, but could have still been allocated.
        if (result.q_id_in_batch >= 0) {
            const auto& qseq = seqs->GetSeqByID(result.q_id_in_batch);

            if (outfmt_ == raptor::OutputFormat::SAM) {
                oss_ << OutputFormatter::UnmappedSAM(qseq, write_custom_tags);
            } else if (outfmt_ == raptor::OutputFormat::PAF) {
                // PAF simply doesn't report unmapped alignments.
            } else if (outfmt_ == raptor::OutputFormat::GFA2) {
                // GFA2 simply doesn't report unmapped alignments.
            } else if (outfmt_ == raptor::OutputFormat::MHAP) {
                // MHAP simply doesn't report unmapped alignments.
            } else if (outfmt_ == raptor::OutputFormat::M4) {
                // M4 simply doesn't report unmapped alignments.
            }
        }
    }
}

}
