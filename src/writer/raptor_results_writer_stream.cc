#include <writer/raptor_results_writer_stream.h>
#include <writer/output_formatter.h>
#include <fstream>
#include <log/log_tools.h>

namespace raptor {

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterStream(std::shared_ptr<std::ostream>& oss_ptr, const mindex::IndexPtr index, OutputFormat outfmt) {
    return std::unique_ptr<raptor::RaptorResultsWriterBase>(new raptor::RaptorResultsWriterStream(oss_ptr, index, outfmt));
}

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterStream(const std::string& out_fn, const mindex::IndexPtr index, OutputFormat outfmt) {
    std::shared_ptr<std::ostream> oss_ptr(&std::cout, [](void*) {});
    if (out_fn.size() > 0 && out_fn != "-") {
        oss_ptr = std::shared_ptr<std::ostream>(new std::ofstream(out_fn));
    }
    return std::unique_ptr<raptor::RaptorResultsWriterBase>(new raptor::RaptorResultsWriterStream(oss_ptr, index, outfmt));
}

RaptorResultsWriterStream::RaptorResultsWriterStream(std::shared_ptr<std::ostream>& oss_ptr,
                                        const mindex::IndexPtr index,
                                        raptor::OutputFormat outfmt)
                                        : oss_ptr_(oss_ptr), index_(index),
                                        outfmt_(outfmt) {

}

RaptorResultsWriterStream::~RaptorResultsWriterStream() {

}

void RaptorResultsWriterStream::WriteHeader(const mindex::HeaderGroupType header_groups) {
    if (outfmt_ == raptor::OutputFormat::SAM) {
        *oss_ptr_ << "@HD\tVN:1.5" << std::endl;

        for (const auto& it_field: header_groups) {
            for (const auto& it_ids: it_field.second) {
                *oss_ptr_ << "@" << it_field.first;
                for (const auto& it_tags: it_ids.second) {
                    *oss_ptr_ << "\t" << it_tags.name << ":" << it_tags.val;
                }
                *oss_ptr_ << "\n";
            }
        }

        if (index_->seqs() != nullptr) {
            for (size_t i = 0; i < index_->seqs()->size(); ++i) {
                std::string header = TrimToFirstSpace(index_->header(i));
                *oss_ptr_ << "@SQ\tSN:" << header << "\tLN:" << index_->len(i) << std::endl;
            }
        }
    } else if (outfmt_ == raptor::OutputFormat::GFA2) {
        *oss_ptr_ << "H\tVN:Z:2.0" << std::endl;

        if (index_->seqs() != nullptr) {
            for (size_t i = 0; i < index_->seqs()->size(); ++i) {
                std::string header = TrimToFirstSpace(index_->header(i));
                *oss_ptr_ << "S\t" << header << "\t" << index_->len(i) << "\t" << "*" << std::endl;
            }
        }
    }
}

void RaptorResultsWriterStream::WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<std::unique_ptr<raptor::RaptorResults>>& results, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    for (auto& result: results) {
        WriteSingleResult(seqs, result, is_alignment_applied, write_custom_tags, one_hit_per_target);
    }
    oss_ptr_->flush();
}

void RaptorResultsWriterStream::WriteSingleResult(const mindex::SequenceFilePtr seqs, const std::unique_ptr<raptor::RaptorResults>& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (result == nullptr) {
        return;
    }

    bool do_output = !result->regions().empty();
    std::string timings_all = OutputFormatter::TimingMapToString(result->timings());
    const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result->regions();
    int64_t q_id = result->q_id();

    // The writing code is generic.
    if (do_output && !regions_to_write.empty()) {
        for (size_t i = 0; i < regions_to_write.size(); i++) {
            auto aln = regions_to_write[i];
            bool is_secondary = (aln->PathId() > 0);
            bool is_supplementary = (aln->SegmentId() > 0);
            auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());

            if (outfmt_ == raptor::OutputFormat::SAM) {
                *oss_ptr_ << OutputFormatter::ToSAM(index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::PAF) {
                *oss_ptr_ << OutputFormatter::ToPAF(index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::GFA2) {
                *oss_ptr_ << OutputFormatter::ToGFA2Edge(index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::MHAP) {
                *oss_ptr_ << OutputFormatter::ToMHAP(index_, qseq, aln);
            } else if (outfmt_ == raptor::OutputFormat::M4) {
                *oss_ptr_ << OutputFormatter::ToM4(index_, qseq, aln);
            }
        }
    } else {
        // If the q_id < 0, it means that the result was initialized, but never
        // updated.
        // This can happen in processing a certain range of reads, and all the other ones
        // outside this range will not be processed, but could have still been allocated.
        if (q_id >= 0) {
            const auto& qseq = seqs->GetSeqByAbsID(q_id);

            if (outfmt_ == raptor::OutputFormat::SAM) {
                *oss_ptr_ << OutputFormatter::UnmappedSAM(qseq, write_custom_tags);
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
