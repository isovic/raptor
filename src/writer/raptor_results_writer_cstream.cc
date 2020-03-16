#include <writer/raptor_results_writer_cstream.h>
#include <writer/output_formatter.h>
#include <fstream>
#include <log/log_tools.h>

namespace raptor {

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(FILE* fp_out, const mindex::IndexPtr index, OutputFormat outfmt) {
    return std::unique_ptr<raptor::RaptorResultsWriterBase>(new raptor::RaptorResultsWriterCStream(fp_out, false, index, outfmt));
}

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(const std::string& out_fn, const mindex::IndexPtr index, OutputFormat outfmt) {
    FILE *fp_out = stdout;
    bool do_fclose = false;
    if (out_fn.size() && out_fn != "-") {
        fp_out = fopen(out_fn.c_str(), "w");
        if (fp_out == nullptr) {
            std::ostringstream err_oss;
            err_oss << "RaptorResultsWriterCStream could not open file '" << out_fn << "'!";
            throw std::runtime_error(err_oss.str());
        }
        do_fclose = true;
    }
    return std::unique_ptr<raptor::RaptorResultsWriterBase>(new raptor::RaptorResultsWriterCStream(fp_out, do_fclose, index, outfmt));
}

RaptorResultsWriterCStream::RaptorResultsWriterCStream(FILE* fp_out, bool do_fclose,
                                        const mindex::IndexPtr index,
                                        raptor::OutputFormat outfmt)
                                        : fp_out_(fp_out)
                                        , should_fclose_(do_fclose)
                                        , index_(index)
                                        , outfmt_(outfmt) {

}

RaptorResultsWriterCStream::~RaptorResultsWriterCStream() {
    fclose(fp_out_);
}

void RaptorResultsWriterCStream::WriteHeader(const mindex::HeaderGroupType header_groups) {
    if (outfmt_ == raptor::OutputFormat::SAM) {
        fprintf(fp_out_, "@HD\tVN:1.5\n");

        for (const auto& it_field: header_groups) {
            for (const auto& it_ids: it_field.second) {
                fprintf(fp_out_, "@%s", it_field.first.c_str());
                for (const auto& it_tags: it_ids.second) {
                    fprintf(fp_out_, "\t%s:%s", it_tags.name.c_str(), it_tags.val.c_str());
                }
                fprintf(fp_out_, "\n");
            }
        }

        if (index_->seqs() != nullptr) {
            for (size_t i = 0; i < index_->seqs()->size(); ++i) {
                std::string header = TrimToFirstSpace(index_->header(i));
                fprintf(fp_out_, "@SQ\tSN:%s\tLN:%d\n", header.c_str(), index_->len(i));
            }
        }
    } else if (outfmt_ == raptor::OutputFormat::GFA2) {
        fprintf(fp_out_, "H\tVN:Z:2.0\n");

        if (index_->seqs() != nullptr) {
            for (size_t i = 0; i < index_->seqs()->size(); ++i) {
                std::string header = TrimToFirstSpace(index_->header(i));
                fprintf(fp_out_, "S\t%s\t%d\t*\n", header.c_str(), index_->len(i));
            }
        }
    }
}

void RaptorResultsWriterCStream::WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<std::unique_ptr<raptor::RaptorResults>>& results, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    for (auto& result: results) {
        WriteSingleResult(seqs, result, is_alignment_applied, write_custom_tags, one_hit_per_target);
    }
    fflush(fp_out_);
}

void RaptorResultsWriterCStream::WriteSingleResult(const mindex::SequenceFilePtr seqs, const std::unique_ptr<raptor::RaptorResults>& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) {
    if (result == nullptr) {
        return;
    }

    bool do_output = !result->regions().empty();
    const std::vector<std::shared_ptr<raptor::RegionBase>>& regions_to_write = result->regions();
    int64_t q_id = result->q_id();
    std::string timings_all;
#ifdef RAPTOR_DEBUG_TIMINGS
    // Don't waste time formatting the timings unless the macro is defined.
    timings_all = OutputFormatter::TimingMapToString(result->timings());
#endif

    // The writing code is generic.
    if (do_output && !regions_to_write.empty()) {
        for (size_t i = 0; i < regions_to_write.size(); i++) {
            auto aln = regions_to_write[i];
            bool is_secondary = (aln->PathId() > 0);
            bool is_supplementary = (aln->SegmentId() > 0);
            auto& qseq = seqs->GetSeqByAbsID(aln->QueryID());

            if (outfmt_ == raptor::OutputFormat::SAM) {
                OutputFormatter::ToSAM(fp_out_, index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::PAF) {
                OutputFormatter::ToPAF(fp_out_, index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::GFA2) {
                OutputFormatter::ToGFA2Edge(fp_out_, index_, qseq, aln, write_custom_tags, timings_all);
            } else if (outfmt_ == raptor::OutputFormat::MHAP) {
                OutputFormatter::ToMHAP(fp_out_, index_, qseq, aln);
            } else if (outfmt_ == raptor::OutputFormat::M4) {
                OutputFormatter::ToM4(fp_out_, index_, qseq, aln);
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
                OutputFormatter::UnmappedSAM(fp_out_, qseq, write_custom_tags);
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
