/*
 * raptor_results_writer_bam.h
 *
 *  Created on: Jun 16, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RESULTS_WRITER_BAM_H_
#define SRC_RAPTOR_RESULTS_WRITER_BAM_H_

#ifdef RAPTOR_COMPILED_WITH_PBBAM

#include <writer/raptor_results_writer_base.h>
#include <stdint.h>
#include <memory>
#include <ostream>
#include <pbbam/BamWriter.h>

namespace raptor {

class RaptorResultsWriterBAM;

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterBAM(const std::string& out_fn, const mindex::IndexPtr index);

class RaptorResultsWriterBAM : RaptorResultsWriterBase {
 public:
  friend std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterBAM(const std::string& out_fn, const mindex::IndexPtr index);

  ~RaptorResultsWriterBAM();

  void WriteHeader(const mindex::HeaderGroupType header_groups);
  void WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<RaptorResults>& results, bool is_alignment_applied, bool write_custom_tag, bool one_hit_per_targets);
  void WriteSingleResult(const mindex::SequenceFilePtr seqs, const RaptorResults& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target);

 private:
  RaptorResultsWriterBAM(const std::string& out_fn, const mindex::IndexPtr index);
  RaptorResultsWriterBAM(const RaptorResultsWriterBAM&) = delete;
  RaptorResultsWriterBAM& operator=(const RaptorResultsWriterBAM&) = delete;

  std::unique_ptr<PacBio::BAM::BamWriter> OpenBAMWriter_(const std::string out_fn, PacBio::BAM::BamHeader& header) const;
  PacBio::BAM::BamRecord ToBAM(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping,
                            int32_t mapq, bool write_custom_tags,
                            const std::string& timings);
  PacBio::BAM::BamRecord ToUnmappedBAM(const mindex::SequencePtr& qseq);

  std::string out_fn_;
  std::unique_ptr<PacBio::BAM::BamWriter> bam_writer_;
  const mindex::IndexPtr index_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RESULT_WRITER_H_ */

#endif
