/*
 * raptor_results_writer.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RESULTS_WRITER_H_
#define SRC_RAPTOR_RESULTS_WRITER_H_

#include <writer/raptor_results_writer_base.h>
#include <stdint.h>
#include <memory>
#include <ostream>

namespace raptor {

class RaptorResultsWriterStream;

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterStream(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

class RaptorResultsWriterStream : RaptorResultsWriterBase {
 public:
  friend std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterStream(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

  ~RaptorResultsWriterStream();

  void WriteHeader(const mindex::HeaderGroupType header_groups);
  void WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<RaptorResults>& results, bool is_alignment_applied, bool write_custom_tag, bool one_hit_per_targets);
  void WriteSingleResult(const mindex::SequenceFilePtr seqs, const RaptorResults& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target);

 private:
  RaptorResultsWriterStream(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);
  RaptorResultsWriterStream(const RaptorResultsWriterStream&) = delete;
  RaptorResultsWriterStream& operator=(const RaptorResultsWriterStream&) = delete;

  std::ostream& oss_;
  const mindex::IndexPtr index_;
  raptor::OutputFormat outfmt_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RESULT_WRITER_H_ */
