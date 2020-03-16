/*
 * raptor_results_writer_stream.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RESULTS_WRITER_CSTREAM_H_
#define SRC_RAPTOR_RESULTS_WRITER_CSTREAM_H_

#include <writer/raptor_results_writer_base.h>
#include <utility/files.hpp>
#include <stdint.h>
#include <memory>

namespace raptor {

class RaptorResultsWriterCStream;

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(FILE* fp_out, const mindex::IndexPtr index, OutputFormat outfmt);
std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(const std::string& out_fn, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

class RaptorResultsWriterCStream : RaptorResultsWriterBase {
 public:
  friend std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(FILE* fp_out, const mindex::IndexPtr index, raptor::OutputFormat outfmt);
  friend std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriterCStream(const std::string& out_fn, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

  ~RaptorResultsWriterCStream();

  void WriteHeader(const mindex::HeaderGroupType header_groups);
  void WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<std::unique_ptr<raptor::RaptorResults>>& results, bool is_alignment_applied, bool write_custom_tag, bool one_hit_per_targets);
  void WriteSingleResult(const mindex::SequenceFilePtr seqs, const std::unique_ptr<raptor::RaptorResults>& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target);

 private:
  RaptorResultsWriterCStream(FILE* fp_out, bool do_fclose, const mindex::IndexPtr index, raptor::OutputFormat outfmt);
  RaptorResultsWriterCStream(const RaptorResultsWriterCStream&) = delete;
  RaptorResultsWriterCStream& operator=(const RaptorResultsWriterCStream&) = delete;

  FILE *fp_out_;
  bool should_fclose_;
  const mindex::IndexPtr index_;
  raptor::OutputFormat outfmt_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RESULT_WRITER_CSTREAM_H_ */
