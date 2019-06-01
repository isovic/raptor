/*
 * raptor_results_writer.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RESULTS_WRITER_H_
#define SRC_RAPTOR_RESULTS_WRITER_H_

#include <stdint.h>
#include <memory>
#include <containers/raptor_results.h>
#include <raptor/index_factory.h>
#include <params/params_raptor.h>
#include <ostream>
#include <sequences/sequence_file.h>

namespace raptor {

class RaptorResultsWriter;

std::unique_ptr<raptor::RaptorResultsWriter> createRaptorResultsWriter(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

class RaptorResultsWriter {
 public:
  friend std::unique_ptr<raptor::RaptorResultsWriter> createRaptorResultsWriter(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);

  ~RaptorResultsWriter();

  void WriteHeader(const mindex::HeaderGroupType header_groups);
  void Write(const mindex::SequenceFilePtr seqs, const std::vector<RaptorResults>& results, bool is_alignment_applied, bool write_custom_tag, bool one_hit_per_targets);
  void WriteSingleResult(const mindex::SequenceFilePtr seqs, const RaptorResults& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target);

 private:
  RaptorResultsWriter(std::ostream& oss, const mindex::IndexPtr index, raptor::OutputFormat outfmt);
  RaptorResultsWriter(const RaptorResultsWriter&) = delete;
  RaptorResultsWriter& operator=(const RaptorResultsWriter&) = delete;

  std::ostream& oss_;
  const mindex::IndexPtr index_;
  raptor::OutputFormat outfmt_;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RESULT_WRITER_H_ */
