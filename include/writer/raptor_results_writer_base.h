/*
 * raptor_results_writer_base.h
 *
 *  Created on: Jun 16, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_RESULTS_WRITER_BASE_H_
#define SRC_RAPTOR_RESULTS_WRITER_BASE_H_

#include <containers/raptor_results.h>
#include <raptor/yield_index.h>
#include <params/params_raptor.h>
#include <sequences/sequence_file.h>

namespace raptor {

class RaptorResultsWriterBase;

using RaptorResultsWriterBasePtr = std::unique_ptr<raptor::RaptorResultsWriterBase>;

class RaptorResultsWriterBase {
 public:
  virtual ~RaptorResultsWriterBase() { }

  virtual void WriteHeader(const mindex::HeaderGroupType header_groups) = 0;
  virtual void WriteBatch(const mindex::SequenceFilePtr seqs, const std::vector<std::unique_ptr<raptor::RaptorResults>>& results, bool is_alignment_applied, bool write_custom_tag, bool one_hit_per_targets) = 0;
  virtual void WriteSingleResult(const mindex::SequenceFilePtr seqs, const std::unique_ptr<raptor::RaptorResults>& result, bool is_alignment_applied, bool write_custom_tags, bool one_hit_per_target) = 0;
};

} /* namespace raptor */

#endif /* SRC_RAPTOR_RESULT_WRITER_BASE_H_ */
