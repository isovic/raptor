/*
 * raptor_results_writer_factory.h
 *
 *  Created on: Jun 16, 2019
 *      Author: isovic
 */

#ifndef SRC_WRITER_RAPTOR_RESULTS_WRITER_FACTORY_H_
#define SRC_WRITER_RAPTOR_RESULTS_WRITER_FACTORY_H_

#include <writer/raptor_results_writer_base.h>
#include <memory>

namespace raptor {

/** @brief A factory function for a concrete writer object.
 *
 */
std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriter(const std::string& out_fn, const mindex::IndexPtr index, raptor::OutputFormat out_fmt);

} /* namespace raptor */

#endif /* SRC_ALIGNER_ALIGNER_FACTORY_H_ */
