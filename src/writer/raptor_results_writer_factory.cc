/*
 * aligner_factory.cc
 *
 *  Created on: Sep 2, 2017
 *      Author: isovic
 */

#include <writer/raptor_results_writer_factory.h>
#include <writer/raptor_results_writer_base.h>
#include <writer/raptor_results_writer_stream.h>
#include <writer/raptor_results_writer_bam.h>
#include <log/log_tools.h>

namespace raptor {

std::unique_ptr<raptor::RaptorResultsWriterBase> createRaptorResultsWriter(const std::string& out_fn, const mindex::IndexPtr index, raptor::OutputFormat out_fmt) {
    switch (out_fmt) {
#ifdef RAPTOR_COMPILED_WITH_PBBAM
    case raptor::OutputFormat::BAM:
        return createRaptorResultsWriterBAM(out_fn, index);
        break;
#endif
    case raptor::OutputFormat::Unknown:
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Unknown output file format.");
        break;
    default:
        return createRaptorResultsWriterStream(out_fn, index, out_fmt);
        break;
    }
    return nullptr;
}

} /* namespace raptor */
