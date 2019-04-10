/*
 * aligner_factory.h
 *
 *  Created on: Sep 2, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNER_FACTORY_H_
#define SRC_ALIGNER_ALIGNER_FACTORY_H_

#include <memory>
#include <aligner/aligner_base.h>

namespace raptor {

enum class AlignerType { Edlib, KSW2Single, KSW2Double, Unknown };

/** @brief A factory function for a concrete aligner object.
 *
 */
std::shared_ptr<raptor::AlignerBase> createAligner(raptor::AlignerType aln_type,
                                                   const raptor::AlignmentOptions &opt);

} /* namespace raptor */

#endif /* SRC_ALIGNER_ALIGNER_FACTORY_H_ */
