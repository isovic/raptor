/*
 * aligner_factory.cc
 *
 *  Created on: Sep 2, 2017
 *      Author: isovic
 */

#include "aligner/aligner_factory.h"
#include "aligner/aligner_base.h"
#include "aligner/aligner_ksw2_single.h"
#include "aligner/aligner_ksw2_double.h"
#include "aligner/aligner_edlib.h"

namespace raptor {

std::shared_ptr<raptor::AlignerBase> createAligner(raptor::AlignerType aln_type,
                                                   const raptor::AlignmentOptions &opt) {
    switch (aln_type) {
    case raptor::AlignerType::KSW2Single:
        return raptor::createAlignerKSW2Single(opt);
        break;
    case raptor::AlignerType::KSW2Double:
        return raptor::createAlignerKSW2Double(opt);
        break;
    case raptor::AlignerType::Edlib:
        return createAlignerEdlib(opt);
        break;

    default:
        break;
    }

    return createAlignerKSW2Single(opt);
}

} /* namespace raptor */
