/*
 * params_aligner.cc
 *
 *  Created on: Mar 04, 2018
 *      Author: Ivan Sovic
 */

#include <params/params_aligner.h>

namespace raptor {

std::shared_ptr<raptor::ParamsAligner> createParamsAligner() {
    return std::shared_ptr<raptor::ParamsAligner>(new raptor::ParamsAligner);
}

ParamsAligner::ParamsAligner()
    :   verbose_level(0),
        debug_qid(-1),
        debug_qname(""),

        // Contains match, mismatch and gap penalties, bandwidth and more.
        aligner_opt(),
        min_identity(65.0),
        max_evalue(1e0),
        use_basic_cigar(false),
        aligner_type(raptor::AlignerType::Edlib),

        no_extend_alignment(false),
        is_rna(false)

      {

}

}  // namespace raptor
