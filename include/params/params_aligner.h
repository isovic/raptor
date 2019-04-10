/*
 * params_aligner.h
 *
 *  Created on: Mar 04, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PARAMS_ALIGNER_H_
#define SRC_RAPTOR_PARAMS_ALIGNER_H_

#include <cstdint>
#include <memory>
#include <string>
#include <aligner/alignment_options.h>
#include <aligner/aligner_factory.h>

namespace raptor {

class ParamsAligner;

typedef std::shared_ptr<raptor::ParamsAligner> ParamsAlignerPtr;

std::shared_ptr<raptor::ParamsAligner> createParamsAligner();

class ParamsAligner {
   public:
    friend std::shared_ptr<raptor::ParamsAligner> createParamsAligner();
    ~ParamsAligner() = default;

    int64_t verbose_level;
    int64_t debug_qid;
    std::string debug_qname;

    // Contains match, mismatch and gap penalties, bandwidth and more.
    raptor::AlignmentOptions aligner_opt;
    double min_identity;
    double max_evalue;
    bool use_basic_cigar;
    raptor::AlignerType aligner_type;

    bool no_extend_alignment;
    bool is_rna;

   private:
    ParamsAligner();
    ParamsAligner(const ParamsAligner&) = delete;
    ParamsAligner& operator=(const ParamsAligner&) = delete;
};

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_ALIGNER_H_ */
