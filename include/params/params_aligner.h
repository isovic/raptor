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

using ParamsAlignerPtr = std::shared_ptr<raptor::ParamsAligner>;

class ParamsAligner {
   public:
    friend std::shared_ptr<raptor::ParamsAligner> createParamsAligner();
    ~ParamsAligner() = default;

    int64_t verbose_level = 0;
    int64_t debug_qid = -1;
    std::string debug_qname = "";

    // Contains match, mismatch and gap penalties, bandwidth and more.
    raptor::AlignmentOptions aligner_opt;
    double min_identity = 65.0;
    double max_evalue = 1e0;
    bool use_basic_cigar = false;
    raptor::AlignerType aligner_type = raptor::AlignerType::Edlib;

    bool no_extend_alignment = false;
    bool is_rna = false;
    bool relabel_secondary_supp = true;
    double min_secondary_to_primary_ratio = 0.80;
    int32_t allowed_suppl_overlap = 0;

   private:
    ParamsAligner() = default;
    ParamsAligner(const ParamsAligner&) = delete;
    ParamsAligner& operator=(const ParamsAligner&) = delete;
};

inline std::shared_ptr<raptor::ParamsAligner> createParamsAligner() {
    return std::shared_ptr<raptor::ParamsAligner>(new raptor::ParamsAligner);
}

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_ALIGNER_H_ */
