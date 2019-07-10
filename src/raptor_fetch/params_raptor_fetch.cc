/*
 * params_raptor_fetch.cc
 *
 *  Created on: Jan 15, 2019
 *      Author: Ivan Sovic
 */

#include <raptor_fetch/params_raptor_fetch.h>

namespace raptor {

std::shared_ptr<raptor::ParamsRaptorFetch> createParamsRaptorFetch() {
    return std::shared_ptr<raptor::ParamsRaptorFetch>(new raptor::ParamsRaptorFetch);
}

ParamsRaptorFetch::ParamsRaptorFetch()
    : verbose_level(0),
      command_line(),
      debug_qid(-1),
      debug_qname(""),
      rdb_path(),
      in_paths(),
      out_prefix(),
      min_cov(0),
      min_len(0),
      min_score(0.0),
      min_idt(0.0),
      use_id_for_output(false),
      job_str("query")

{}

}  // namespace raptor
