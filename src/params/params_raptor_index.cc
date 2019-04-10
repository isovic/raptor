/*
 * params_raptor_index.cc
 *
 *  Created on: Dec 13, 2018
 *      Author: Ivan Sovic
 */

#include <params/params_raptor_index.h>

namespace raptor {

std::shared_ptr<raptor::ParamsRaptorIndex> createParamsRaptorIndex() {
    return std::shared_ptr<raptor::ParamsRaptorIndex>(new raptor::ParamsRaptorIndex);
}

ParamsRaptorIndex::ParamsRaptorIndex()
    : subprogram(),
      command_line(),
      index_params(mindex::createIndexParams()),
      verbose_level(0),
      ref_paths(),
    //   ref_fmt(SequenceFormat::Unknown),
      out_path(),
      batch_size(400),
      num_threads(1),
      keep_lowercase(false)

{}

}  // namespace raptor
