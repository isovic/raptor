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
    : subprogram(),
      command_line(),
      verbose_level(0),
      in_path()

{}

}  // namespace raptor
