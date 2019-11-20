/*
 * args_raptor_sova.h
 *
 *  Created on: Mar 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_RAPTOR_HOLY_H_
#define SRC_ARGS_RAPTOR_HOLY_H_

#include <cstdint>
#include <memory>
#include <string>
#include <raptor_sova/params_raptor_sova.h>

namespace raptor {
namespace sova {

int ProcessArgsRaptorSova(int argc, char **argv, std::shared_ptr<raptor::sova::ParamsRaptorSova> parameters);

} /* namespace sova */
} /* namespace raptor */

#endif /* SRC_ARGS_H_ */
