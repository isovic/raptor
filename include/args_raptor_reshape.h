/*
 * args_raptor_reshape.h
 *
 *  Created on: Mar 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_RAPTOR_RESHAPE_H_
#define SRC_ARGS_RAPTOR_RESHAPE_H_

#include <cstdint>
#include <memory>
#include <string>
#include <raptor_reshape/params_raptor_reshape.h>

namespace raptor {

int ProcessArgsRaptorReshape(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorReshape> parametert);

} /* namespace raptor */

#endif /* SRC_ARGS_H_ */
