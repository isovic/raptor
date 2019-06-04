/*
 * args_raptor_fetch.h
 *
 *  Created on: Jun 04, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_RAPTOR_FETCH_H_
#define SRC_ARGS_RAPTOR_FETCH_H_

#include <cstdint>
#include <memory>
#include <string>
#include <raptor_fetch/params_raptor_fetch.h>

namespace raptor {

int ProcessArgsRaptorFetch(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorFetch> parameters);

} /* namespace raptor */

#endif /* SRC_ARGS_H_ */
