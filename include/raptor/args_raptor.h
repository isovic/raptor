/*
 * args_raptor.h
 *
 *  Created on: Mar 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_RAPTOR_H_
#define SRC_ARGS_RAPTOR_H_

#include <cstdint>
#include <memory>
#include <string>
#include <params/params_raptor.h>

namespace raptor {

void VerboseShortHelpRaptor(int argc, char **argv);
void VerboseShortHelpRaptorAndExit(int argc, char **argv, int ret_val = 1);
int ProcessArgsRaptor(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptor> parameters);

} /* namespace raptor */

#endif /* SRC_ARGS_H_ */
