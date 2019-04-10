/*
 * args.h
 *
 *  Created on: Mar 07, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_ARGS_H_
#define SRC_ARGS_H_

#include <cstdint>
#include <memory>
#include <string>
#include <params/params_raptor.h>
#include <params/params_raptor_index.h>
#include <raptor_reshape/params_raptor_reshape.h>
#include <raptor_fetch/params_raptor_fetch.h>

namespace raptor {

void VerboseShortHelpAndExit(int argc, char **argv);
int ProcessArgsRaptor(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptor> parameters);
int ProcessArgsRaptorIndex(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorIndex> parameters);
int ProcessArgsRaptorReshape(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorReshape> parametert);
int ProcessArgsRaptorFetch(int argc, char **argv, std::shared_ptr<raptor::ParamsRaptorFetch> parametert);

} /* namespace raptor */

#endif /* SRC_ARGS_H_ */
