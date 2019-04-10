/*
 * parameters_raptor.cc
 *
 *  Created on: May 30, 2017
 *      Author: isovic
 */

#include <params/params_raptor.h>

namespace raptor {

std::shared_ptr<raptor::ParamsRaptor> createParamsRaptor() {
    return std::shared_ptr<raptor::ParamsRaptor>(new raptor::ParamsRaptor());
}

} /* namespace raptor */
