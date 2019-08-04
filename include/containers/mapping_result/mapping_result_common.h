/*
 * mapping_result_common.h
 *
 *  Created on: Aug 04, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_MAPPING_RESULT_COMMON_H_
#define SRC_CONTAINERS_MAPPING_RESULT_COMMON_H_

#include <memory>
#include <vector>
#include <containers/mapping_result/mapping_result_base.h>
#include <containers/region/region_base.h>
#include <containers/region/region_aligned.h>
#include <containers/region/region_mapped.h>

namespace raptor {

void RelabelSupplementary(std::vector<std::shared_ptr<raptor::RegionBase>>& ret);

}

#endif
