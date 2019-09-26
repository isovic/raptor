/*
 * index_factory.h
 *
 *  Created on: Sep 26, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_INDEX_INDEX_FACTORY_H_
#define SRC_INDEX_INDEX_FACTORY_H_

#include <index/index_base.h>
#include <index/index_params.h>
#include <index/index_types.h>
#include <memory>

namespace mindex {

/** @brief A factory function for a concrete index object.
 *
 */
std::unique_ptr<mindex::IndexBase> createIndex(mindex::IndexType index_type, std::shared_ptr<mindex::IndexParams> params);

} /* namespace mindex */

#endif /* SRC_ALIGNER_ALIGNER_FACTORY_H_ */
