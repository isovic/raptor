/*
 * index_factory.cc
 *
 *  Created on: Sep 26, 2019
 *      Author: Ivan Sovic
 */

#include <index/index_factory.h>
#include <index/minimizer_index.h>
#include <index/dense_index.h>
#include <log/log_tools.h>

namespace mindex {

std::unique_ptr<mindex::IndexBase> createIndex(mindex::IndexType index_type, std::shared_ptr<mindex::IndexParams> params) {
    switch (index_type) {
    case mindex::IndexType::Minimizer:
        return mindex::createMinimizerIndex(params);
        break;
    case mindex::IndexType::Dense:
        return mindex::createDenseIndex(params);
        break;
    default:
        FATAL_REPORT(ERR_UNEXPECTED_VALUE, "Unknown index type.");
        break;
    }
    return nullptr;
}

} /* namespace mindex */
