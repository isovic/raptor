/*
 * yield_index.h
 *
 *  Created on: May 30, 2017
 *      Author: isovic
 */

#ifndef SRC_RAPTOR_INDEX_FACTORY_H_
#define SRC_RAPTOR_INDEX_FACTORY_H_

#include <string>
#include <memory>
#include <index/minimizer_index.h>
#include <params/params_raptor.h>

namespace raptor {

std::shared_ptr<mindex::IndexBase> YieldIndex(const std::vector<std::string>& ref_paths, mindex::SequenceFormat ref_fmt,
                                                const std::string &index_path, bool rebuild_index,
                                                bool index_on_the_fly, bool auto_rebuild_index,
                                                int64_t rdb_block_id,
                                                mindex::IndexType index_type,
                                                std::shared_ptr<mindex::ParamsIndex> index_params);

}  // namespace raptor

#endif /* SRC_RAPTOR_INDEX_FACTORY_H_ */
