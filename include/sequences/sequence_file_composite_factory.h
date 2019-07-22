/*
 * sequence_file_composite_factory.
 *
 *  Created on: Jul 22, 2019
 *      Author: isovic
 */

#ifndef SRC_INDEX_SEQUENCE_FILE_COMPOSITE_FACTORY_H_
#define SRC_INDEX_SEQUENCE_FILE_COMPOSITE_FACTORY_H_

#include <memory>
#include <vector>

#include <sequences/sequence_file_composite_base.h>
#include <sequences/sequence_file_composite_fofn.h>
#include <sequences/sequence_file_composite_pbxml.h>

namespace mindex {

mindex::SequenceFileCompositeBasePtr createSequenceFileCompositeFactory(const std::vector<std::string>& in_paths, mindex::SequenceFormat in_fmt);

}

#endif
