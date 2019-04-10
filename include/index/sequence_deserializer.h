/*
 * sequence_deserializer.h
 *
 *  Created on: Jan 16, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_DESERIALIZER_H_
#define SRC_SEQUENCES_SEQUENCE_DESERIALIZER_H_

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <index/sequence.h>
#include <index/sequence_file_enums.h>
#include <index/sequence_file_handlers.h>

namespace mindex {

class SequenceDeserializer {
  public:
    static mindex::SequencePtr DeserializeSequence(mindex::SequenceFileHandlersPtr& fp_handler, mindex::SequenceFormat in_fmt, bool convert_to_uppercase);
    static mindex::SequencePtr DeserializeSequenceFromFastx(mindex::SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase);
    static mindex::SequencePtr DeserializeSequenceFromGFA1(mindex::SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase);
    static mindex::SequencePtr DeserializeSequenceFromGFA2(mindex::SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase);
};

}

#endif
