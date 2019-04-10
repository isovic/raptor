/*
 * sequence.h
 *
 *  Created on: Jan 13, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_SERIALIZER_H_
#define SRC_SEQUENCES_SEQUENCE_SERIALIZER_H_

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <index/sequence.h>
#include <index/sequence_file.h>

namespace mindex {

class SequenceSerializer {
    public:
    static bool SerializeSequence(std::ostream& ofs, const mindex::SequencePtr& seq, mindex::SequenceFormat to_fmt);
    static bool SerializeSequenceToFasta(std::ostream& ofs, const mindex::SequencePtr& seq);
    static bool SerializeSequenceToFastq(std::ostream& ofs, const mindex::SequencePtr& seq);
    static bool SerializeSequenceFile(std::ostream& ofs, const mindex::SequenceFilePtr& seq_file, mindex::SequenceFormat to_fmt);
};

}

#endif
