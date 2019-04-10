/*
 * sequence_serializer.cc
 *
 *  Created on: Jan 13, 2019
 *      Author: Ivan Sovic
 */

#include <index/sequence_serializer.h>

namespace mindex {

bool SequenceSerializer::SerializeSequence(std::ostream& ofs, const mindex::SequencePtr& seq, mindex::SequenceFormat to_fmt) {
    if (seq == nullptr) {
        return false;
    }
    if (to_fmt == mindex::SequenceFormat::Fasta) {
        return SerializeSequenceToFasta(ofs, seq);
    } else if (to_fmt == mindex::SequenceFormat::Fastq) {
        return SerializeSequenceToFastq(ofs, seq);

    } else if (to_fmt == mindex::SequenceFormat::SAM) {

    } else if (to_fmt == mindex::SequenceFormat::GFA) {

    } else if (to_fmt == mindex::SequenceFormat::GFA1) {

    } else if (to_fmt == mindex::SequenceFormat::GFA2) {

    } else {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Unsupported SequenceFormat. Skipping.");
    }

    return false;
}

bool SequenceSerializer::SerializeSequenceToFasta(std::ostream& ofs, const mindex::SequencePtr& seq) {
    if (seq == nullptr) {
        return false;
    }
    ofs << ">" << seq->header() << std::endl;
    ofs.write(reinterpret_cast<const char *>(&seq->data()[0]), seq->data().size());
    ofs << std::endl;
    return true;
}

bool SequenceSerializer::SerializeSequenceToFastq(std::ostream& ofs, const mindex::SequencePtr& seq) {
    if (seq == nullptr) {
        return false;
    }
    if (seq->data().size() != seq->qual().size()) {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Sequence '%s' does not have the same length of the data and qual arrays. Skipping.\n", seq->header().c_str());
        return false;
    }
    ofs << "@" << seq->header() << std::endl;
    ofs.write(reinterpret_cast<const char *>(&seq->data()[0]), seq->data().size());
    ofs << std::endl;
    ofs << "+" << std::endl;
    ofs.write(reinterpret_cast<const char *>(&seq->qual()[0]), seq->qual().size());
    ofs << std::endl;
    return true;
}

bool SequenceSerializer::SerializeSequenceFile(std::ostream& ofs, const mindex::SequenceFilePtr& seq_file, mindex::SequenceFormat to_fmt) {
    for (size_t i = 0; i < seq_file->seqs().size(); ++i) {
        const auto& seq = seq_file->seqs()[i];
        bool rv = SerializeSequence(ofs, seq, to_fmt);
        if (rv == false) {
            WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Could not serialize sequence file at sequence %ld", seq->header().c_str());
            return false;
        }
    }
    return true;
}

}
