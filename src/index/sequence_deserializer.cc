/*
 * sequence_deserializer.cc
 *
 *  Created on: Jan 16, 2019
 *      Author: Ivan Sovic
 */

#include <index/sequence_deserializer.h>
#include <index/sequence.h>
#include <log/log_tools.h>
#include <utility/files.hpp>
#include <sstream>

namespace mindex {

mindex::SequencePtr SequenceDeserializer::DeserializeSequence(mindex::SequenceFileHandlersPtr& fp_handler, mindex::SequenceFormat in_fmt, bool convert_to_uppercase) {

    mindex::SequencePtr ret;

    if (in_fmt == mindex::SequenceFormat::Fasta) {
        ret = SequenceDeserializer::DeserializeSequenceFromFastx(fp_handler, convert_to_uppercase);

    } else if (in_fmt == mindex::SequenceFormat::Fastq) {
        ret = SequenceDeserializer::DeserializeSequenceFromFastx(fp_handler, convert_to_uppercase);

    } else if (in_fmt == mindex::SequenceFormat::SAM) {

    } else if (in_fmt == mindex::SequenceFormat::GFA) {

    } else if (in_fmt == mindex::SequenceFormat::GFA1) {
        ret = SequenceDeserializer::DeserializeSequenceFromGFA1(fp_handler, convert_to_uppercase);

    } else if (in_fmt == mindex::SequenceFormat::GFA2) {
        ret = SequenceDeserializer::DeserializeSequenceFromGFA2(fp_handler, convert_to_uppercase);

    } else {
        WARNING_REPORT(ERR_UNEXPECTED_VALUE, "Unsupported SequenceFormat. Skipping.");
        return nullptr;
    }

    return ret;
}

mindex::SequencePtr SequenceDeserializer::DeserializeSequenceFromFastx(SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase) {
    int64_t id = -1;
    int64_t abs_id = -1;

    int32_t l = kseq_read(fp_handler->fp_kseq);

    if (l < 0) {
        return nullptr;
    }

    auto seq = mindex::createSequence();

    std::string header;
    if (fp_handler->fp_kseq->name.l > 0) {
        header += std::string(fp_handler->fp_kseq->name.s);
    }
    if (fp_handler->fp_kseq->comment.l > 0) {
        header += std::string(" ") + std::string(fp_handler->fp_kseq->comment.s);
    }

    seq->data().insert(seq->data().end(), (int8_t *) fp_handler->fp_kseq->seq.s, (int8_t *) (fp_handler->fp_kseq->seq.s + fp_handler->fp_kseq->seq.l));
    seq->header(header);
    seq->id(id);
    seq->abs_id(abs_id);

    if (convert_to_uppercase) {
        seq->ToUppercase();
    }

    return seq;
}

mindex::SequencePtr SequenceDeserializer::DeserializeSequenceFromGFA1(SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase) {
    int64_t id = -1;
    int64_t abs_id = -1;

    mindex::SequencePtr seq = nullptr;

    std::string line;
    while (!raptor::ReadGZLine(fp_handler->fp_gzip, line)) {
        if (line.size() == 0) {
            continue;
        }
        if (line[0] == 'S') {
            std::istringstream ss(line);
            std::string keyword, header, seq_data;
            ss >> keyword >> header >> seq_data;

            if (seq_data.size() == 0 || seq_data == "*") {
                continue;
            }

            seq = mindex::createSequence();
            seq->data().insert(seq->data().end(), (int8_t *) &seq_data[0], (int8_t *) (&seq_data[0] + seq_data.size()));
            seq->header(header);
            seq->id(id);
            seq->abs_id(abs_id);

            if (convert_to_uppercase) {
                seq->ToUppercase();
            }

            // We found a valid sequence. Yield.
            break;
        }
    }

    return seq;
}

mindex::SequencePtr SequenceDeserializer::DeserializeSequenceFromGFA2(SequenceFileHandlersPtr& fp_handler, bool convert_to_uppercase) {
    int64_t id = -1;
    int64_t abs_id = -1;

    mindex::SequencePtr seq = nullptr;

    std::string line;

    while (raptor::ReadGZLine(fp_handler->fp_gzip, line)) {
        if (line.size() == 0) {
            continue;
        }
        if (line[0] == 'S') {
            std::istringstream ss(line);
            std::string keyword, header, seq_data;
            int64_t seq_len = 0;
            ss >> keyword >> header >> seq_len >> seq_data;

            // fprintf (stderr, "Line: '%s'\n", line.c_str());
            // fprintf (stderr, "keyword = %s, header = %s, seq_len = %d\n", keyword.c_str(), header.c_str(), seq_len);

            if (seq_data.size() == 0 || seq_data == "*") {
                continue;
            }

            seq = mindex::createSequence();
            seq->data().insert(seq->data().end(), (int8_t *) &seq_data[0], (int8_t *) (&seq_data[0] + seq_data.size()));
            seq->header(header);
            seq->id(id);
            seq->abs_id(abs_id);

            if (convert_to_uppercase) {
                seq->ToUppercase();
            }

            // We found a valid sequence. Yield.
            break;
        }
    }

    return seq;
}

}
