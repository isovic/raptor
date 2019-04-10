/*
 * alignment_position.h
 *
 *  Created on: Dec 1, 2017
 *      Author: isovic
 */

#ifndef SRC_ALIGNER_ALIGNMENT_POSITION_H_
#define SRC_ALIGNER_ALIGNMENT_POSITION_H_

#include <cstdint>

namespace raptor {

class AlignmentPosition {
   public:
    AlignmentPosition() : qstart(0), qend(0), tstart(0), tend(0) {}
    AlignmentPosition(int32_t _qstart, int32_t _qend, int32_t _tstart, int32_t _tend)
        : qstart(_qstart), qend(_qend), tstart(_tstart), tend(_tend) {}
    AlignmentPosition(const AlignmentPosition& op)
        : AlignmentPosition(op.qstart, op.qend, op.tstart, op.tend) {}

    int32_t qstart, qend;  // Query and target alignment start and end positions. End position
    int32_t tstart, tend;  // is inclusive (the position of the last base).
};

}  // namespace raptor

#endif
