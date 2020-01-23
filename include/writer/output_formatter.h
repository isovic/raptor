/*
 * output_formatter.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_WRITER_OUTPUT_FORMATTER_H_
#define SRC_WRITER_OUTPUT_FORMATTER_H_

#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <containers/region/region_base.h>
#include <raptor/yield_index.h>
#include <types/typedefs.h>
#include <utility/revcmp.hpp>
#include <aligner/aligner_util.hpp>
#include <sequences/sequence.h>
#include <utility/stringutil.h>

namespace raptor {

class OutputFormatter {
public:
    static std::string TimingMapToString(const std::unordered_map<std::string, double>& timings);

    static std::string UnmappedSAM(const mindex::SequencePtr& qseq, bool write_custom_tags);

    static std::string ToSAM(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                             const std::shared_ptr<raptor::RegionBase> mapping,
                             bool write_custom_tags,
                             const std::string& timings);
    static std::string ToPAF(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                             const std::shared_ptr<raptor::RegionBase> mapping,
                             bool write_custom_tags,
                             const std::string& timings);
    static std::string ToMHAP(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                              const std::shared_ptr<raptor::RegionBase> mapping);
    static std::string ToM4(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                            const std::shared_ptr<raptor::RegionBase> mapping);
    static std::string ToGFA2Edge(const mindex::IndexPtr index, const mindex::SequencePtr& qseq,
                                const std::shared_ptr<raptor::RegionBase> mapping,
                                bool write_custom_tags,
                                const std::string& timings);

private:
    static double CalcIdentity_(
                int32_t q_start, int32_t q_end,
                int32_t t_start, int32_t t_end,
                int32_t cov_bases_q, int32_t cov_bases_t,
                int32_t edit_dist, int32_t match_bases);

};

} /* namespace raptor */

#endif /* SRC_WRITER_OUTPUT_FORMATTER_H_ */
