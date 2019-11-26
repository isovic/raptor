/*
 * params_raptor_sova.h
 *
 *  Created on: Nov 12, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PARAMS_RAPTOR_HOLY_H_
#define SRC_RAPTOR_PARAMS_RAPTOR_HOLY_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <sequences/sequence_file_enums.h>
#include <writer/output_formatter_enums.h>
#include <index/params_index.h>
#include <index/index_types.h>
#include <raptor_sova/params_sova_mapper.h>
#include <params/params_aligner.h>

namespace raptor {
namespace sova {

struct ParamsRaptorSova {
    /*
      Command-line related options.
    */
    std::string subprogram = "";
    std::string command_line = "";

    /*
      In/out options.
    */
    std::vector<std::string> ref_paths;
    std::vector<std::string> query_paths;
    std::string out_path = "";

    mindex::SequenceFormat ref_fmt = mindex::SequenceFormat::Unknown;
    mindex::SequenceFormat infmt = mindex::SequenceFormat::Unknown;
    raptor::OutputFormat outfmt = OutputFormat::Unknown;

    double batch_size = 100.0;
    mindex::BatchLoadType batch_type = mindex::BatchLoadType::MB;
    mindex::IndexType index_type = mindex::IndexType::Undefined;
    bool keep_lowercase = false;
    bool ref_and_reads_path_same = false;  // We need to know this for overlapping.
    int64_t rdb_block_ref = -1;
    int64_t rdb_block_query = -1;

    /*
      Index options.
    */
    std::shared_ptr<mindex::ParamsIndex> index_params = mindex::createParamsIndex();
    bool index_on_the_fly = false;
    bool calc_only_index = false;
    bool add_symmetric_arcs = false;

    /*
      Mapping options.
    */
    std::shared_ptr<raptor::sova::ParamsSovaMapper> mapper_params = raptor::sova::createParamsSovaMapper();

    /*
     * Filtering options.
    */
    bool one_hit_per_target = false;

    /*
      Misc.
    */
    int64_t num_threads = 1;
    int64_t verbose_level = 1;
    int64_t start_read = 0;
    int64_t num_reads_to_process = -1;

    /*
      Debug options.
    */
    int64_t debug_qid = -1;
    std::string debug_qname = "";

    /*
     These need to be classified in categories:
    */
    int64_t bestn = 0;
    bool strict_format = false;

    std::string composite;
};

inline std::shared_ptr<raptor::sova::ParamsRaptorSova> createParamsRaptorSova() {
    return std::shared_ptr<raptor::sova::ParamsRaptorSova>(new raptor::sova::ParamsRaptorSova());
}

} /* namespace sova */
}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_ */
