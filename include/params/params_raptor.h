/*
 * params_raptor.h
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_PARAMS_RAPTOR_H_
#define SRC_PARAMS_RAPTOR_H_

#include <cstdint>
#include <memory>
#include <string>
#include <sequences/sequence_file_enums.h>
#include <graph/segment_graph_parser_enums.h>
#include <writer/output_formatter_enums.h>
#include <index/params_index.h>
#include <index/index_types.h>
#include <params/params_mapper.h>
#include <params/params_aligner.h>

namespace raptor {

class ParamsRaptor {
  public:
    friend std::shared_ptr<raptor::ParamsRaptor> createParamsRaptor();
    ~ParamsRaptor() = default;

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
    std::string graph_path = "";
    std::string out_path = "";

    mindex::SequenceFormat ref_fmt = mindex::SequenceFormat::Unknown;
    mindex::SequenceFormat infmt = mindex::SequenceFormat::Unknown;
    raptor::GraphFormat graph_fmt = GraphFormat::Unknown;
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

    bool rebuild_index = false;
    bool auto_rebuild_index = false;
    bool index_on_the_fly = false;
    bool calc_only_index = false;
    bool add_symmetric_arcs = false;

    /*
      Mapping options.
    */
    std::shared_ptr<raptor::ParamsMapper> mapper_params = raptor::createParamsMapper();

    /*
      Aligner options.
    */
    std::shared_ptr<raptor::ParamsAligner> aligner_params = raptor::createParamsAligner();

    /*
     * Filtering options.
    */
    double min_identity = 65.0;
    double max_evalue = 1e0;
    int32_t min_mapq = 0;
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
    int64_t verbose_output = 0;

    /*
     These need to be classified in categories:
    */
    int64_t bestn = 0;
    double bestn_threshold = 0.01;
    int64_t min_map_len = 0;
    bool do_align = false;
    bool strict_format = false;
    bool do_diff = false;

    std::string composite;

  private:
    ParamsRaptor() = default;
    ParamsRaptor(const ParamsRaptor&) = delete;
    ParamsRaptor& operator=(const ParamsRaptor&) = delete;
};

inline std::shared_ptr<raptor::ParamsRaptor> createParamsRaptor() {
    return std::shared_ptr<raptor::ParamsRaptor>(new raptor::ParamsRaptor());
}

} /* namespace raptor */

#endif /* SRC_PARAMS_RAPTOR_H_ */
