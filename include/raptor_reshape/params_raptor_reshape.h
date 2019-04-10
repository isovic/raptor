/*
 * params_raptor_reshape.h
 *
 *  Created on: Dec 22, 2017
 *      Author: Ivan Sovic
 *
 * Raptor-repack is intended to reshape the input data.
 * For example, the data can be split into blocks, and compressed,
 * where each block can be output separately.
 */

#ifndef SRC_RAPTOR_PARAMS_RAPTOR_REPACK_H_
#define SRC_RAPTOR_PARAMS_RAPTOR_REPACK_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <index/sequence_file_enums.h>

namespace raptor {

class ParamsRaptorReshape;

std::shared_ptr<raptor::ParamsRaptorReshape> createParamsRaptorReshape();

class ParamsRaptorReshape {
   public:
    friend std::shared_ptr<raptor::ParamsRaptorReshape> createParamsRaptorReshape();
    ~ParamsRaptorReshape() = default;

    /*
      Command-line related options.
    */
    std::string subprogram = "";
    std::string command_line = "";

    int64_t verbose_level;
    std::vector<std::string> in_paths;
    mindex::SequenceFormat in_fmt;
    mindex::SequenceFormat out_fmt;
    std::string out_prefix = "";
    double in_batch_size;
    bool keep_lowercase = false;
    bool rename_seqs = false;
    double block_size = 400.0;
    bool split_blocks = false;
    bool symlink_files = false;

   private:
    ParamsRaptorReshape();
    ParamsRaptorReshape(const ParamsRaptorReshape&) = delete;
    ParamsRaptorReshape& operator=(const ParamsRaptorReshape&) = delete;
};

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_ */
