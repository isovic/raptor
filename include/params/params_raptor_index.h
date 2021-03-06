/*
 * params_raptor_index.h
 *
 *  Created on: Dec 13, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_
#define SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <index/params_index.h>
#include <sequences/sequence_file_enums.h>

namespace raptor {

class ParamsRaptorIndex {
   public:
    friend std::shared_ptr<raptor::ParamsRaptorIndex> createParamsRaptorIndex();
    ~ParamsRaptorIndex() = default;

    /*
      Command-line related options.
    */
    std::string subprogram = "";
    std::string command_line = "";

    std::shared_ptr<mindex::ParamsIndex> index_params = nullptr;
    int64_t verbose_level = 0;
    std::vector<std::string> ref_paths;
    // SequenceFormat ref_fmt;
    std::string out_path = "";
    double batch_size = 400;
    int64_t num_threads = 1;
    bool keep_lowercase = false;

   private:
    ParamsRaptorIndex() {
        index_params = mindex::createParamsIndex();
    }
    ParamsRaptorIndex(const ParamsRaptorIndex&) = delete;
    ParamsRaptorIndex& operator=(const ParamsRaptorIndex&) = delete;
};

inline std::shared_ptr<raptor::ParamsRaptorIndex> createParamsRaptorIndex() {
    return std::shared_ptr<raptor::ParamsRaptorIndex>(new raptor::ParamsRaptorIndex);
}

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_ */
