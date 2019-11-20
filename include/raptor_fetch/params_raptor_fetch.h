/*
 * params_raptor_fetch.h
 *
 *  Created on: Dec 22, 2017
 *      Author: Ivan Sovic
 *
 * Raptor-fetch is intended to perform random access to the Raptor DB.
 * reshape the input data.
 */

#ifndef SRC_RAPTOR_PARAMS_RAPTOR_FETCH_H_
#define SRC_RAPTOR_PARAMS_RAPTOR_FETCH_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace raptor {

class ParamsRaptorFetch {
   public:
    friend std::shared_ptr<raptor::ParamsRaptorFetch> createParamsRaptorFetch();
    ~ParamsRaptorFetch() = default;

    int64_t verbose_level = 0;
    std::string command_line;
    int64_t debug_qiq = -1;
    std::string debug_qname;

    std::string rdb_path;
    std::vector<std::string> in_paths;
    std::string out_prefix;
    int32_t min_cov = 0;
    int64_t min_len = 0;
    double min_score = 0.0;
    double min_idt = 0.0;
    bool use_id_for_output = false;
    std::string job_str = "fetch";

   private:
    ParamsRaptorFetch() = default;
    ParamsRaptorFetch(const ParamsRaptorFetch&) = delete;
    ParamsRaptorFetch& operator=(const ParamsRaptorFetch&) = delete;
};

inline std::shared_ptr<raptor::ParamsRaptorFetch> createParamsRaptorFetch() {
    return std::shared_ptr<raptor::ParamsRaptorFetch>(new raptor::ParamsRaptorFetch);
}

}  // namespace raptor

#endif /* SRC_RAPTOR_PARAMS_RAPTOR_INDEX_H_ */
