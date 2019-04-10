/*
 * mapping_env.h
 *
 *  Created on: Nov 27, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_MAPPING_ENV_H_
#define SRC_MAPPING_ENV_H_

#include <cstdint>
#include <sstream>
#include <string>
#include <memory>

namespace raptor {

class MappingEnv;

std::shared_ptr<raptor::MappingEnv> createMappingEnv();

std::shared_ptr<MappingEnv> createMappingEnv(int32_t _t_id, int64_t _index_t_start,
                                                 int32_t _t_len, bool _t_rev,
                                                 int32_t _q_id, int32_t _q_len,
                                                 bool _q_rev);
class MappingEnv {
public:
    friend std::shared_ptr<raptor::MappingEnv> createMappingEnv();

    friend std::shared_ptr<MappingEnv> createMappingEnv(int32_t _t_id, int64_t _index_t_start,
                                                            int32_t _t_len, bool _t_rev,
                                                            int32_t _q_id, int32_t _q_len,
                                                            bool _q_rev);

    std::string Verbose() const {
        std::ostringstream oss;
        oss << "t_id = " << t_id << ", index_t_start = "
            << index_t_start << ", t_len = "
            << t_len << ", t_rev = "
            << t_rev << ", q_id = "
            << q_id << ", q_len = "
            << q_len << ", q_rev = "
            << q_rev;
        return oss.str();
    }

    int32_t t_id;           // Target ID.
    int64_t index_t_start;  // Start position of the target/reference sequence in the index; e.g. index may contain all sequences concatenated in a linear array.
    int32_t t_len;          // Length of the target sequence.
    bool t_rev;             // Is target sequence reversed?
    int32_t q_id;           // ID of the query sequence.
    int32_t q_len;          // Length of the query sequence.
    bool q_rev;             // Is query sequence reversed? This should always be false, but we'll keep this for symmetry.

private:
    MappingEnv(const MappingEnv&) = delete;

    MappingEnv& operator=(const MappingEnv&) = delete;

    MappingEnv();

    MappingEnv(int32_t _t_id, int64_t _index_t_start,
               int32_t _t_len, bool _t_rev,
               int32_t _q_id, int32_t _q_len,
               bool _q_rev);
};

}

#endif
