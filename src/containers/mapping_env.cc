#include "containers/mapping_env.h"

namespace raptor {

std::shared_ptr<raptor::MappingEnv> createMappingEnv() {
    return std::shared_ptr<raptor::MappingEnv>(new raptor::MappingEnv());
}

std::shared_ptr<raptor::MappingEnv> createMappingEnv(int32_t _t_id, int64_t _index_t_start,
                                                 int32_t _t_len, bool _t_rev,
                                                 int32_t _q_id, int32_t _q_len,
                                                 bool _q_rev) {
    return std::shared_ptr<raptor::MappingEnv>(
                new raptor::MappingEnv(_t_id, _index_t_start, _t_len,
                _t_rev, _q_id, _q_len, _q_rev));
}

MappingEnv::MappingEnv() :  t_id(0), index_t_start(0), t_len(0), t_rev(0),
                            q_id(0), q_len(0), q_rev(0) {

}

MappingEnv::MappingEnv(int32_t _t_id, int64_t _index_t_start,
                        int32_t _t_len, bool _t_rev,
                        int32_t _q_id, int32_t _q_len, bool _q_rev) :
                        t_id(_t_id), index_t_start(_index_t_start),
                        t_len(_t_len), t_rev(_t_rev), q_id(_q_id),
                        q_len(_q_len), q_rev(_q_rev) {
}

}
