/*
 * target_hits.hpp
 *
 *  Created on: Sep 03, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_TARGET_HITS_H_
#define SRC_TARGET_HITS_H_

#include <stdint.h>

#include <containers/region/region_mapped.h>
#include <containers/mapping_env.h>

namespace raptor {

template<class T>
class TargetHits;

// template<class T>
// std::shared_ptr<raptor::TargetHits<T>> createTargetHits();

// template<class T>
// std::shared_ptr<raptor::TargetHits<T>> createTargetHits(std::shared_ptr<raptor::MappingEnv> env);

template<class T>
class TargetHits {
 public:
    // friend std::shared_ptr<raptor::TargetHits<T>> createTargetHits();

    // friend std::shared_ptr<raptor::TargetHits<T>> createTargetHits(std::shared_ptr<raptor::MappingEnv> env);
    TargetHits() :  env_(nullptr), hits_(),
                    cov_bases_t_(0), cov_bases_q_(0), num_seeds_(0),
                    score_(0) { }

    TargetHits(std::shared_ptr<raptor::MappingEnv> _env) :
                env_(_env), hits_(),
                cov_bases_t_(0), cov_bases_q_(0), num_seeds_(0),
                score_(0) { }

    TargetHits(std::shared_ptr<raptor::MappingEnv> _env, const std::vector<T>& _hits) :
                env_(_env), hits_(_hits),
                cov_bases_t_(0), cov_bases_q_(0), num_seeds_(0),
                score_(0) { }

    // TargetHits(const TargetHits<T>& op) :
    //             env_(op.env()), hits_(op.hits()), cov_bases_t_(op.cov_bases_t()),
    //             cov_bases_q_(op.cov_bases_q()), num_seeds_(op.num_seeds()) { }

    // TargetHits(const std::shared_ptr<TargetHits<T>>& op) :
    //             env_(op->env()), hits_(op->hits()), cov_bases_t_(op->cov_bases_t()),
    //             cov_bases_q_(op->cov_bases_q()), num_seeds_(op->num_seeds()) { }

    ~TargetHits() { }

    void AppendHits(const std::vector<T>& new_hits) {
        hits_.insert(hits_.end(), new_hits.begin(), new_hits.end());
    }

    std::string Verbose() const {
        std::ostringstream oss;
        oss << "Target hits:" << std::endl;
        oss << "  Env: " << env_->Verbose() << std::endl;
        oss << "  Hits:" << std::endl;
        for (size_t i = 0; i < hits().size(); i++) {
            oss << "[" << i << "] " << hits()[i].Verbose() << std::endl;
        }
        return oss.str();
    }

    std::string VerbosePointers() const {
        std::ostringstream oss;
        oss << "Target hits:" << std::endl;
        oss << "  Env: " << env_->Verbose() << std::endl;
        oss << "  Hits:" << std::endl;
        for (size_t i = 0; i < hits().size(); i++) {
            oss << "[" << i << "] " << hits()[i]->Verbose() << std::endl;
        }
        return oss.str();
    }

    // Getters.
    const std::shared_ptr<raptor::MappingEnv> env() const {
        return env_;
    }

    const std::vector<T>& hits() const {
        return hits_;
    }

    std::vector<T>& hits() {
        return hits_;
    }

    int32_t cov_bases_t() const {
        return cov_bases_t_;
    }

    int32_t cov_bases_q() const {
        return cov_bases_q_;
    }

    int32_t num_seeds() const {
        return num_seeds_;
    }

    int32_t score() const {
        return score_;
    }

    // Setters.
    void env(std::shared_ptr<raptor::MappingEnv> env) {
        env_ = env;
    }

    void hits(const std::vector<T>& _hits) {
        hits_ = _hits;
    }

    void cov_bases_t(int32_t _cov_bases_t) {
        cov_bases_t_ = _cov_bases_t;
    }

    void cov_bases_q(int32_t _cov_bases_q) {
        cov_bases_q_ = _cov_bases_q;
    }

    void num_seeds(int32_t _num_seeds) {
        num_seeds_ = _num_seeds;
    }

    void score(int32_t _score) {
        score_ = _score;
    }

private:
    std::shared_ptr<raptor::MappingEnv> env_;
    std::vector<T> hits_;
    int32_t cov_bases_t_;
    int32_t cov_bases_q_;
    int32_t num_seeds_;
    int32_t score_;            // A score to rank different chains by. E.g. if the hits were chained, than this can be the DP score of chaining.

};

// template<class T>
// std::shared_ptr<raptor::TargetHits<T>> createTargetHits() {
//     return std::shared_ptr<raptor::TargetHits<T>>(new raptor::TargetHits<T>());
// }

// template<class T>
// std::shared_ptr<raptor::TargetHits<T>> createTargetHits(std::shared_ptr<raptor::MappingEnv> env) {
//     return std::shared_ptr<raptor::TargetHits<T>>(new raptor::TargetHits<T>(env));
// }

}

#endif /* RANGE_H_ */
