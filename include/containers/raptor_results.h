/*
 * raptor_results.h
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_RAPTOR_RESULTS_H_
#define SRC_CONTAINERS_RAPTOR_RESULTS_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <containers/region/region_base.h>

namespace raptor {

class RaptorResults;

std::unique_ptr<raptor::RaptorResults> createRaptorResults(
                        int64_t q_id,
                        const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                        const std::unordered_map<std::string, double>& timings,
                        int32_t mapq);

class RaptorResults {
  public:
    friend std::unique_ptr<raptor::RaptorResults> createRaptorResults(
                        int64_t q_id,
                        const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                        const std::unordered_map<std::string, double>& timings,
                        int32_t mapq);

    ~RaptorResults() { }

    int64_t q_id() const {
        return q_id_;
    }
    const std::vector<std::shared_ptr<raptor::RegionBase>>& regions() const {
        return regions_;
    }
    const std::unordered_map<std::string, double>& timings() const {
        return timings_;
    }
    int32_t mapq() const {
        return mapq_;
    }

  private:
    RaptorResults()
                    : q_id_(-1),
                    regions_(),
                    timings_(),
                    mapq_(0)
    { }
    RaptorResults(
                int64_t q_id,
                const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                const std::unordered_map<std::string, double>& timings,
                int32_t mapq)
                :   q_id_(q_id),
                    regions_(regions),
                    timings_(timings),
                    mapq_(mapq)
    { }

    int64_t q_id_;
    std::vector<std::shared_ptr<raptor::RegionBase>> regions_;
    std::unordered_map<std::string, double> timings_;
    int32_t mapq_;
};

inline std::unique_ptr<raptor::RaptorResults> createRaptorResults(
                        int64_t q_id,
                        const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                        const std::unordered_map<std::string, double>& timings,
                        int32_t mapq) {
    return std::unique_ptr<raptor::RaptorResults>(new raptor::RaptorResults(q_id, regions, timings, mapq));

}

} /* namespace raptor */

#endif /* SRC_CONTAINERS_RAPTOR_RESULTS_H_ */
