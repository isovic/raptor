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
                        const std::unordered_map<std::string, double>& timings);

class RaptorResults {
  public:
    friend std::unique_ptr<raptor::RaptorResults> createRaptorResults(
                        int64_t q_id,
                        const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                        const std::unordered_map<std::string, double>& timings);

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

  private:
    RaptorResults()
                    : q_id_(-1),
                    regions_(),
                    timings_()
    { }
    RaptorResults(
                int64_t q_id,
                const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                const std::unordered_map<std::string, double>& timings)
                :   q_id_(q_id),
                    regions_(regions),
                    timings_(timings)
    { }

    int64_t q_id_;
    std::vector<std::shared_ptr<raptor::RegionBase>> regions_;
    std::unordered_map<std::string, double> timings_;
};

inline std::unique_ptr<raptor::RaptorResults> createRaptorResults(
                        int64_t q_id,
                        const std::vector<std::shared_ptr<raptor::RegionBase>>& regions,
                        const std::unordered_map<std::string, double>& timings) {
    return std::unique_ptr<raptor::RaptorResults>(new raptor::RaptorResults(q_id, regions, timings));

}

} /* namespace raptor */

#endif /* SRC_CONTAINERS_RAPTOR_RESULTS_H_ */
