/*
 * raptor_results.h
 *
 *  Created on: May 30, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_RAPTOR_RESULTS_H_
#define SRC_CONTAINERS_RAPTOR_RESULTS_H_

#include <containers/mapping_result/linear_mapping_result.h>
#include <containers/mapping_result/graph_mapping_result.h>
#include <containers/mapping_result/aligned_mapping_result.h>
#include <containers/path_alignment.h>

namespace raptor {

class RaptorResults {
 public:
  RaptorResults()
                : q_id_in_batch(-1),
                  regions(),
                  timings(),
                  mapq(0)
  { }

  ~RaptorResults() { }

  int64_t q_id_in_batch;
  std::vector<std::shared_ptr<raptor::RegionBase>> regions;
  std::unordered_map<std::string, double> timings;
  int32_t mapq;
};

} /* namespace raptor */

#endif /* SRC_CONTAINERS_RAPTOR_RESULTS_H_ */
