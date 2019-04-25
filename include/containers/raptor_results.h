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
                  mapping_result(nullptr),
                  graph_mapping_result(nullptr),
                  aln_result(nullptr) { }

  ~RaptorResults() { }

  int64_t q_id_in_batch;
  // Intermediate results.
  std::shared_ptr<raptor::LinearMappingResult> mapping_result;
  std::shared_ptr<raptor::GraphMappingResult> graph_mapping_result;
  std::shared_ptr<raptor::AlignedMappingResult> aln_result;
  // Final results are extracted as linear regions through the graph.
  std::vector<std::shared_ptr<raptor::RegionBase>> regions;
};

} /* namespace raptor */

#endif /* SRC_CONTAINERS_RAPTOR_RESULTS_H_ */
