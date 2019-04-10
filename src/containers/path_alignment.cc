/*
 * path_alignment.cc
 *
 *  Created on: Jan 03, 2018
 *      Author: Ivan Sovic
 */

#include <containers/path_alignment.h>

namespace raptor {

std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path) {
    return std::shared_ptr<raptor::PathAlignment>(new raptor::PathAlignment(_path));
}

std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path,
                                                       int64_t _path_score,
                                                       const std::vector<std::shared_ptr<raptor::RegionAligned>>& _alns,
                                                       std::shared_ptr<raptor::AlignmentResult> _entire_alignment) {
    return std::shared_ptr<raptor::PathAlignment>(new raptor::PathAlignment(_path, _path_score, _alns, _entire_alignment));
}

PathAlignment::PathAlignment(const std::shared_ptr<raptor::LocalPath> _path)
    :   path_(_path),
        path_score_(0),
        alns_(),
        entire_alignment_(nullptr),
        timings_() {
}

PathAlignment::PathAlignment(const std::shared_ptr<raptor::LocalPath> _path,
                             int64_t _path_score,
                             const std::vector<std::shared_ptr<raptor::RegionAligned>>& _alns,
                             std::shared_ptr<raptor::AlignmentResult> _entire_alignment)
    :   path_(_path),
        path_score_(_path_score),
        alns_(_alns),
        entire_alignment_(_entire_alignment),
        timings_() {

}

}
