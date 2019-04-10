/*
 * path_alignment.h
 *
 *  Created on: Jan 03, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTANIERS_PATH_ALIGNMENT_H_
#define SRC_CONTANIERS_PATH_ALIGNMENT_H_

#include <stdint.h>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <aligner/alignment_result.h>
#include <containers/region/region_aligned.h>
#include <graph/anchor_graph.h>
#include <graph/local_path.h>

namespace raptor {

class PathAlignment;

typedef std::shared_ptr<raptor::PathAlignment> PathAlignmentPtr;

std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path);
std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path, int64_t _path_score,
                                                       const std::vector<std::shared_ptr<raptor::RegionAligned>>& _aln,
                                                       std::shared_ptr<raptor::AlignmentResult> _entire_alignment);

/*
 * Holds all alignment information for a single LocalPath.
 * This includes: a pointer to the path in question, the total alignment
 * score (needed for comparison of two paths), and a vector of
 * aligned regions of the paths (each of which can be to a different target).
*/
class PathAlignment {
 public:
    friend std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path);
    friend std::shared_ptr<raptor::PathAlignment> createPathAlignment(const std::shared_ptr<raptor::LocalPath> _path, int64_t _path_score,
                                                                  const std::vector<std::shared_ptr<raptor::RegionAligned>>& _aln,
                                                                  std::shared_ptr<raptor::AlignmentResult> _entire_alignment);

    ~PathAlignment() { }

    /*
     * Getters.
    */
    const std::shared_ptr<raptor::LocalPath> path() const {
        return path_;
    }
    const std::vector<std::shared_ptr<raptor::RegionAligned>>& alns() const {
        return alns_;
    }
    int64_t path_score() const {
        return path_score_;
    }
    const std::shared_ptr<raptor::AlignmentResult> entire_alignment() const {
        return entire_alignment_;
    }
    const std::unordered_map<std::string, double>& timings() const {
        return timings_;
    }

    /*
     * Setters.
    */
    // void path(std::shared_ptr<raptor::LocalPath> _path) {
    //     path_ = _path;
    // }
    void alns(const std::vector<std::shared_ptr<raptor::RegionAligned>>& _alns) {
        alns_ = _alns;
    }
    void path_score(int64_t _path_score) {
        path_score_ = _path_score;
    }
    void entire_alignment(std::shared_ptr<raptor::AlignmentResult> _entire_alignment) {
        entire_alignment_ = _entire_alignment;
    }
    void add_timing(const std::string& name, double val) {
        timings_[name] = val;
    }

 private:
    PathAlignment(std::shared_ptr<raptor::LocalPath> _path);

    PathAlignment(const std::shared_ptr<raptor::LocalPath> _path, int64_t _path_score,
                  const std::vector<std::shared_ptr<raptor::RegionAligned>>& _alns,
                  std::shared_ptr<raptor::AlignmentResult> _entire_alignment);

    const std::shared_ptr<raptor::LocalPath> path_;
    int64_t path_score_;
    std::vector<std::shared_ptr<raptor::RegionAligned>> alns_;
    std::shared_ptr<raptor::AlignmentResult> entire_alignment_;
    std::unordered_map<std::string, double> timings_;
};

}

#endif
