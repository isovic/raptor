/*
 * linear_mapping_result.h
 *
 *  Created on: Sep 03, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_MAPPING_RESULT_H_
#define SRC_CONTAINERS_MAPPING_RESULT_H_

#include <stdint.h>
#include <memory>
#include <string>
#include <vector>
#include <containers/target_hits.hpp>
#include <raptor/yield_index.h>
#include <types/typedefs.h>
#include <containers/mapping_result/mapping_result_base.h>

namespace raptor {

class LinearMappingResult;

/*
 * Using a shared_ptr instead of unique_ptr because this structure might
 * get juggled around and moved from vector to vector.
 * Setting it as a shared_ptr reduced the wrong usage at the expense of
 * a few bytes of space.
*/
std::shared_ptr<raptor::LinearMappingResult> createMappingResult(
                        int64_t qseq_len,
                        int64_t qseq_id,
                        std::string
                        qseq_header,
                        mindex::IndexPtr index);

class LinearMappingResult : public raptor::MappingResultBase {
public:
    friend std::shared_ptr<raptor::LinearMappingResult> createMappingResult(
                        int64_t qseq_len,
                        int64_t qseq_id,
                        std::string
                        qseq_header,
                        mindex::IndexPtr index);

    // Interface implementation.
    std::vector<std::shared_ptr<raptor::RegionBase>> CollectRegions(bool one_hit_per_target) const;
    int64_t QueryId() const;
    int64_t QueryLen() const;
    std::string QueryHeader() const;
    const mindex::IndexPtr Index() const;
    MapperReturnValueBase ReturnValue() const;
    const std::unordered_map<std::string, double>& Timings() const;

    std::string Verbose() const;
    std::string WriteAsCSV(const char separator) const;

    // Custom getters / setters.
    int64_t qseq_len() const {
        return qseq_len_;
    }

    int64_t qseq_id() const {
        return qseq_id_;
    }

    void target_anchors(const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& new_anchors) {
        target_anchors_ = new_anchors;
    }

    void target_hits(const std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>>& new_target_hits) {
        target_hits_ = new_target_hits;
    }

    const std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors() const {
        return target_anchors_;
    }

    const std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>>& target_hits() const {
        return target_hits_;
    }

    void ClearTargetHits() {
        target_hits_.clear();
    }

    void return_value(MapperReturnValueBase _return_value) {
        return_value_ = _return_value;
    }

    MapperReturnValueBase return_value() const {
        return return_value_;
    }

    std::unordered_map<std::string, double>& timings() {
        return timings_;
    }

    const std::unordered_map<std::string, double>& timings() const {
        return timings_;
    }

    void timings(const std::unordered_map<std::string, double>& _timings) {
        timings_ = _timings;
    }

    std::vector<std::shared_ptr<raptor::TargetAnchorType>>& target_anchors() {
        return target_anchors_;
    }

    std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>>& target_hits() {
        return target_hits_;
    }

private:
    LinearMappingResult(
            int64_t qseq_len,
            int64_t qseq_id,
            std::string
            qseq_header,
            mindex::IndexPtr index);
    LinearMappingResult(const LinearMappingResult&) = delete;
    LinearMappingResult& operator=(const LinearMappingResult&) = delete;

    int64_t qseq_id_;
    int64_t qseq_len_;
    std::string qseq_header_;
    mindex::IndexPtr index_;
    raptor::MapperReturnValueBase return_value_;

    // Anchors summarize a colinear chain of target hits.
    std::vector<std::shared_ptr<raptor::TargetAnchorType>> target_anchors_;
    // Corresponding target hits (seeds) that make up the anchors above.
    std::vector<std::shared_ptr<raptor::TargetHits<mindex::SeedHitPacked>>> target_hits_;
    // Time measurements, useful for debugging.
    std::unordered_map<std::string, double> timings_;
};

}

#endif
