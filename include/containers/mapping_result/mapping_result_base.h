/*
 * mapping_result_base.h
 *
 *  Created on: Nov 02, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_MAPPING_RESULT_BASE_H_
#define SRC_CONTAINERS_MAPPING_RESULT_BASE_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <containers/region/region_base.h>
#include <raptor/index_factory.h>

namespace raptor {

enum class MapperReturnValueBase {
    OK,
    ParamsNotSet,
    IndexVectorEmpty,
    IndexIsNull,
    QlenIsZero,
    QseqTooShort,
    NotRunYet,
    Failed
};

class MappingResultBase {
public:
    virtual ~MappingResultBase() {
    }

    virtual std::vector<std::shared_ptr<raptor::RegionBase>> CollectRegions(bool one_hit_per_target) const = 0;
    virtual int64_t QueryId() const = 0;
    virtual int64_t QueryLen() const = 0;
    virtual std::string QueryHeader() const = 0;
    virtual const mindex::IndexPtr Index() const = 0;
    virtual MapperReturnValueBase ReturnValue() const = 0;
    virtual const std::unordered_map<std::string, double>& Timings() const = 0;
};

}

#endif
