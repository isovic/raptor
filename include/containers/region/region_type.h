/*
 * mapping_base.h
 *
 *  Created on: Jul 19, 2019
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_REGION_TYPE_H_
#define SRC_CONTAINERS_REGION_TYPE_H_

namespace raptor {

enum class RegionType {
    Primary,
    SupplementaryPrimary,           // Extends the primary alignment.
    Secondary,
    SupplementarySecondary,          // Extends any non primary alignment.
    Undefined
};

}

#endif
