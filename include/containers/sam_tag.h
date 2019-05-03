/*
 * mapping_base.h
 *
 *  Created on: Dec 3, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_CONTAINERS_SAM_TAG_H_
#define SRC_CONTAINERS_SAM_TAG_H_

#include <string>

namespace raptor {

class SamTag {
    public:
    SamTag() : name(), type(), val() {
    }

    SamTag(const std::string& _name, const std::string& _type, const std::string& _val)
                    : name(_name), type(_type), val(_val) {
    }

    ~SamTag() = default;

    // The following three fields should be consistent with SAM/BAM conventions.
    // They can simply be joined with a ':', and they would automatically be
    // compliant with the SAM/BAM formats.

    // The SAM/BAM compliant naming. It should include the data type too. E.g. cg:z
    std::string name;
    // The SAM/BAM compliant naming of a tag type. E.g. i for integer.
    std::string type;
    // Values packed in a string, regardless of the underlying data type, for simplicity.
    std::string val;

    std::string FormatAsSAM() const {
        return name + ":" + type + ":" + val;
    }
};

}

#endif
