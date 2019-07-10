/*
 * overlap_file_handlers.h
 *
 *  Created on: Jan 06, 2019
 *      Author: Ivan Sovic
 */

#ifndef INCLUDE_OVERLAP_FILE_HANDLERS_H_
#define INCLUDE_OVERLAP_FILE_HANDLERS_H_

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <zlib.h>
#include <log/log_tools.h>

namespace raptor {

class OverlapFileHandlersGZ;

typedef std::unique_ptr<raptor::OverlapFileHandlersGZ> OverlapFileHandlersGZPtr;

std::unique_ptr<raptor::OverlapFileHandlersGZ> createOverlapFileHandlersGZ(const std::string& in_path);

class OverlapFileHandlersGZ {
public:
    friend std::unique_ptr<raptor::OverlapFileHandlersGZ> createOverlapFileHandlersGZ(const std::string& in_path);

    ~OverlapFileHandlersGZ() {
        gzclose(fp);
        in_path = "";
        fp = NULL;
    }

    std::string in_path;
    gzFile fp;

private:
    OverlapFileHandlersGZ(const std::string& _in_path)
                : in_path(_in_path), fp(NULL) {

        fp = gzopen(_in_path.c_str(), "r");

        if (fp == NULL) {
            FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", _in_path.c_str());
        }

        in_path = _in_path;
    }

    OverlapFileHandlersGZ(const OverlapFileHandlersGZ&) = delete;
    OverlapFileHandlersGZ& operator=(const OverlapFileHandlersGZ&) = delete;

};

inline std::unique_ptr<raptor::OverlapFileHandlersGZ> createOverlapFileHandlersGZ(const std::string& in_path) {
    return std::unique_ptr<raptor::OverlapFileHandlersGZ>(new raptor::OverlapFileHandlersGZ(in_path));
}

}

#endif
