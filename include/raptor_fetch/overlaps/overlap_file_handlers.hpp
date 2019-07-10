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
#include <cstdio>
#include <log/log_tools.h>

#include <iostream>

namespace raptor {

class OverlapFileHandlers;

typedef std::unique_ptr<raptor::OverlapFileHandlers> OverlapFileHandlersPtr;

std::unique_ptr<raptor::OverlapFileHandlers> createOverlapFileHandlers(const std::string& in_path);

class OverlapFileHandlers {
public:
    friend std::unique_ptr<raptor::OverlapFileHandlers> createOverlapFileHandlers(const std::string& in_path);

    ~OverlapFileHandlers() {
        fclose(fp);
        in_path = "";
        fp = NULL;
    }

    bool ReadLine(std::string& line) {
        if (fp == NULL) {
            return false;
        }
        const size_t BUFF_SIZE = 4096;
        char line_c[BUFF_SIZE + 1];
        line = "";
        size_t num_read = 0;
        // while((num_read = fgets(line_c, BUFF_SIZE, fp))) {
        while((num_read = fread(&line_c[0], sizeof(char), BUFF_SIZE, fp))) {
            // line = std::string(line_c);
            line_c[num_read] = '\0';
            size_t len = num_read;
            bool found = false;
            for (size_t i = 0; i < num_read; ++i) {
                if (line_c[i] == '\n' || line_c[i] == '\0') {
                    line_c[i] = '\0';
                    size_t curr_pos = ftell(fp);
                    fseek(fp, curr_pos - (num_read - i) + 1, SEEK_SET);
                    found = true;
                    len = i;
                    break;
                }
            }
            line += std::string(line_c, len);
            if (found == true || num_read < BUFF_SIZE) {
                return true;
            }
        }
        return false;
    }

    std::string in_path;
    FILE* fp;

private:
    OverlapFileHandlers(const std::string& _in_path)
                : in_path(_in_path), fp(NULL) {

        fp = fopen(_in_path.c_str(), "r");

        if (fp == NULL) {
            FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", _in_path.c_str());
        }

        in_path = _in_path;
    }

    OverlapFileHandlers(const OverlapFileHandlers&) = delete;
    OverlapFileHandlers& operator=(const OverlapFileHandlers&) = delete;

};

inline std::unique_ptr<raptor::OverlapFileHandlers> createOverlapFileHandlers(const std::string& in_path) {
    return std::unique_ptr<raptor::OverlapFileHandlers>(new raptor::OverlapFileHandlers(in_path));
}

}

#endif
