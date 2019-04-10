/*
 * sequence_file_handlers.h
 *
 *  Created on: Jun 15, 2018
 *      Author: Ivan Sovic
 */

#ifndef SRC_SEQUENCES_SEQUENCE_FILE_HANDLERS_H_
#define SRC_SEQUENCES_SEQUENCE_FILE_HANDLERS_H_

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <zlib.h>
#include <sequences/kseq.h>
#include <log/log_tools.h>

namespace mindex {

class SequenceFileHandlers;

typedef std::unique_ptr<mindex::SequenceFileHandlers> SequenceFileHandlersPtr;

// std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers();
std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers(const std::string& _in_path);

KSEQ_INIT(gzFile, gzread)

class SequenceFileHandlers {
public:
    // friend std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers();
    friend std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers(const std::string& _in_path);

    bool Open(const std::string& _in_path) {
        in_path_ = _in_path;

        fp_gzip = gzopen(_in_path.c_str(), "r");

        if (fp_gzip == NULL) {
            FATAL_REPORT(ERR_OPENING_FILE, "File path: '%s'.", _in_path.c_str());
            return false;
        }

        fp_kseq = kseq_init(fp_gzip);
        return true;
    }

    ~SequenceFileHandlers() {
        if (fp_kseq != NULL) {
            // FATAL_REPORT(ERR_CLOSING_FILE, "Offending variable: fp_kseq.");
            kseq_destroy(fp_kseq);
            fp_kseq = NULL;
        }

        if (fp_gzip != NULL) {
            gzclose(fp_gzip);
            fp_gzip = NULL;
        }

        in_path_ = "";
    }

    bool Seek(int64_t abs_pos) {
        if (IsOpen_() == false) {
            return false;
        }
        gzseek(fp_gzip, abs_pos, SEEK_SET);
        kseq_destroy(fp_kseq);
        fp_kseq = kseq_init(fp_gzip);
        return true;
    }

    int64_t Tell() const {
        if (IsOpen_() == false) {
            return -1;
        }
        return gztell(fp_gzip);
    }

    const std::string& in_path() const {
        return in_path_;
    }

    kseq_t *fp_kseq;
    gzFile fp_gzip;

private:
    SequenceFileHandlers() : fp_kseq(NULL), fp_gzip(NULL), in_path_() { }

    // SequenceFileHandlers(const std::string& _in_path)
    //             : in_path(_in_path), fp_kseq(NULL), fp_gzip(NULL) {

    // }

    SequenceFileHandlers(const SequenceFileHandlers&) = delete;
    SequenceFileHandlers& operator=(const SequenceFileHandlers&) = delete;

    std::string in_path_;

    bool IsOpen_() const {
        if (in_path_.empty() || fp_kseq == nullptr || fp_gzip == nullptr) {
            return false;
        }
        return true;
    }

};

// inline std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers() {
//     return std::unique_ptr<mindex::SequenceFileHandlers>(new mindex::SequenceFileHandlers());
// }

inline std::unique_ptr<mindex::SequenceFileHandlers> createSequenceFileHandlers(const std::string& _in_path) {
    auto ret = std::unique_ptr<mindex::SequenceFileHandlers>(new mindex::SequenceFileHandlers());
    ret->Open(_in_path);
    return std::move(ret);
}

}

#endif
