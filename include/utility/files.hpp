#ifndef SRC_UTILITY_FILES_H_
#define SRC_UTILITY_FILES_H_

#include <unistd.h>
#include <zlib.h>
#include <algorithm>
#include <string>

namespace raptor {

inline bool FileExists(const std::string& fname) {
    /* from
    * http://stackoverflow.com/questions/230062/whats-the-best-way-to-check-if-a-file-exists-in-c-cross-platform
    */
    if (access(fname.c_str(), F_OK) != -1) {
        return true;
    }
    return false;
}

inline bool ReadGZLine(gzFile gzip_fp, std::string &ret) {
    const int32_t BUFF_SIZE = 4096;
    int8_t buff[BUFF_SIZE + 1];

    int32_t read_bytes = 0;
    bool ln_end_found = false;
    int32_t ln_end_pos = 0;

    ret = std::string();

    // bool is_eof = false;
    while ((read_bytes = gzread(gzip_fp, &buff[0], BUFF_SIZE)) > 0) {
        buff[read_bytes] = '\0';
        for (int32_t i=0; i<read_bytes; i++) {
            if (buff[i] == '\n') {
                buff[i] = '\0';
                gzseek(gzip_fp, -(read_bytes - i - 1), SEEK_CUR);
                ln_end_found = true;
                ln_end_pos = i;
                break;
            }
        }
        std::string buff_str((char *) &buff[0]);
        ret += buff_str;

        if (ln_end_found) break;
    }

    int32_t ret_len = ret.size();

    // This allows for empty lines to be output as well (but not if there is an EOF, i.e. read_bytes == 0).
    if (ret_len > 0 || read_bytes > 0) return true;

    return false;
}

inline std::string GetFileExt(const std::string& path) {
    int32_t pos = path.find_last_of(".");
    std::string ext = path.substr(pos + 1);
    if (ext == "gz") {
        ext = path.substr(path.find_last_of(".", (pos - 1)) + 1);
    }
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

}

#endif
