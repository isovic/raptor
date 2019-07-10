#ifndef _IPA_INCLUDE_OVERLAP_PARSER_H_
#define _IPA_INCLUDE_OVERLAP_PARSER_H_

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <raptor_fetch/overlaps/overlap.h>
#include <utility/stringutil.h>

namespace raptor {

inline std::unique_ptr<raptor::Overlap> ParseOverlap(const std::string& line) {
    std::string a_name, b_name;
    float score = 0.0f, identity = 0.0f;
    bool a_rev = false, b_rev = false;
    coord_t a_start = 0, a_end = 0, a_len = 0;
    coord_t b_start = 0, b_end = 0, b_len = 0;
    coord_t clip_a_start = 0, clip_a_end = 0, clip_b_start = 0, clip_b_end = 0;
    std::string type;

    char a_name_buff[1000] = "\0", b_name_buff[1000] = "\0";
    char type_buff[1000] = "\0";
    int32_t a_rev_int = 0, b_rev_int = 0;

    // std::vector<std::string> tokens = raptor::TokenizeToWhitespaces(line);
    // a_name = tokens[0];
    // b_name = tokens[1];
    // score = std::stof(tokens[2]);
    // identity = std::stof(tokens[3]);
    // a_rev_int = std::stoi(tokens[4]);
    // a_start = std::stoi(tokens[5]);
    // a_end = std::stoi(tokens[6]);
    // a_len = std::stoi(tokens[7]);
    // b_rev_int = std::stoi(tokens[8]);
    // b_start = std::stoi(tokens[9]);
    // b_end = std::stoi(tokens[10]);
    // b_len = std::stoi(tokens[11]);
    // a_rev = a_rev_int != 0;
    // b_rev = b_rev_int != 0;
    // if (tokens.size() < 16) {
    //     clip_a_start = 0;
    //     clip_a_end = a_len;
    //     clip_b_start = 0;
    //     clip_b_end = b_len;
    // } else {
    //     clip_a_start = std::stoi(tokens[12]);
    //     clip_a_end = std::stoi(tokens[13]);
    //     clip_b_start = std::stoi(tokens[14]);
    //     clip_b_end = std::stoi(tokens[15]);
    // }
    // if (tokens.size() >= 17) {
    //     type = tokens[16];
    // }

    int32_t num_fields = sscanf(line.c_str(),
                            "%s %s %f %f "
                            "%d %d %d %d "
                            "%d %d %d %d "
                            "%d %d %d %d %s",
                            a_name_buff, b_name_buff, &score, &identity,
                            &a_rev_int, &a_start, &a_end, &a_len,
                            &b_rev_int, &b_start, &b_end, &b_len,
                            &clip_a_start, &clip_a_end,
                            &clip_b_start, &clip_b_end,
                            type_buff
                            );

    a_name = std::string(a_name_buff);
    b_name = std::string(b_name_buff);
    a_rev = a_rev_int != 0;
    b_rev = b_rev_int != 0;
    type = std::string(type_buff);

    if (num_fields < 16) {
        clip_a_start = 0;
        clip_a_end = a_len;
        clip_b_start = 0;
        clip_b_end = b_len;
    }

    // std::istringstream iss(line);
    // std::string param;
    // std::vector<std::string> params;
    // while (iss >> param) {
    //     params.emplace_back(param);
    // }

    // a_name = params[0];
    // b_name = params[1];
    // score = std::stof(params[2]);
    // identity = std::stof(params[3]);

    std::unique_ptr<raptor::Overlap> ret = raptor::createOverlap(
        a_name, b_name, score, identity,
        a_rev, a_start, a_end, a_len,
        b_rev, b_start, b_end, b_len,
        clip_a_start, clip_a_end, clip_b_start, clip_b_end, type);

    return ret;
}

}

#endif
