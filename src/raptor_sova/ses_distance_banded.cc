#include <raptor_sova/ses_distance_banded.h>

#include <algorithm>
#include <cmath>
#include <iostream>

namespace raptor {
namespace ses {

int32_t BandedSESDistance(const std::string& q, const std::string& t, double maxd_frac, double bandw_frac) {
    return BandedSESDistance(q.c_str(), q.size(), t.c_str(), t.size(), maxd_frac, bandw_frac);
}

int32_t BandedSESDistance(const char* q, size_t qlen, const char* t, size_t tlen, double maxd_frac, double bandw_frac) {
    int32_t N = qlen;
    int32_t M = tlen;
    int32_t d_max = std::min(N + M, static_cast<int32_t>(maxd_frac * (N + M)));
    int32_t band_w = qlen * bandw_frac;

    int32_t zero_offset = d_max + 1;
    std::vector<int32_t> v(2 * d_max + 3, MINUS_INF);

    { // Initialization is required. Outside of the main loop, so that the
        // MINUS_INF trick can be used.
        int32_t d = 0, k = 0, x = 0, y = 0;
        while (x < N && y < M && q[x] == t[y]) {
            x += 1;
            y += 1;
        }
        v[zero_offset] = x;

        if (x >= N && y >= M) {
            return d;
        }
    }

    std::vector<int32_t> u(2 * d_max + 3, MINUS_INF);

    int32_t band_tolerance = band_w / 2 + 1;
    int32_t min_k = -1, max_k = 1;  // +- 1 because we handled the '0' case above.
    int32_t best_u = 0;

    for (int32_t d = 1; d < d_max; d++) {
        if ((max_k - min_k) > band_w) {
            break;
        }

        // printf ("min_k = %d, max_k = %d\n", min_k, max_k);
        for (int32_t k = min_k; k <= max_k; k += 2) {
            int32_t kz = k + zero_offset;
            int32_t x = std::max(v[kz + 1], v[kz - 1] + 1);

            int32_t y = x - k;
            // printf ("x = %d, y = %d, N = %d, M = %d\n", x, y, N, M);
            while (x < N && y < M && q[x] == t[y]) {
                x += 1;
                y += 1;
            }

            v[kz] = x;
            u[kz] = x + y;
            best_u = std::max(u[kz], best_u);

            if (x >= N || y >= M) {
                return d;
            }
        }

        int32_t new_min_k = max_k;
        int32_t new_max_k = min_k;
        for (int32_t k = min_k; k <= max_k; k+=2) {
            if (u[k + zero_offset] >= (best_u - band_tolerance)) {
                new_min_k = std::min(k, new_min_k);
                new_max_k = std::max(k, new_max_k);
            }
        }

        min_k = new_min_k - 1;
        max_k = new_max_k + 1;
    }

    return -1;
}

SesResults BandedSESDistanceAdvanced(const char* q, size_t qlen, const char* t, size_t tlen, double maxd_frac, double bandw_frac, int32_t match_score, int32_t indel_penalty) {
    SesResults ret;

    int32_t N = qlen;
    int32_t M = tlen;
    int32_t d_max = std::min(N + M, static_cast<int32_t>(maxd_frac * (N + M)));
    int32_t band_w = qlen * bandw_frac;
    int32_t score = 0;

    int32_t zero_offset = d_max + 1;
    std::vector<int32_t> v(2 * d_max + 3, MINUS_INF);

    { // Initialization is required. Outside of the main loop, so that the
        // MINUS_INF trick can be used.
        int32_t d = 0, k = 0, x = 0, y = 0;
        while (x < N && y < M && q[x] == t[y]) {
            ++x;
            ++y;
        }
        v[zero_offset] = x;

        ret.max_score = x * match_score;
        ret.max_q = x;
        ret.max_t = y;
        ret.max_score_diffs = 0;
        ret.diffs = 0;

        if (x >= N || y >= M) {
            ret.valid = true;
            ret.last_score = ret.max_score;
            ret.last_q = x;
            ret.last_t = y;
            return ret;
        }
    }

    std::vector<int32_t> u(2 * d_max + 3, MINUS_INF);
    std::vector<int32_t> scores(2 * d_max + 3, 0);

    int32_t band_tolerance = band_w / 2 + 1;
    int32_t min_k = -1, max_k = 1;  // +- 1 because we handled the '0' case above.
    int32_t best_u = 0;

    scores[zero_offset] = ret.max_score;

    for (int32_t d = 1; d < d_max; d++) {
        ret.diffs = d;
        if ((max_k - min_k) > band_w) {
            ret.valid = false;
            break;
        }

        for (int32_t k = min_k; k <= max_k; k += 2) {
            int32_t kz = k + zero_offset;
            auto next = v[kz + 1];
            auto prev = v[kz - 1] + 1;
            // int32_t x = std::max(v[kz + 1], v[kz - 1] + 1);
            int32_t x = next;
            scores[kz] = scores[kz + 1];
            if (prev > next) {
                x = prev;
                scores[kz] = scores[kz - 1];
            }
            int32_t y = x - k;
            scores[kz] += indel_penalty;
            while (x < N && y < M && q[x] == t[y]) {
                ++x;
                ++y;
                scores[kz] += match_score;
            }

            v[kz] = x;
            u[kz] = x + y;
            best_u = std::max(u[kz], best_u);

            if (scores[kz] > ret.max_score) {
                ret.max_score = scores[kz];
                ret.max_q = x;
                ret.max_t = y;
                ret.max_score_diffs = d;
                // std::cerr << "    [x = " << x << ", y = " << y << ", k = " << k << ", min_k = " << min_k << ", max_k = " << max_k << ", band_w = " << band_w << ", d = " << d << ", d_max = " << d_max << "]: max_score = " << ret.max_score << "\n";
            }

            if (x >= N || y >= M) {
                ret.valid = true;
                ret.last_score = scores[kz];
                ret.last_q = x;
                ret.last_t = y;
                break;
            }
        }

        if (ret.valid) {
            break;
        }

        int32_t new_min_k = max_k;
        int32_t new_max_k = min_k;
        for (int32_t k = min_k; k <= max_k; k+=2) {
            if (u[k + zero_offset] >= (best_u - band_tolerance)) {
                new_min_k = std::min(k, new_min_k);
                new_max_k = std::max(k, new_max_k);
            }
        }

        min_k = new_min_k - 1;
        max_k = new_max_k + 1;
    }

    return ret;
}

}
}
