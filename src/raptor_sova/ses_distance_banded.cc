#include <raptor_sova/ses_distance_banded.h>

#include <algorithm>
#include <cmath>

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

}
}
