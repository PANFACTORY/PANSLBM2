#define _USE_AVX_DEFINES
#include <stdio.h>
#include <cmath>
#include "../src/particle/d3q15.h"

using namespace PANSLBM2;

int main() {
    D3Q15<double> p(2, 2);
    for (int idx = 0; idx < p.nxyz; ++idx) {
        p.f0[idx] = pow(100, idx);
        for (int c = 1; c < p.nc; ++c) {
            p.f[p.IndexF(idx, c)] = pow(100, idx) + c;
        }
    }
    __m256d __f[p.nc];
    p.LoadF(0, __f);
    for (int c = 0; c < p.nc; ++c) {
        printf("%lf %lf %lf %lf\n", __f[c][3], __f[c][2], __f[c][1], __f[c][0]);
    }
    p.StoreF(0, __f);
    for (int idx = 0; idx < p.nxyz; ++idx) {
        printf("%lf", p.f0[idx]);
        for (int c = 1; c < p.nc; ++c) {
            printf(" %lf", p.f[p.IndexF(idx, c)]);
        }
        printf("\n");
    }
}