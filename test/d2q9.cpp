#define _USE_AVX_DEFINES
#include <stdio.h>
#include "../src/particle/d2q9.h"

using namespace PANSLBM2;

int main() {
    //  Test ShuffleToSoA
    alignas(32) double f[32] = { 
        11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 
        21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0,
        31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 
        41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0
    };
    __m256d __f[8];
    D2Q9<double>::ShuffleToSoA(f, __f);
    for (int c = 0; c < 8; ++c) {
        printf("%lf %lf %lf %lf\n", __f[c][3], __f[c][2], __f[c][1], __f[c][0]);
    }

    //  Test ShuffleToAoS
    D2Q9<double>::ShuffleToAoS(f, __f);
    for (int idx = 0; idx < 4; ++idx) {
        printf("%lf %lf %lf %lf %lf %lf %lf %lf\n", 
            f[8*idx + 0], f[8*idx + 1], f[8*idx + 2], f[8*idx + 3], 
            f[8*idx + 4], f[8*idx + 5], f[8*idx + 6], f[8*idx + 7]);
    }
}