#define _USE_AVX_DEFINES
#include <stdio.h>
#include "../src/particle/d3q15.h"

using namespace PANSLBM2;

int main() {
    //  Test ShuffleToSoA
    alignas(32) double f[56] = { 
          0.1,   0.2,   0.3,   0.4,   0.5,   0.6,   0.7,   0.8,   0.9,    1.0,    1.1,    1.2,    1.3,    1.4,
         10.0,  20.0,  30.0,  40.0,  50.0,  60.0,  70.0,  80.0,  90.0,  100.0,  110.0,  120.0,  130.0,  140.0,
         -0.1,  -0.2,  -0.3,  -0.4,  -0.5,  -0.6,  -0.7,  -0.8,  -0.9,   -1.0,   -1.1,   -1.2,   -1.3,   -1.4,
        -10.0, -20.0, -30.0, -40.0, -50.0, -60.0, -70.0, -80.0, -90.0, -100.0, -110.0, -120.0, -130.0, -140.0
    };
    __m256d __f[14];
    D3Q15<double>::ShuffleToSoA(f, __f);
    for (int c = 0; c < 14; ++c) {
        printf("%lf %lf %lf %lf\n", __f[c][3], __f[c][2], __f[c][1], __f[c][0]);
    }

    //  Test ShuffleToAoS
    D3Q15<double>::ShuffleToAoS(f, __f);
    for (int idx = 0; idx < 4; ++idx) {
        printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
            f[14*idx +  0], f[14*idx +  1], f[14*idx +  2], f[14*idx +  3], f[14*idx +  4], f[14*idx +  5], f[14*idx +  6], 
            f[14*idx +  7], f[14*idx +  8], f[14*idx +  9], f[14*idx + 10], f[14*idx + 11], f[14*idx + 12], f[14*idx + 13]
        );
    }
}