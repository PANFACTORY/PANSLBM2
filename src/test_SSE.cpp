#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <xmmintrin.h>

void verify(const int n, const float a, const float* x, const float* y) {
    float* vY = new float[n];

    for (int i = 0; i < n; i++) {
        vY[i] = a*x[i];
    }

    for (int i = 0; i < n; i++) {
        if (y[i] != vY[i]) {
            printf("error, i = %d, y = %f, vY = %f\n", i, y[i], vY[i]);
        }
    }

    delete[] vY;
}

int main() {
    const int loop = 65536, n = 65536;

    alignas(16) volatile float a = (float)rand();
    alignas(16) float x[n], y[n];
    
    for (int i = 0; i < n; i++) {
        x[i] = (float)rand();
    }

clock_t startTime = clock();

    for (int j = 0; j < loop; j++) {
        /*for (int i = 0; i < n; i++) {
            y[i] = a*x[i];
        }*/

        __m128 va = _mm_load1_ps((const float*)&a);
        for (int i = 0; i < n; i += sizeof(va)/sizeof(float)) {
            __m128 vx = _mm_loadu_ps(&x[i]);
            __m128 vr = _mm_mul_ps(va, vx);
            _mm_storeu_ps(&y[i], vr);
        }
    }

clock_t stopTime = clock();

    float etime = (float)(stopTime - startTime)/CLOCKS_PER_SEC;
    printf("%15.7f sec", etime);

    verify(n, a, x, y);
    return 0;
}