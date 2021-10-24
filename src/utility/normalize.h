#pragma once
#include <cmath>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

namespace PANSLBM2 {
    template<class T>
    void Normalize(T *_v, int _size) {
        T vmax_buffer = T(), vmax;
        for (int idx = 0; idx < _size; ++idx) {
            if (vmax_buffer < fabs(_v[idx])) {
                vmax_buffer = fabs(_v[idx]);
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&vmax_buffer, &vmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
        vmax = vmax_buffer;
#endif
        for (int idx = 0; idx < _size; ++idx) {
            _v[idx] /= vmax;
        }
    }
}