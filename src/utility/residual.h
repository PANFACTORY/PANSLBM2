#pragma once
#include <cmath>
#include "mpi.h"

namespace PANSLBM2 {
    template<class T>
    T residual(const T *_ux, const T *_uy, const T *_uxp, const T *_uyp, int _nxy) {
        T norm_buffer[2] = { T() }, norm[2] = { T() };
        for (int idx = 0; idx < _nxy; ++idx) {
            norm_buffer[0] += pow(_ux[idx] - _uxp[idx], 2.0) + pow(_uy[idx] - _uyp[idx], 2.0);
            norm_buffer[1] += pow(_ux[idx], 2.0) + pow(_uy[idx], 2.0);
        }
        MPI_Allreduce(norm_buffer, norm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(norm[0]/norm[1]);
    }
}