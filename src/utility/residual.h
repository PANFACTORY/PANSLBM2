#pragma once
#include <cmath>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

namespace PANSLBM2 {
    template<class T>
    T Residual(const T *_ux, const T *_uy, const T *_uxp, const T *_uyp, int _nxy) {
        T norm_buffer[2] = { T() }, norm[2] = { T() };
        for (int idx = 0; idx < _nxy; ++idx) {
            norm_buffer[0] += pow(_ux[idx] - _uxp[idx], 2.0) + pow(_uy[idx] - _uyp[idx], 2.0);
            norm_buffer[1] += pow(_ux[idx], 2.0) + pow(_uy[idx], 2.0);
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(norm_buffer, norm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(norm[0]/norm[1]);
#else
        return sqrt(norm_buffer[0]/norm_buffer[1]);
#endif
    }

    template<class T>
    T Residual(const T *_ux, const T *_uy, const T *_uz, const T *_uxp, const T *_uyp, const T *_uzp, int _nxyz) {
        T norm_buffer[2] = { T() }, norm[2] = { T() };
        for (int idx = 0; idx < _nxyz; ++idx) {
            norm_buffer[0] += pow(_ux[idx] - _uxp[idx], 2.0) + pow(_uy[idx] - _uyp[idx], 2.0) + pow(_uz[idx] - _uzp[idx], 2.0);
            norm_buffer[1] += pow(_ux[idx], 2.0) + pow(_uy[idx], 2.0) + pow(_uz[idx], 2.0);
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(norm_buffer, norm, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(norm[0]/norm[1]);
#else
        return sqrt(norm_buffer[0]/norm_buffer[1]);
#endif
    }
}