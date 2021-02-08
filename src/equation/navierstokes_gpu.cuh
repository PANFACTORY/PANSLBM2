//*****************************************************************************
//  Title       :   src/equation/navierstokes_gpu.cuh
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/04
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   Update macroscopic values, rho, ux, uy with gpu
        //*********************************************************************
        template<class T>
        __global__ void kernel_UpdateMacro(int _nc, T **_ft_d, T *_rho_d, T *_ux_d, T *_uy_d) {
            int i = threadIdx.x;
            _rho_d[i] = T();
            _ux_d[i] = T();
            _uy_d[i] = T();
            for (int j = 0; j < _nc; j++) {
                _rho_d[i] += _ft_d[j][i];
            }
            _ux_d[i] /= _rho_d[i];
            _uy_d[i] /= _rho_d[i];
        }

        template<class T, template<class>class P>
        void UpdateMacroGPU(P<T>& _p, T *_rho, T *_ux, T *_uy) {
            assert(P<T>::nd == 2);

            T *ft_d[P<T>::nc];
            for (int j = 0; j < P<T>::nc; j++) {
                cudaMalloc((void**)&ft_d[j], _p.np*sizeof(T));
                cudaMemcpy(ft_d[j], _p.ft[j], _p.np*sizeof(T), cudaMemcpyHostToDevice);
            }

            T *rho_d, *ux_d, *uy_d;
            cudaMalloc((void**)&rho_d, _p.np*sizeof(T));
            cudaMalloc((void**)&ux_d, _p.np*sizeof(T));
            cudaMalloc((void**)&uy_d, _p.np*sizeof(T));

            kernel_UpdateMacro<<< 1, 1 >>>(P<T>::nc, ft_d, rho_d, ux_d, uy_d);

            cudaMemcpy(_rho, rho_d, _p.np*sizeof(T), cudaMemcpyDeviceToHost);
            cudaFree(rho_d);
            cudaMemcpy(_ux, ux_d, _p.np*sizeof(T), cudaMemcpyDeviceToHost);
            cudaFree(ux_d);
            cudaMemcpy(_uy, uy_d, _p.np*sizeof(T), cudaMemcpyDeviceToHost);
            cudaFree(uy_d);
        }
    }
}