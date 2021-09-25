//*****************************************************************************
//  Title       :   src/equation_avx/navierstokes_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/01/24
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace NS_AVX {
        //  Function of updating macroscopic values of NS for 2D
        template<class P>
        void Macro(__m256d &__rho, __m256d &__ux, __m256d &__uy, __m256d __f0, const __m256d *__f) {
            __rho = __f0;
            __ux = _mm256_setzero_pd();
            __uy = _mm256_setzero_pd();
            for (int c = 1; c < P::nc; ++c) {
                __rho = _mm256_add_pd(__rho, __f[c]);
                __ux = _mm256_add_pd(__ux, _mm256_mul_pd(__f[c], P::__cx[c]));
                __uy = _mm256_add_pd(__uy, _mm256_mul_pd(__f[c], P::__cy[c]));
            }
            __ux = _mm256_div_pd(__ux, __rho);
            __uy = _mm256_div_pd(__uy, __rho);
        }

        //  Function of getting equilibrium of NS for 2D
        template<class P>
        __m256d Equilibrium(__m256d __rho, __m256d __ux, __m256d __uy, int _c) {
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0), __45 = _mm256_set1_pd(4.5), __15 = _mm256_set1_pd(1.5);
            __m256d __1m15uu = _mm256_sub_pd(__1, _mm256_mul_pd(__15, _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy))));
            __m256d __cu = _mm256_add_pd(_mm256_mul_pd(P::__cx[_c], __ux), _mm256_mul_pd(P::__cy[_c], __uy));
            __m256d __3cup45cucu = _mm256_add_pd(_mm256_mul_pd(__3, __cu), _mm256_mul_pd(__45, _mm256_mul_pd(__cu, __cu)));
            return _mm256_mul_pd(P::__ei[_c], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, __3cup45cucu)));
        }

        //  Function of Update macro, Collide and Stream of NS for 2D
        template<template<class>class P>
        void MacroCollideStream(P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity, bool _issave = false) {
            double omega = 1.0/(3.0*_viscosity + 0.5);
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(1.0 - omega);

            for (int idx = 0; idx < _p.nxyz; idx += P<double>::packsize) {
                //  Pack f0 and f
                __m256d __f0 = _mm256_loadu_pd(&_p.f0[idx]), __f[P<double>::nc] = { 0 };
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_loadu_pd(&_p.f[P<double>::IndexF(idx, c)]);
                }

                //  Update macro
                __m256d __rho, __ux, __uy;
                Macro<P<double> >(__rho, __ux, __uy, __f0, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                }

                //  Collide
                _mm256_storeu_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f0), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    _mm256_storeu_pd(&_p.f[P<double>::IndexF(idx, c)], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, c))));
                }
            }
        }
    }
}