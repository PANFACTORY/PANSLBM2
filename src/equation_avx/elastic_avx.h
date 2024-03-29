//*****************************************************************************
//  Title       :   src/equation_avx/elastic_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace EL {
        template<class T, template<class>class P>void Macro(T, T &, T &, T &, T &, T &, T &, const T *, int);       //  Function of updating macroscopic values of EL for 2D       
        template<class T, template<class>class P>void Macro(T, T &, T &, T &, T &, T &, T &, const T *, T, int);    //  Function of updating macroscopic values of EL with topology optimization for 2D
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T, T, T, T, T);                        //  Function of getting equilibrium of EL for 2D

        //  Function of updating macroscopic values of EL for 2D
        template<class P>
        void Macro(const __m256d &__rho, __m256d &__ux, __m256d & __uy, __m256d &__sxx, __m256d &__sxy, __m256d &__syx, __m256d &__syy, const __m256d *__f) {
            __ux = _mm256_setzero_pd();
            __uy = _mm256_setzero_pd();
            __sxx = _mm256_setzero_pd();
            __sxy = _mm256_setzero_pd();
            __syx = _mm256_setzero_pd();
            __syy = _mm256_setzero_pd();
            for (int c = 1; c < P::nc; ++c) {
                __ux = _mm256_add_pd(__ux, _mm256_mul_pd(P::__cx[c], __f[c]));
                __uy = _mm256_add_pd(__uy, _mm256_mul_pd(P::__cy[c], __f[c]));
                __sxx = _mm256_sub_pd(__sxx, _mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cx[c]), __f[c]));
                __sxy = _mm256_sub_pd(__sxy, _mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cy[c]), __f[c]));
                __syx = _mm256_sub_pd(__syx, _mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cx[c]), __f[c]));
                __syy = _mm256_sub_pd(__syy, _mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cy[c]), __f[c]));
            }
            __m256d __invrho = _mm256_div_pd(_mm256_set1_pd(1.0), __rho);
            __ux = _mm256_mul_pd(__ux, __invrho);
            __uy = _mm256_mul_pd(__uy, __invrho);
        }

        //  Function of updating macroscopic values of EL for 2D
        template<class P>
        void Macro(const __m256d &__rho, __m256d &__ux, __m256d & __uy, __m256d &__sxx, __m256d &__sxy, __m256d &__syx, __m256d &__syy, const __m256d *__f, const __m256d &__gamma) {
            Macro<P>(__rho, __ux, __uy, __sxx, __sxy, __syx, __syy, __f);
            __sxx = _mm256_mul_pd(__sxx, __gamma);
            __sxy = _mm256_mul_pd(__sxy, __gamma);
            __syx = _mm256_mul_pd(__syx, __gamma);
            __syy = _mm256_mul_pd(__syy, __gamma);
        }

        //  Function of getting equilibrium of EL for 2D
        template<class P>
        void Equilibrium(__m256d *__feq, const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__sxx, const __m256d &__sxy, const __m256d &__syx, const __m256d &__syy) {
            __m256d __trs = _mm256_add_pd(__sxx, __syy);
            for (int c = 0; c < P::nc; ++c) {
                __m256d __cu = _mm256_add_pd(_mm256_mul_pd(P::__cx[c], __ux), _mm256_mul_pd(P::__cy[c], __uy));
                __m256d __csc = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cx[c]), __sxx), _mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cy[c]), __sxy)), _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cx[c]), __syx), _mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cy[c]), __syy)));          
                __feq[c] = _mm256_mul_pd(__ei[c], _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_mul_pd(__rho, __cu)), _mm256_mul_pd(_mm256_set1_pd(4.5), __csc)), _mm256_mul_pd(_mm256_set1_pd(1.5), __trs)));
            }
        }

        //  Function of Update macro and Collide of NS for 2D
        template<template<class>class P>
        void MacroCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double *_sxx, double *_sxy, double *_syx, double *_syy, double _tau, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f and rho
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]);

                //  Update macro
                __m256d __ux, __uy, __sxx, __sxy, __syx, __syy;
                Macro<P<double> >(__rho, __ux, __uy, __sxx, __sxy, __syx, __syy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_sxx[idx], __sxx);
                    _mm256_storeu_pd(&_sxy[idx], __sxy);
                    _mm256_storeu_pd(&_syx[idx], __syx);
                    _mm256_storeu_pd(&_syy[idx], __syy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy, __sxx, __sxy, __syx, __syy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ux, uy, sxx, sxy, syx, syy;
                Macro<double, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _sxx[idx] = sxx;
                    _sxy[idx] = sxy;
                    _syx[idx] = syx;
                    _syy[idx] = syy;
                }

                //  Collide
                Equilibrium<double, P>(feq, _rho[idx], ux, uy, sxx, sxy, syx, syy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro and Collide of NS for 2D
        template<template<class>class P>
        void MacroExtendedCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double *_sxx, double *_sxy, double *_syx, double *_syy, double _tau, const T *_gamma, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f and rho
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __gamma = _mm256_loadu_pd(&_gamma[idx]); 

                //  Update macro
                __m256d __ux, __uy, __sxx, __sxy, __syx, __syy;
                Macro<P<double> >(__rho, __ux, __uy, __sxx, __sxy, __syx, __syy, __f, __gamma);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_sxx[idx], __sxx);
                    _mm256_storeu_pd(&_sxy[idx], __sxy);
                    _mm256_storeu_pd(&_syx[idx], __syx);
                    _mm256_storeu_pd(&_syy[idx], __syy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy, __sxx, __sxy, __syx, __syy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ux, uy, sxx, sxy, syx, syy;
                Macro<double, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, _p.f, _gamma[idx], idx);

                //  Save macro if need
                if (_issave) {
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _sxx[idx] = sxx;
                    _sxy[idx] = sxy;
                    _syx[idx] = syx;
                    _syy[idx] = syy;
                }

                //  Collide
                Equilibrium<double, P>(feq, _rho[idx], ux, uy, sxx, sxy, syx, syy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }
    }
}