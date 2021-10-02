//*****************************************************************************
//  Title       :   src/equation_avx/adjointelastic_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/14
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace AEL {
        template<class T, template<class>class P>void Macro(T &, T &, T &, T &, T &, T &, T &, const T *, const T *, int);  //  Function of updating macroscopic values of AEL for 2D
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T, T, T, T, T, T);                             //  Function of getting equilibrium of AEL for 2D

        //  Function of updating macroscopic values of AEL for 2D
        template<class P>
        void Macro(__m256d &__irho, __m256d &__imx, __m256d &__imy, __m256d &__isxx, __m256d &__isxy, __m256d &__isyx, __m256d &__isyy, const __m256d *__f) {
            __irho = _mm256_mul_pd(P::__ei[0], __f[0]);
            __imx = _mm256_setzero_pd();
            __imy = _mm256_setzero_pd();
            __isxx = _mm256_setzero_pd();
            __isxy = _mm256_setzero_pd();
            __isyx = _mm256_setzero_pd();
            __isyy = _mm256_setzero_pd();
            for (int c = 1; c < P::nc; ++c) {
                __irho = _mm256_add_pd(__irho, _mm256_mul_pd(P::__ei[c], __f[c]));
                __imx = _mm256_add_pd(__imx, _mm256_mul_pd(_mm256_mul_pd(P::__ei[c], P::__cx[c]), __f[c]));
                __imy = _mm256_add_pd(__imy, _mm256_mul_pd(_mm256_mul_pd(P::__ei[c], P::__cy[c]), __f[c]));
                __isxx = _mm256_add_pd(__isxx, _mm256_mul_pd(P::__ei[c], _mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cx[c]), __f[c])));
                __isxy = _mm256_add_pd(__isxy, _mm256_mul_pd(P::__ei[c], _mm256_mul_pd(_mm256_mul_pd(P::__cx[c], P::__cy[c]), __f[c])));
                __isyx = _mm256_add_pd(__isyx, _mm256_mul_pd(P::__ei[c], _mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cx[c]), __f[c])));
                __isyy = _mm256_add_pd(__isyy, _mm256_mul_pd(P::__ei[c], _mm256_mul_pd(_mm256_mul_pd(P::__cy[c], P::__cy[c]), __f[c])));
            }
        }

        //  Function of getting equilibrium of AEL for 2D
        template<class P>
        __m256d Equilibrium(__m256d __irho, __m256d __imx, __m256d __imy, __m256d __isxx, __m256d __isxy, __m256d __isyx, __m256d __isyy, __m256d __gamma, int _c) {
            __m256d __3 = _mm256_set1_pd(3.0), __45 = _mm256_set1_pd(4.5), __15 = _mm256_set1_pd(1.5);
            __m256d __mic = _mm256_add_pd(_mm256_mul_pd(__imx, P::__cx[_c]), _mm256_mul_pd(__imy, P::__cy[_c]));
            __m256d __csc = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(P::__cx[_c], P::__cx[_c]), __isxx), _mm256_mul_pd(_mm256_mul_pd(P::__cx[_c], P::__cy[_c]), __isxy)), _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(P::__cy[_c], P::__cx[_c]), __isyx), _mm256_mul_pd(_mm256_mul_pd(P::__cy[_c], P::__cy[_c]), __isyy)));
            __m256d __trs = _mm256_mul_pd(__irho, _mm256_add_pd(_mm256_mul_pd(P::__cx[_c], P::__cx[_c]), _mm256_mul_pd(P::__cy[_c], P::__cy[_c])));
            return _mm256_add_pd(_mm256_mul_pd(__3, __mic), _mm256_sub_pd(_mm256_mul_pd(__45, _mm256_mul_pd(__gamma, __csc)), _mm256_mul_pd(__15, _mm256_mul_pd(__gamma, __trs))));
        }

        //  Function of Update macro and Collide of AEL for 2D
        template<template<class>class P>
        void MacroCollide(P<double>& _p, double *_irho, double *_imx, double *_imy, double *_isxx, double *_isxy, double *_isyx, double *_isyy, double _tau, const double *_gamma, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega);
            #pragma omp parallel for private(feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f and rho
                __m256d __f[P<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __gamma = _mm256_loadu_pd(&_gamma[idx]); 

                //  Update macro
                __m256d __irho, __imx, __imy, __isxx, __isxy, __isyx, __isyy;
                Macro<P<double> >(__irho, __imx, __imy, __isxx, __isxy, __isyx, __isyy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_irho[idx], __irho);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_isxx[idx], __isxx);
                    _mm256_storeu_pd(&_isxy[idx], __isxy);
                    _mm256_storeu_pd(&_isyx[idx], __isyx);
                    _mm256_storeu_pd(&_isyy[idx], __isyy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[0]), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__irho, __imx, __imy, __isxx, __isxy, __isyx, __isyy, __gamma, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__irho, __imx, __imy, __isxx, __isxy, __isyx, __isyy, __gamma, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double irho, imx, imy, isxx, isxy, isyx, isyy;
                Macro<double, P>(irho, imx, imy, isxx, isxy, isyx, isyy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _irho[idx] = irho;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _isxx[idx] = isxx;
                    _isxy[idx] = isxy;
                    _isyx[idx] = isyx;
                    _isyy[idx] = isyy;
                }

                //  Collide
                Equilibrium<double, P>(feq, irho, imx, imy, isxx, isxy, isyx, isyy, _gamma[idx]);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }
    }
}