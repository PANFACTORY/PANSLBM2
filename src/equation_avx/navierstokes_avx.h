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
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   Update macroscopic values, rho, ux, uy
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _particle, double* _rho, double* _ux, double* _uy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;
            
            __m256d __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = { 0 }, __ux = { 0 }, __uy = { 0 };

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __rho = _mm256_add_pd(__rho, __ft);
                    __ux = _mm256_add_pd(__ux, _mm256_mul_pd(__ft, __cx[j]));
                    __uy = _mm256_add_pd(__uy, _mm256_mul_pd(__ft, __cy[j]));
                }

                _mm256_storeu_pd(&_rho[i], __rho);
                _mm256_storeu_pd(&_ux[i], _mm256_div_pd(__ux, __rho));
                _mm256_storeu_pd(&_uy[i], _mm256_div_pd(__uy, __rho));
            }

            for (int i = ne; i < _particle.np; i++) {
                _rho[i] = 0.0;
                _ux[i] = 0.0;
                _uy[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _rho[i] += _particle.ft[j][i];
                    _ux[i] += P<double>::cx[j]*_particle.ft[j][i];
                    _uy[i] += P<double>::cy[j]*_particle.ft[j][i];
                }
                _ux[i] /= _rho[i];
                _uy[i] /= _rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _viscosity, P<double>& _particle, double* _rho, double* _ux, double* _uy) {
            assert(P<double>::nd == 2);
            
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5), iomega = 1.0 - omega;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega);

            double a = 1.0, b = 3.0, c = 4.5, d = 1.5;
            __m256d __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b), __c = _mm256_broadcast_sd((const double*)&c), __d = _mm256_broadcast_sd((const double*)&d);

            __m256d __cx[P<double>::nc], __cy[P<double>::nc], __ej[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
                __ej[j] = _mm256_broadcast_sd((const double*)&P<double>::ei[j]);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = _mm256_loadu_pd(&_rho[i]), __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]);
                __m256d __1m15uu = _mm256_sub_pd(__a, _mm256_mul_pd(__d, _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy))));

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __cju = _mm256_add_pd(_mm256_mul_pd(__cx[j], __ux), _mm256_mul_pd(__cy[j], __uy));
                    __m256d __3cup45cucu = _mm256_add_pd(_mm256_mul_pd(__b, __cju), _mm256_mul_pd(__c, _mm256_mul_pd(__cju, __cju)));
                    __m256d __feq = _mm256_mul_pd(__ej[j], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, __3cup45cucu)));
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);                    
                    _mm256_storeu_pd(&_particle.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }

            for (int i = ne; i < _particle.np; i++) {
                double uu = 1.0 - 1.5*(_ux[i]*_ux[i] + _uy[i]*_uy[i]);
                for (int j = 0; j < P<double>::nc; j++) {
                    double ciu = P<double>::cx[j]*_ux[i] + P<double>::cy[j]*_uy[i];
                    double feq = P<double>::ei[j]*_rho[i]*(3.0*ciu + 4.5*ciu*ciu + uu);
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _rho, T _ux, T _uy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy; 
                T uu = _ux*_ux + _uy*_uy;
                _particle.ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }
    }
}