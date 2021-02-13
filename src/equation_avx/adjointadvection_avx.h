//*****************************************************************************
//  Title       :   src/equation_avx/adjointadvection_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx


namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection   :   External force with heat generation
        //*********************************************************************
        template<template<class>class P>
        void ExternalForceHeatgeneration(P<double>& _particle, double *_tem, double *_beta) {
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;
            
            double a = 1.0;
            __m256d __a = _mm256_broadcast_sd((const double*)&a), __dx = _mm256_broadcast_sd((const double*)&_particle.dx);
            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = _mm256_loadu_pd(&_tem[i]), __beta = _mm256_loadu_pd(&_beta[i]);
                _mm256_storeu_pd(&_tem[i], _mm256_div_pd(_mm256_add_pd(__tem, _mm256_mul_pd(__dx, __beta)), _mm256_add_pd(__a, _mm256_mul_pd(__dx, __beta))));
            }
            
            for (int i = ne; i < _particle.np; i++) {
                _tem[i] = (_tem[i] + _particle.dx*_beta[i])/(1.0 + _particle.dx*_beta[i]);
            }
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint navier-stokes 2D    :   External force with heat exchange
        //*********************************************************************
        template<template<class>class P>
        void ExternalForceHeatexchange(double _diffusivity, P<double>& _particlef, double *_rho, double *_ux, double *_uy, double *_tem, double *_iqx, double *_iqy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particlef.np/packsize)*packsize;

            double omega = 1.0/(3.0*_diffusivity*_particlef.dt/(_particlef.dx*_particlef.dx) + 0.5), a = 3.0;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __a = _mm256_broadcast_sd((const double*)&a);

            __m256d __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = _mm256_loadu_pd(&_rho[i]), __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __tem = _mm256_loadu_pd(&_tem[i]), __iqx = _mm256_loadu_pd(&_iqx[i]), __iqy = _mm256_loadu_pd(&_iqy[i]); 

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particlef.ft[j][i]);
                    _mm256_storeu_pd(&_particlef.ft[j][i], _mm256_add_pd(__ft, _mm256_mul_pd(__a, _mm256_mul_pd(__tem, _mm256_mul_pd(__omega, _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(__cx[j], __ux), __iqx), _mm256_mul_pd(_mm256_sub_pd(__cy[j], __uy), __iqy)), __rho))))));
                }
            }

            for (int i = ne; i < _particlef.np; i++) {
                for (int j = 0; j < P<double>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_tem[i]*omega*((P<double>::cx[j] - _ux[i])*_iqx[i] + (P<double>::cy[j] - _uy[i])*_iqy[i])/_rho[i];
                }
            }
        }
    }

    namespace AAD {
        //*********************************************************************
        //  Adjoint advection 2D    :   Update macroscopic values
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _particle, double *_item, double *_iqx, double *_iqy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j] ,cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __item = { 0 }, __iqx = { 0 }, __iqy = { 0 };
                
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __item = _mm256_add_pd(__item, _mm256_mul_pd(__ei[j], __ft));
                    __iqx = _mm256_add_pd(__iqx, _mm256_mul_pd(__ei[j], _mm256_mul_pd(__cx[j], __ft)));
                    __iqy = _mm256_add_pd(__iqy, _mm256_mul_pd(__ei[j], _mm256_mul_pd(__cy[j], __ft)));
                }

                _mm256_storeu_pd(&_item[i], __item);
                _mm256_storeu_pd(&_iqx[i], __iqx);
                _mm256_storeu_pd(&_iqy[i], __iqy);
            }

            for (int i = ne; i < _particle.np; i++) {
                _item[i] = 0.0;
                _iqx[i] = 0.0;
                _iqy[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _item[i] += P<double>::ei[j]*_particle.ft[j][i];
                    _iqx[i] += P<double>::ei[j]*P<double>::cx[j]*_particle.ft[j][i];
                    _iqy[i] += P<double>::ei[j]*P<double>::cy[j]*_particle.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _diffusivity, P<double>& _particle, double *_ux, double *_uy, double *_item, double *_iqx, double *_iqy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5), iomega = 1.0 - omega, a = 3.0;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega), __a = _mm256_broadcast_sd((const double*)&a);

            for (int i = 0; i < ne; i += packsize) {
                __m256d __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __item = _mm256_loadu_pd(&_item[i]), __iqx = _mm256_loadu_pd(&_iqx[i]), __iqy = _mm256_loadu_pd(&_iqy[i]);
                __m256d __feq = _mm256_add_pd(__item, _mm256_mul_pd(__a, _mm256_add_pd(_mm256_mul_pd(__ux, __iqx), _mm256_mul_pd(__uy, __iqy)))); 
                
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    _mm256_storeu_pd(&_particle.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }

            for (int i = ne; i < _particle.np; i++) {
                double feq = _item[i] + 3.0*(_ux[i]*_iqx[i] + _uy[i]*_iqy[i]);
                for (int j = 0; j < P<double>::nc; j++) {
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   External force
        //*********************************************************************
        template<template<class>class P>
        void ExternalForceHeatexchange(P<double>& _particle, double *_item, double *_beta) {
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double a = 1.0;
            __m256d __a = _mm256_broadcast_sd((const double*)&a), __dx = _mm256_broadcast_sd((const double*)&_particle.dx);
            for (int i = 0; i < ne; i += packsize) {
                __m256d __item = _mm256_loadu_pd(&_item[i]), __beta = _mm256_loadu_pd(&_beta[i]);
                _mm256_storeu_pd(&_item[i], _mm256_div_pd(_mm256_sub_pd(__item, __beta), _mm256_add_pd(__a, _mm256_mul_pd(__dx, __beta))));
            }
            
            for (int i = ne; i < _particle.np; i++) {
                _item[i] = (_item[i] - _beta[i])/(1.0 + _particle.dx*_beta[i]);
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _ux, T _uy, T _item, T _iqx, T _iqy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                _particle.ft[j][_i] = _item + 3.0*(_ux*_iqx + _uy*_iqy);
            }
        }
    }
}