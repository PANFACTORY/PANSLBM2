//*****************************************************************************
//  Title       :   src/equation_avx/elastic_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace EL {
        //*********************************************************************
        //  Elastic 2D  :   Update macroscopic values
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _p, double *_rho, double *_ux, double *_uy, double *_sxx, double *_sxy, double *_syx, double *_syy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_p.np/packsize)*packsize;

            __m256d __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = _mm256_loadu_pd(&_rho[i]), __ux = { 0 }, __uy = { 0 }, __sxx = { 0 }, __sxy = { 0 }, __syx = { 0 }, __syy = { 0 };

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_p.ft[j][i]);
                    __ux = _mm256_add_pd(__ux, _mm256_mul_pd(__cx[j], __ft));
                    __uy = _mm256_add_pd(__uy, _mm256_mul_pd(__cy[j], __ft));
                    __sxx = _mm256_sub_pd(__sxx, _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cx[j]), __ft));
                    __sxy = _mm256_sub_pd(__sxy, _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cy[j]), __ft));
                    __syx = _mm256_sub_pd(__syx, _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cx[j]), __ft));
                    __syy = _mm256_sub_pd(__syy, _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cy[j]), __ft));
                }

                _mm256_storeu_pd(&_ux[i], _mm256_div_pd(__ux, __rho));
                _mm256_storeu_pd(&_uy[i], _mm256_div_pd(__uy, __rho));
                _mm256_storeu_pd(&_sxx[i], __sxx);
                _mm256_storeu_pd(&_sxy[i], __sxy);
                _mm256_storeu_pd(&_syx[i], __syx);
                _mm256_storeu_pd(&_syy[i], __syy);
            }

            for (int i = ne; i < _p.np; i++) {
                _ux[i] = 0.0;
                _uy[i] = 0.0;
                _sxx[i] = 0.0;
                _sxy[i] = 0.0;
                _syx[i] = 0.0;
                _syy[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _ux[i] += P<double>::cx[j]*_p.ft[j][i];
                    _uy[i] += P<double>::cy[j]*_p.ft[j][i];
                    _sxx[i] -= P<double>::cx[j]*P<double>::cx[j]*_p.ft[j][i];
                    _sxy[i] -= P<double>::cx[j]*P<double>::cy[j]*_p.ft[j][i];
                    _syx[i] -= P<double>::cy[j]*P<double>::cx[j]*_p.ft[j][i];
                    _syy[i] -= P<double>::cy[j]*P<double>::cy[j]*_p.ft[j][i]; 
                }
                _ux[i] /= _rho[i];
                _uy[i] /= _rho[i];
            }
        }

        //*********************************************************************
        //  Elastic 2D  :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _elasticy, P<double>& _p, double *_rho, double *_ux, double *_uy, double *_sxx, double *_sxy, double *_syx, double *_syy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_p.np/packsize)*packsize;

            double omega = 1.0/(3.0*_diffusivity*_p.dt/(_p.dx*_p.dx) + 0.5), iomega = 1.0 - omega, a = 3.0, b = 4.5, c = 1.5;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega), __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b), __c = _mm256_broadcast_sd((const double*)&c);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = _mm256_loadu_pd(&_rho[i]), __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __sxx = _mm256_loadu_pd(&_sxx[i]), __sxy = _mm256_loadu_pd(&_sxy[i]), __syx = _mm256_loadu_pd(&_syx[i]), __syy = _mm256_loadu_pd(&_syy[i]);
                __m256d __trs = _mm256_add_pd(__sxx, __syy);
                
                for (int j = 0; j < P<T>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_p.ft[j][i]);
                    __m256d __ciu = _mm256_add_pd(_mm256_mul_pd(__cx[j], __ux), _mm256_mul_pd(__cy[j], __uy));
                    __m256d __csc = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(__cx[j], __cx[j]), __sxx), _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cy[j]), __sxy)), _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(__cy[j], __cx[j]), __syx), _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cy[j]), __syy)));
                    
                    __m256d __feq = _mm256_mul_pd(__ei[j], _mm256_add_pd(_mm256_sub_pd(_mm256_mul_pd(__a, _mm256_mul_pd(__rho, __ciu)), _mm256_mul_pd(__b, __csc)), _mm256_mul_pd(__c, __trs)));
                    _mm256_storeu_pd(&_p.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }

            for (int i = ne; i < _p.np; i++) {
                for (int j = 0; j < P<double>::nc; j++) {
                    T cu = P<double>::cx[j]*_ux[i] + P<double>::cy[j]*_uy[i];
                    T csc = P<double>::cx[j]*_sxx[i]*P<double>::cx[j] + P<double>::cx[j]*_sxy[i]*P<double>::cy[j] + P<double>::cy[j]*_syx[i]*P<double>::cx[j] + P<double>::cy[j]*_syy[i]*P<double>::cy[j];
                    T trs = _sxx[i] + _syy[i];
                    T feq = P<double>::ei[j]*(3.0*_rho[i]*cu - 4.5*csc + 1.5*trs);
                    _p.ftp1[j][i] = (1.0 - omega)*_p.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Elastic 2D  :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _p, T _rho, T _ux, T _uy, T _sxx, T _sxy, T _syx, T _syy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _p.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T cu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy;
                T csc = P<T>::cx[j]*_sxx*P<T>::cx[j] + P<T>::cx[j]*_sxy*P<T>::cy[j] + P<T>::cy[j]*_syx*P<T>::cx[j] + P<T>::cy[j]*_syy*P<T>::cy[j];
                T trs = _sxx + _syy;
                _p.ft[j][_i] = P<T>::ei[j]*(3.0*_rho*cu - 4.5*csc + 1.5*trs);
            }
        }
    }
}