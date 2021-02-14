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
    namespace EL {
        //*********************************************************************
        //  Elastic 2D  :   Expand macroscopic values
        //*********************************************************************
        template<template<class>class P>
        void ExpandMacro(P<double>& _p, double *_sxx, double *_sxy, double *_syx, double *_syy, double *_gamma) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_p.np/packsize)*packsize;

            for (int i = 0; i < ne; i += packsize) {
                __m256d __sxx = _mm256_loadu_pd(&_sxx[i]), __sxy = _mm256_loadu_pd(&_sxy[i]), __syx = _mm256_loadu_pd(&_syx[i]), __syy = _mm256_loadu_pd(&_syy[i]), __gamma = _mm256_loadu_pd(&_gamma[i]);
                _mm256_storeu_pd(&_sxx[i], _mm256_mul_pd(__sxx, __gamma));
                _mm256_storeu_pd(&_sxy[i], _mm256_mul_pd(__sxy, __gamma));
                _mm256_storeu_pd(&_syx[i], _mm256_mul_pd(__syx, __gamma));
                _mm256_storeu_pd(&_syy[i], _mm256_mul_pd(__syy, __gamma));
            }

            for (int i = ne; i < _p.np; i++) {
                _sxx[i] *= _gamma[i];
                _sxy[i] *= _gamma[i];
                _syx[i] *= _gamma[i];
                _syy[i] *= _gamma[i];
            }
        }
    }

    namespace AEL {
        //*********************************************************************
        //  Adjoint elastic 2D  :   Update macroscopic values
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _p, double *_irho, double *_imx, double *_imy, double *_isxx, double *_isxy, double *_isyx, double *_isyy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_p.np/packsize)*packsize;

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __irho = { 0 }, __imx = { 0 }, __imy = { 0 }, __isxx = { 0 }, __isxy = { 0 }, __isyx = { 0 }, __isyy = { 0 };
                
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_p.ft[j][i]);
                    __irho = _mm256_add_pd(__irho, _mm256_mul_pd(__ei[j], __ft));
                    __imx = _mm256_add_pd(__imx, _mm256_mul_pd(_mm256_mul_pd(__ei[j], __cx[j]), __ft));
                    __imy = _mm256_add_pd(__imy, _mm256_mul_pd(_mm256_mul_pd(__ei[j], __cy[j]), __ft));
                    __isxx = _mm256_add_pd(__isxx, _mm256_mul_pd(__ei[j], _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cx[j]), __ft)));
                    __isxy = _mm256_add_pd(__isxy, _mm256_mul_pd(__ei[j], _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cy[j]), __ft)));
                    __isyx = _mm256_add_pd(__isyx, _mm256_mul_pd(__ei[j], _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cx[j]), __ft)));
                    __isyy = _mm256_add_pd(__isyy, _mm256_mul_pd(__ei[j], _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cy[j]), __ft)));
                }

                _mm256_storeu_pd(&_irho[i], __irho);
                _mm256_storeu_pd(&_imx[i], __imx);
                _mm256_storeu_pd(&_imy[i], __imy);
                _mm256_storeu_pd(&_isxx[i], __isxx);
                _mm256_storeu_pd(&_isxy[i], __isxy);
                _mm256_storeu_pd(&_isyx[i], __isyx);
                _mm256_storeu_pd(&_isyy[i], __isyy);
            }

            for (int i = ne; i < _p.np; i++) {
                _irho[i] = 0.0;
                _imx[i] = 0.0;
                _imy[i] = 0.0;
                _isxx[i] = 0.0;
                _isxy[i] = 0.0;
                _isyx[i] = 0.0;
                _isyy[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _irho[i] += P<double>::ei[j]*_p.ft[j][i];
                    _imx[i] += P<double>::ei[j]*P<double>::cx[j]*_p.ft[j][i];
                    _imy[i] += P<double>::ei[j]*P<double>::cy[j]*_p.ft[j][i];
                    _isxx[i] += P<double>::ei[j]*P<double>::cx[j]*P<double>::cx[j]*_p.ft[j][i];
                    _isxy[i] += P<double>::ei[j]*P<double>::cx[j]*P<double>::cy[j]*_p.ft[j][i];
                    _isyx[i] += P<double>::ei[j]*P<double>::cy[j]*P<double>::cx[j]*_p.ft[j][i];
                    _isyy[i] += P<double>::ei[j]*P<double>::cy[j]*P<double>::cy[j]*_p.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint elastic 2D  :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _elasticy, P<double>& _p, double *_irho, double *_imx, double *_imy, double *_isxx, double *_isxy, double *_isyx, double *_isyy, double *_gamma) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_p.np/packsize)*packsize;

            double omega = 1.0/(3.0*_elasticy*_p.dt/(_p.dx*_p.dx) + 0.5), iomega = 1.0 - omega, a = 3.0, b = 4.5, c = 1.5;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega), __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b), __c = _mm256_broadcast_sd((const double*)&c);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __irho = _mm256_loadu_pd(&_irho[i]), __imx = _mm256_loadu_pd(&_imx[i]), __imy = _mm256_loadu_pd(&_imy[i]), __gamma = _mm256_loadu_pd(&_gamma[i]);
                __m256d __isxx = _mm256_loadu_pd(&_isxx[i]), __isxy = _mm256_loadu_pd(&_isxy[i]), __isyx = _mm256_loadu_pd(&_isyx[i]), __isyy = _mm256_loadu_pd(&_isyy[i]);
                
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_p.ft[j][i]);
                    __m256d __mic = _mm256_add_pd(_mm256_mul_pd(__imx, __cx[j]), _mm256_mul_pd(__imy, __cy[j]));
                    __m256d __csc = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(__cx[j], __cx[j]), __isxx), _mm256_mul_pd(_mm256_mul_pd(__cx[j], __cy[j]), __isxy)), _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(__cy[j], __cx[j]), __isyx), _mm256_mul_pd(_mm256_mul_pd(__cy[j], __cy[j]), __isyy)));
                    __m256d __trs = _mm256_mul_pd(__irho, _mm256_add_pd(_mm256_mul_pd(__cx[j], __cx[j]), _mm256_mul_pd(__cy[j], __cy[j])));
                    __m256d __feq = _mm256_add_pd(_mm256_mul_pd(__a, __mic), _mm256_sub_pd(_mm256_mul_pd(__b, _mm256_mul_pd(__gamma, __csc)), _mm256_mul_pd(__c, _mm256_mul_pd(__gamma, __trs))));
                    _mm256_storeu_pd(&_p.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }
            
            for (int i = ne; i < _p.np; i++) {
                for (int j = 0; j < P<double>::nc; j++) {
                    double imc = _imx[i]*P<double>::cx[j] + _imy[i]*P<double>::cy[j];
                    double cisc = P<double>::cx[j]*_isxx[i]*P<double>::cx[j] + P<double>::cx[j]*_isxy[i]*P<double>::cy[j] + P<double>::cy[j]*_isyx[i]*P<double>::cx[j] + P<double>::cy[j]*_isyy[i]*P<double>::cy[j];
                    double irhocc = _irho[i]*(P<double>::cx[j]*P<double>::cx[j] + P<double>::cy[j]*P<double>::cy[j]);
                    double feq = 3.0*imc + 4.5*_gamma[i]*cisc - 1.5*_gamma[i]*irhocc;
                    _p.ftp1[j][i] = (1.0 - omega)*_p.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint elastic 2D  :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _p, T _irho, T _imx, T _imy, T _isxx, T _isxy, T _isyx, T _isyy, T _gamma) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _p.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T imc = _imx*P<T>::cx[j] + _imy*P<T>::cy[j];
                T cisc = P<T>::cx[j]*_isxx*P<T>::cx[j] + P<T>::cx[j]*_isxy*P<T>::cy[j] + P<T>::cy[j]*_isyx*P<T>::cx[j] + P<T>::cy[j]*_isyy*P<T>::cy[j];
                T irhocc = _irho*(P<T>::cx[j]*P<T>::cx[j] + P<T>::cy[j]*P<T>::cy[j]);
                _p.ft[j][_i] = 3.0*imc + 4.5*_gamma*cisc -1.5*_gamma*irhocc;
            }
        }
    }
}