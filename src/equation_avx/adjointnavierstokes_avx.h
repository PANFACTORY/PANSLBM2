//*****************************************************************************
//  Title       :   src/equation_avx/adjointnavierstokes_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/11
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>
#include <numeric>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<template<class>class P>
        void ExternalForceBrinkman(P<double>& _particle, double *_rho, double *_ux, double *_uy, double *_alphax, double *_alphay) {
            assert(P<T>::nd == 2);
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double one = 1.0;
            __mm256d __one = _mm256_broadcast_sd((const double)&one), __dx = _mm256_broadcast_sd((const double)&_particle.dx);
            for (int i = 0; i < ne; i += packsize) {
                __m256d __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __rho = _mm256_loadu_pd(&_rho[i]), __alphax = _mm256_loadu_pd(&_alphax[i]), __alphay = _mm256_loadu_pd(&_alphay[i]);
                _mm256_storeu_pd(&_ux[i], _mm256_div_pd(__ux, _mm256_add_pd(__one, _mm256_div_pd(_mm256_mul_pd(__dx, __alphax), __rho))));
                _mm256_storeu_pd(&_uy[i], _mm256_div_pd(__uy, _mm256_add_pd(__one, _mm256_div_pd(_mm256_mul_pd(__dx, __alphay), __rho))));
            }

            for (int i = ne; i < _particle.np; i++) {
                _ux[i] /= 1.0 + _particle.dx*_alphax[i]/_rho[i];
                _uy[i] /= 1.0 + _particle.dx*_alphay[i]/_rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Check convergence
        //*********************************************************************
        template<template<class>class P>
        bool CheckConvergence(P<double>& _particle, double _eps, double *_uxt, double *_uyt, double *_uxtm1, double *_uytm1) {
            assert(P<T>::nd == 2);
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            __m256d __unorm = { 0 }, __dunorm = { 0 };
            for (int i = 0; i < ne; i += packsize) {
                __mm256d __uxt = _mm256_loadu_pd(&_uxt[i]), __uyt = _mm256_loadu_pd(&_uyt[i]), __uxtm1 = _mm256_loadu_pd(&_uxtm1[i]), __uytm1 = _mm256_loadu_pd(&_uytm1[i]);
                __unorm = _mm256_add_pd(__unorm, _mm256_add_pd(_mm256_mul_pd(__uxt, __uxt), _mm256_mul_pd(__uyt, __uyt)));
                __mm256d __uxtmuxtm1 = _mm256_sub_pd(__uxt, __uxtm1), __uytmuytm1 = _mm256_sub_pd(__uyt, __uytm1); 
                __dunorm = _mm256_add_pd(__dunorm, _mm256_add_pd(_mm256_mul_pd(__uxtmuxtm1, __uxtmuxtm1), _mm256_mul_pd(__uytmuytm1, __uytmuytm1)));
            }

            double tmp[packsize] = { 0 }, tmpd[packsize] = { 0 };
            _mm256_store_pd(tmp, __unorm);
            _mm256_store_pd(tmpd, __dunorm);

            double unorm = std::accumulate(tmp, tmp + packsize, 0.0), dunorm = std::accumulate(tmpd, tmpd + packsize, 0.0);
            for (int i = ne; i < _particle.np; i++) {
                unorm += _uxt[i]*_uxt[i] + _uyt[i]*_uyt[i];
                dunorm += (_uxt[i] - _uxtm1[i])*(_uxt[i] - _uxtm1[i]) + (_uyt[i] - _uytm1[i])*(_uyt[i] - _uytm1[i]);
            }
            return sqrt(dunorm/unorm) < _eps ? true : false;
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Update macroscopic values, rho*, u*, v*
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _ip, T* _iux, T* _iuy, T* _imx, T* _imy) {
            assert(P<T>::nd == 2);
            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double a = 1.0, b = 3.0, c = 4.5, d = 1.5;
            __m256d __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b), __c = _mm256_broadcast_sd((const double*)&c), __d = _mm256_broadcast_sd((const double*)&d);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __rho = { 0 }, __ux = { 0 }, __uy = { 0 }, __ip = { 0 }, __iux = { 0 }, __iuy = { 0 }, __imx = { 0 }, __imy = { 0 };
                __m256d __uu = _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy));

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]), __ciu = _mm256_add_pd(_mm256_mul_pd(__cx[j], __ux), _mm256_mul_pd(__cy[j], __uy));
                    __ip = _mm256_add_pd(__ip, _mm256_mul_pd(__ft, _mm256_mul_pd(__ei[j], _mm256_add_pd(__a, _mm256_add_pd(_mm256_mul_pd(__b, __ciu), _mm256_suv_pd(_mm256_mul_pd(__c, _mm256_mul_pd(__ciu, __ciu)), _mm256_mul_pd(__d, __uu)))))));
                    __iux = _mm256_add_pd(__iux, _mm256_mul_pd(__ft, _mm256_mul_pd(__ei[j], _mm256_add_pd(__cx[j], _mm256_suv_pd(_mm256_mul_pd(b, _mm256_mul_pd(__ciu, __cx[j])), __ux)))));
                    __iuy = _mm256_add_pd(__iuy, _mm256_mul_pd(__ft, _mm256_mul_pd(__ei[j], _mm256_add_pd(__cy[j], _mm256_suv_pd(_mm256_mul_pd(b, _mm256_mul_pd(__ciu, __cy[j])), __uy)))));
                    __imx = _mm256_add_pd(__imx, _mm256_mul_pd(__ft, _mm256_mul_pd(__ei[j], __cx[j])));
                    __imy = _mm256_add_pd(__imy, _mm256_mul_pd(__ft, _mm256_mul_pd(__ei[j], __cy[j])));
                }
            }

            for (int i = ne; i < _particle.np; i++) {
                _ip[i] = T();
                _iux[i] = T();
                _iuy[i] = T();
                _imx[i] = T();
                _imy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i]; 
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    _ip[i] += _particle.ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                    _iux[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cx[j] + 3.0*ciu*P<T>::cx[j] - _ux[i]);
                    _iuy[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cy[j] + 3.0*ciu*P<T>::cy[j] - _uy[i]);
                    _imx[i] += _particle.ft[j][i]*P<T>::ei[j]*P<T>::cx[j];
                    _imy[i] += _particle.ft[j][i]*P<T>::ei[j]*P<T>::cy[j];
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _ux, T* _uy, T* _ip, T* _iux, T* _iuy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T feq = _ip[i] + 3.0*(_iux[i]*(P<T>::cx[j] - _ux[i]) + _iuy[i]*(P<T>::cy[j] - _uy[i]));
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _rho, T* _iux, T* _iuy, T* _imx, T* _imy, T* _alphax, T* _alphay) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _imx[i] /= 1.0 + _particle.dx*_alphax[i]/_rho[i];
                _imy[i] /= 1.0 + _particle.dx*_alphay[i]/_rho[i];
                _iux[i] -= _particle.dx*_alphax[i]*_imx[i]/_rho[i];
                _iuy[i] -= _particle.dx*_alphay[i]*_imy[i]/_rho[i];
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _ux, T _uy, T _ip, T _iux, T _iuy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                _particle.ft[j][_i] = _ip + 3.0*(_iux*(P<T>::cx[j] - _ux) + _iuy*(P<T>::cy[j] - _uy));
            }
        }
    }  
}