//*****************************************************************************
//  Title       :   src/equation_avx/advection_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection 2D    :   Update macroscopic values, T and flux
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _particle, double *_tem, double *_qx, double *_qy, double *_ux, double *_uy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            __m256d __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = { 0 }, __qx = { 0 }, __qy = { 0 }, __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]);

                for (int j = 0; j < P<double>::nc; j++) {
                    __mm256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __tem = _mm256_add_pd(__tem, __ft);
                    __qx = _mm256_add_pd(__qx, _mm256_mul_pd(__cx[j], __ft));
                    __qy = _mm256_add_pd(__qy, _mm256_mul_pd(__cy[j], __ft));
                }

                _mm256_storeu_pd(&_tem[i], __tem);
                _mm256_storeu_pd(&_qx[i], _mm256_sub_pd(__qx, _mm256_mul_pd(__tem, __ux)));
                _mm256_storeu_pd(&_qy[i], _mm256_sub_pd(__qy, _mm256_mul_pd(__tem, __uy)));
            }

            for (int i = ne; i < _particle.np; i++) {
                _tem[i] = 0.0;
                _qx[i] = 0.0;
                _qy[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _temperature[i] += _particle.ft[j][i];
                    _qx[i] += P<double>::cx[j]*_particle.ft[j][i];
                    _qy[i] += P<double>::cy[j]*_particle.ft[j][i];
                }
                _qx[i] -= _temperature[i]*_ux[i];
                _qy[i] -= _temperature[i]*_uy[i];
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Update macroscopic values, T and flux
        //*********************************************************************
        template<template<class>class P>
        void UpdateMacro(P<double>& _particle, double *_tem, double *_qx, double *_qy, double *_qz, double *_ux, double *_uy, double *_uz) {
            assert(P<double>::nd == 3);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            __m256d __cx[P<double>::nc], __cy[P<double>::nc], __cz[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double cxd = P<double>::cx[j], cyd = P<double>::cy[j], czd = P<double>::cz[j]; 
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
                __cz[j] = _mm256_broadcast_sd((const double*)&czd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = { 0 }, __qx = { 0 }, __qy = { 0 }, __qz = { 0 }, __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __uz = _mm256_loadu_pd(&_uz[i]);
                
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __tem = _mm256_add_pd(__tem, __ft);
                    __qx = _mm256_add_pd(__qx, _mm256_mul_pd(__cx[j], __ft));
                    __qy = _mm256_add_pd(__qy, _mm256_mul_pd(__cy[j], __ft));
                    __qz = _mm256_add_pd(__qz, _mm256_mul_pd(__cz[j], __ft));         
                }
                
                _mm256_storeu_pd(&_tem[i], __tem);
                _mm256_storeu_pd(&_qx[i], _mm256_sub_pd(__qx, _mm256_mul_pd(__tem, __ux)));
                _mm256_storeu_pd(&_qy[i], _mm256_sub_pd(__qy, _mm256_mul_pd(__tem, __uy)));
                _mm256_storeu_pd(&_qz[i], _mm256_sub_pd(__qz, _mm256_mul_pd(__tem, __uz)));
            }

            for (int i = ne; i < _particle.np; i++) {
                _temperature[i] = 0.0;
                _qx[i] = 0.0;
                _qy[i] = 0.0;
                _qz[i] = 0.0;
                for (int j = 0; j < P<double>::nc; j++) {
                    _temperature[i] += _particle.ft[j][i];
                    _qx[i] += P<double>::cx[j]*_particle.ft[j][i];
                    _qy[i] += P<double>::cy[j]*_particle.ft[j][i];
                    _qz[i] += P<double>::cz[j]*_particle.ft[j][i];
                }
                _qx[i] -= _temperature[i]*_ux[i];
                _qy[i] -= _temperature[i]*_uy[i];
                _qz[i] -= _temperature[i]*_uz[i];
            }
        }

        //*********************************************************************
        //  Advection 2D    :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _diffusivity, P<double>& _particle, double *_tem, double *_ux, double *_uy) {
            assert(P<double>::nd == 2);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5), iomega = 1.0 - omega, a = 1.0, b = 3.0;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega), __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = _mm256_loadu_pd(&_tem[i]), __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]);

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __m256d __ciu = _mm256_add_pd(_mm256_mul_pd(__cx[j], __ux), _mm256_mul_pd(__cy[j], __uy));
                    __m256d __feq = _mm256_mul_pd(__ei, _mm256_mul_pd(__tem, _mm256_add_pd(__a, _mm256_mul_pd(__b, __ciu))));
                    _mm256_storeu_pd(&_particle.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }
    
            for (int i = ne; i < _particle.np; i++) {
                for (int j = 0; j < P<double>::nc; j++) {
                    T ciu = P<double>::cx[j]*_ux[i] + P<double>::cy[j]*_uy[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<double>::ei[j]*_tem[i]*(1.0 + 3.0*ciu);
                }
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Collision term
        //*********************************************************************
        template<template<class>class P>
        void Collision(double _diffusivity, P<double>& _particle, double *_tem, double *_ux, double *_uy, double *_uz) {
            assert(P<double>::nd == 3);

            const int packsize = 32/sizeof(double), ne = (_particle.np/packsize)*packsize;

            double omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5), iomega = 1.0 - omega, a = 1.0, b = 3.0;
            __m256d __omega = _mm256_broadcast_sd((const double*)&omega), __iomega = _mm256_broadcast_sd((const double*)&iomega), __a = _mm256_broadcast_sd((const double*)&a), __b = _mm256_broadcast_sd((const double*)&b);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc], __cz[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j], czd = P<double>::cz[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
                __cz[j] = _mm256_broadcast_sd((const double*)&czd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = _mm256_loadu_pd(&_tem[i]), __ux = _mm256_loadu_pd(&_ux[i]), __uy = _mm256_loadu_pd(&_uy[i]), __uz = _mm256_loadu_pd(&_uz[i]);

                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particle.ft[j][i]);
                    __m256d __ciu = _mm256_add_pd(_mm256_mul_pd(__cx[j], __ux), _mm256_add_pd(_mm256_mul_pd(__cy[j], __uy), _mm256_mul_pd(__cz[j], __uz)));
                    __m256d __feq = _mm256_mul_pd(__ei, _mm256_mul_pd(__tem, _mm256_add_pd(__a, _mm256_mul_pd(__b, __ciu))));
                    _mm256_storeu_pd(&_particle.ftp1[j][i], _mm256_add_pd(_mm256_mul_pd(__iomega, __ft), _mm256_mul_pd(__omega, __feq)));
                }
            }

            for (int i = ne; i < _particle.np; i++) {
                for (int j = 0; j < P<double>::nc; j++) {
                    T ciu = P<double>::cx[j]*_ux[i] + P<double>::cy[j]*_uy[i] + P<double>::cz[j]*_uz[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<double>::ei[j]*_temperature[i]*(1.0 + 3.0*ciu);
                }
            }
        }

        //*********************************************************************
        //  Advection 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _temperature, T _ux, T _uy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy;
                _particle.ft[j][_i] = P<T>::ei[j]*_temperature*(1.0 + 3.0*ciu);
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _temperature, T _ux, T _uy, T _uz) {
            assert(P<T>::nd == 3 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy + P<T>::cz[j]*_uz;
                _particle.ft[j][_i] = P<T>::ei[j]*_temperature*(1.0 + 3.0*ciu);
            }
        }
    }

    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with natural convection
        //*********************************************************************
        template<template<class>class P, template<class>class Q>
        void ExternalForceNaturalConvection(double _gx, double _gy, double _tem0, P<double>& _particlef, Q<double>& _particleg) {
            assert(P<double>::nd == 2 && Q<double>::nd == 2 && _particlef.nx == _particleg.nx && _particlef.ny == _particleg.ny);

            const int packsize = 32/sizeof(double), ne = (_particlef.np/packsize)*packsize;

            double a = 3.0;
            __m256d __a = _mm256_broadcast_sd((const double*)&one), __dx = _mm256_broadcast_sd((const double*)&_particlef.dx);
            __m256d __gx = _mm256_broadcast_sd((const double*)&_gx), __gy = _mm256_broadcast_sd((const double*)&_gy), __tem0 = _mm256_broadcast_sd((const double*)&_tem0);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = { 0 };
                
                for (int j = 0; j < Q<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particleg.ft[j][i]);
                    __tem = _mm256_add_pd(__tem, __ft);
                }
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particlef.ft[j][i]), __cig = _mm256_add_pd(_mm256_mul_pd(__cx[j], __gx), _mm256_mul_pd(__cy[j], __gy));
                    _mm256_storeu_pd(&_particlef.ft[j][i], _mm256_add_pd(__ft, _mm256_mul_pd(__a, _mm256_mul_pd(__dx, _mm256_mul_pd(__ei[j], _mm256_mul_pd(__cig, _mm256_sub_pd(__tem, __tem0)))))));
                }
            }

            for (int i = ne; i < _particlef.np; i++) {
                double temperature = 0.0;
                for (int j = 0; j < Q<double>::nc; j++) {
                    temperature += _particleg.ft[j][i];
                }
                for (int j = 0; j < P<double>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_particlef.dx*P<double>::ei[j]*(Q<double>::cx[j]*_gx + Q<double>::cy[j]*_gy)*(temperature - _temperature0);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   External force with natural convection
        //*********************************************************************
        template<template<class>class P, template<class>class Q>
        void ExternalForceNaturalConvection(double _gx, double _gy, double _gz, double _tem0, P<double>& _particlef, Q<double>& _particleg) {
            assert(P<double>::nd == 3 && Q<double>::nd == 3 && _particlef.nx == _particleg.nx && _particlef.ny == _particleg.ny && _particlef.nz == _particleg.nz);
            
            const int packsize = 32/sizeof(double), ne = (_particlef.np/packsize)*packsize;

            double a = 3.0;
            __m256d __a = _mm256_broadcast_sd((const double*)&one), __dx = _mm256_broadcast_sd((const double*)&_particlef.dx);
            __m256d __gx = _mm256_broadcast_sd((const double*)&_gx), __gy = _mm256_broadcast_sd((const double*)&_gy), __gz = _mm256_broadcast_sd((const double*)&_gz), __tem0 = _mm256_broadcast_sd((const double*)&_tem0);

            __m256d __ei[P<double>::nc], __cx[P<double>::nc], __cy[P<double>::nc], __cz[P<double>::nc];
            for (int j = 0; j < P<double>::nc; j++) {
                double eid = P<double>::ei[j], cxd = P<double>::cx[j], cyd = P<double>::cy[j], czd = P<double>::cz[j]; 
                __ei[j] = _mm256_broadcast_sd((const double*)&eid);
                __cx[j] = _mm256_broadcast_sd((const double*)&cxd);
                __cy[j] = _mm256_broadcast_sd((const double*)&cyd);
                __cz[j] = _mm256_broadcast_sd((const double*)&czd);
            }

            for (int i = 0; i < ne; i += packsize) {
                __m256d __tem = { 0 };
                
                for (int j = 0; j < Q<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particleg.ft[j][i]);
                    __tem = _mm256_add_pd(__tem, __ft);
                }
                for (int j = 0; j < P<double>::nc; j++) {
                    __m256d __ft = _mm256_loadu_pd(&_particlef.ft[j][i]), __cig = _mm256_add_pd(_mm256_mul_pd(__cx[j], __gx), _mm256_add_pd(_mm256_mul_pd(__cy[j], __gy), _mm256_mul_pd(__cz[j], __gz)));
                    _mm256_storeu_pd(&_particlef.ft[j][i], _mm256_add_pd(__ft, _mm256_mul_pd(__a, _mm256_mul_pd(__dx, _mm256_mul_pd(__ei[j], _mm256_mul_pd(__cig, _mm256_sub_pd(__tem, __tem0)))))));
                }
            }
            
            for (int i = ne; i < _particlef.np; i++) {
                double temperature = 0.0;
                for (int j = 0; j < Q<double>::nc; j++) {
                    temperature += _particleg.ft[j][i];
                }
                for (int j = 0; j < P<double>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_particlef.dx*P<double>::ei[j]*(Q<double>::cx[j]*_gx + Q<double>::cy[j]*_gy + Q<double>::cz[j]*_gz)*(temperature - _temperature0);
                }
            }
        }
    }
}