//*****************************************************************************
//  Title       :   src/equation_avx/adjointadvection_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace AAD {
        //  Function of getting sensitivity of temperature at heat source for D2Q9 along x edge
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongXEdge(
            Q<double>& _q, int _i, int _directionx, const double *_ux, const double *_uy, const double *_ig, double *_dfds, const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        double qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directionx == -1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 + 3.0*_ux[idx])*(-6.0 + 4.0*_ig[IndexG(idx, 1)] + _ig[IndexG(idx, 5)] + _ig[IndexG(idx, 8)])
                                + 3.0*_uy[idx]*(_ig[IndexG(idx, 5)] - _ig[IndexG(idx, 8)])
                            )/(36.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        } else if (_directionx == 1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 - 3.0*_ux[idx])*(-6.0 + 4.0*_ig[IndexG(idx, 3)] + _ig[IndexG(idx, 6)] + _ig[IndexG(idx, 7)])
                                + 3.0*_uy[idx]*(_ig[IndexG(idx, 6)] - _ig[IndexG(idx, 7)])
                            )/(36.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D2Q9 along y edge
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongYEdge(
            Q<double>& _q, int _j, int _directiony, const double *_ux, const double *_uy, const double *_ig, double *_dfds, const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        double qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directiony == -1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 + 3.0*_uy[idx])*(-6.0 + 4.0*_ig[IndexG(idx, 2)] + _ig[IndexG(idx, 5)] + _ig[IndexG(idx, 6)])
                                + 3.0*_ux[idx]*(_ig[IndexG(idx, 5)] - _ig[IndexG(idx, 6)])
                            )/(36.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        } else if (_directiony == 1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 - 3.0*_uy[idx])*(-6.0 + 4.0*_ig[IndexG(idx, 4)] + _ig[IndexG(idx, 7)] + _ig[IndexG(idx, 8)])
                                + 3.0*_ux[idx]*(_ig[IndexG(idx, 8)] - _ig[IndexG(idx, 7)])
                            )/(36.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15 along x face
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongXFace(
            Q<double>& _q, int _i, int _directionx, const double *_ux, const double *_uy, const double *_uz, const double *_ig, double *_dfds, const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            double qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionx == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_ux[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 1)] + _ig[IndexG(idx, 7)] + _ig[IndexG(idx, 9)] + _ig[IndexG(idx, 10)] + _ig[IndexG(idx, 12)])
                                    + 3.0*_uy[idx]*(_ig[IndexG(idx, 7)] - _ig[IndexG(idx, 9)] + _ig[IndexG(idx, 10)] - _ig[IndexG(idx, 12)])
                                    + 3.0*_uz[idx]*(_ig[IndexG(idx, 7)] + _ig[IndexG(idx, 9)] - _ig[IndexG(idx, 10)] - _ig[IndexG(idx, 12)])
                                )/(72.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directionx == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_ux[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 4)] + _ig[IndexG(idx, 8)] + _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 13)] + _ig[IndexG(idx, 14)])
                                    + 3.0*_uy[idx]*(_ig[IndexG(idx, 8)] - _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 13)] - _ig[IndexG(idx, 14)])
                                    + 3.0*_uz[idx]*(_ig[IndexG(idx, 8)] - _ig[IndexG(idx, 11)] - _ig[IndexG(idx, 13)] + _ig[IndexG(idx, 14)])
                                )/(72.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15 along y face
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongYFace(
            Q<double>& _q, int _j, int _directiony, const double *_ux, const double *_uy, const double *_uz, const double *_ig, double *_dfds, const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            double qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directiony == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_uy[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 2)] + _ig[IndexG(idx, 7)] + _ig[IndexG(idx, 8)] + _ig[IndexG(idx, 10)] + _ig[IndexG(idx, 13)])
                                    + 3.0*_uz[idx]*(_ig[IndexG(idx, 7)] + _ig[IndexG(idx, 8)] - _ig[IndexG(idx, 10)] - _ig[IndexG(idx, 13)])
                                    + 3.0*_ux[idx]*(_ig[IndexG(idx, 7)] - _ig[IndexG(idx, 8)] + _ig[IndexG(idx, 10)] - _ig[IndexG(idx, 13)])
                                )/(72.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directiony == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_uy[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 5)] + _ig[IndexG(idx, 9)] + _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 12)] + _ig[IndexG(idx, 14)])
                                    + _uz[idx]*(_ig[IndexG(idx, 9)] - _ig[IndexG(idx, 11)] - _ig[IndexG(idx, 12)] + _ig[IndexG(idx, 14)])
                                    + _ux[idx]*(_ig[IndexG(idx, 9)] - _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 12)] - _ig[IndexG(idx, 14)])
                                )/(72.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }
            
        //  Function of getting sensitivity of temperature at heat source for D3Q15 along z face
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongZFace(
            Q<double>& _q, int _k, int _directionz, const double *_ux, const double *_uy, const double *_uz, const double *_ig, double *_dfds, const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            double qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionz == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_uz[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 3)] + _ig[IndexG(idx, 7)] + _ig[IndexG(idx, 8)] + _ig[IndexG(idx, 9)] + _ig[IndexG(idx, 14)])
                                    + _ux[idx]*(_ig[IndexG(idx, 7)] - _ig[IndexG(idx, 8)] + _ig[IndexG(idx, 9)] - _ig[IndexG(idx, 14)])
                                    + _uy[idx]*(_ig[IndexG(idx, 7)] + _ig[IndexG(idx, 8)] - _ig[IndexG(idx, 9)] - _ig[IndexG(idx, 14)])
                                )/(72.0*(1.0 - 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directionz == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_uz[idx])*(-12.0 + 8.0*_ig[IndexG(idx, 6)] + _ig[IndexG(idx, 10)] + _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 12)] + _ig[IndexG(idx, 13)])
                                    + _ux[idx]*(_ig[IndexG(idx, 10)] - _ig[IndexG(idx, 11)] + _ig[IndexG(idx, 12)] - _ig[IndexG(idx, 13)])
                                    + _uy[idx]*(_ig[IndexG(idx, 10)] - _ig[IndexG(idx, 11)] - _ig[IndexG(idx, 12)] + _ig[IndexG(idx, 13)])
                                )/(72.0*(1.0 + 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }   
    }

    namespace AAD {
        template<class T, template<class>class Q>void Macro(T &, T &, T &, const T *, const T *, int);                              //  Function of updating macroscopic values of AAD for 2D
        template<class T, template<class>class Q>void Macro(T &, T &, T &, T &, const T *, const T *, int);                         //  Function of updating macroscopic values of AAD for 3D
        template<class T, template<class>class Q>void Equilibrium(T *, T, T, T, T, T);                                              //  Function of getting equilibrium of AAD for 2D
        template<class T, template<class>class Q>void Equilibrium(T *, T, T, T, T, T, T, T);                                        //  Function of getting equilibrium of AAD for 3D
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T, T, T, T, T, T *, T *, T, int);           //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T, T, T, T, T, T, T, T, T *, T *, T, int);  //  Function of applying external force with Brinkman model and advection of AAD for 3D
        template<class T, template<class>class Q>void ExternalForceHeatExchange(T, T *, T *, T, int);                               //  Function of applying external force with heat exchange of AAD for 2D/3D       
        template<class T, template<class>class Q>void ExternalForceNaturalConvection(T, T, T, T, T *, T *, int);                    //  Function of applying external force with natural convection of AAD for 2D
        template<class T, template<class>class Q>void ExternalForceNaturalConvection(T, T, T, T, T, T, T *, T *, int);              //  Function of applying external force with natural convection of AAD for 3D
        template<class T, template<class>class P>void ExternalForceMassFlow(T, T, T, T, T, T *, T *, int);                          //  Function of applying external force with mass flow maximization of AAD for 2D
        template<class T, template<class>class P>void ExternalForceMassFlow(T, T, T, T, T, T, T, T *, T *, int);                    //  Function of applying external force with mass flow maximization of AAD for 3D
        template<class T, template<class>class Q>void ExternalForceHeatCompliance(T, T *, T *, int);                                //  Function of applying external force with heat compliance of AAD for 2D/3D

        //  Function of updating macroscopic values of AAD for 2D
        template<class Q>
        void Macro(__m256d &__item, __m256d &__iqx, __m256d &__iqy, const __m256d *__g) {
            __item = _mm256_mul_pd(Q::__ei[0], __g[0]);
            __iqx = _mm256_setzero_pd();
            __iqy = _mm256_setzero_pd();
            for (int c = 1; c < Q::nc; ++c) {
                __m256d __gei = _mm256_mul_pd(Q::__ei[c], __g[c]);
                __item = _mm256_add_pd(__item, __gei);
                __iqx = _mm256_add_pd(__iqx, _mm256_mul_pd(__gei, Q::__cx[c]));
                __iqy = _mm256_add_pd(__iqy, _mm256_mul_pd(__gei, Q::__cy[c]));
            }
        }

        //  Function of updating macroscopic values of AAD for 3D
        template<class Q>
        void Macro(__m256d &__item, __m256d &__iqx, __m256d &__iqy, __m256d &__iqz, const __m256d *__g) {
            __item = _mm256_mul_pd(Q::__ei[0], __g[0]);
            __iqx = _mm256_setzero_pd();
            __iqy = _mm256_setzero_pd();
            __iqz = _mm256_setzero_pd();
            for (int c = 1; c < Q::nc; ++c) {
                __m256d __gei = _mm256_mul_pd(Q::__ei[c], __g[c]);
                __item = _mm256_add_pd(__item, __gei);
                __iqx = _mm256_add_pd(__iqx, _mm256_mul_pd(__gei, Q::__cx[c]));
                __iqy = _mm256_add_pd(__iqy, _mm256_mul_pd(__gei, Q::__cy[c]));
                __iqz = _mm256_add_pd(__iqz, _mm256_mul_pd(__gei, Q::__cz[c]));
            }
        }

        //  Function of getting equilibrium of AAD for 2D
        template<class Q>
        void Equilibrium(__m256d *__geq, const __m256d &__item, const __m256d &__iqx, const __m256d &__iqy, const __m256d &__ux, const __m256d &__uy) {
            __m256d __coef = _mm256_add_pd(__item, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_mul_pd(__iqx, __ux), _mm256_mul_pd(__iqy, __uy))));
            for (int c = 0; c < Q::nc; ++c) {
                __geq[c] = __coef;
            }
        }

        //  Function of getting equilibrium of AAD for 3D
        template<class Q>
        void Equilibrium(__m256d *__geq, const __m256d &__item, const __m256d &__iqx, const __m256d &__iqy, const __m256d &__iqz, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz) {
            __m256d __coef = _mm256_add_pd(__item, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__iqx, __ux), _mm256_mul_pd(__iqy, __uy)), _mm256_mul_pd(__iqz, __uz))));
            for (int c = 0; c < Q::nc; ++c) {
                __geq[c] = __coef;
            } 
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__imx, const __m256d &__imy, const __m256d &__tem, const __m256d &__iqx, const __m256d &__iqy, const __m256d &__omegag, __m256d *__f, const __m256d &__alpha) {
            __m256d __coef = _mm256_div_pd(_mm256_set1_pd(3.0), _mm256_add_pd(__rho, __alpha));
            __m256d __coefx = _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(__tem, __iqx), __omegag), _mm256_mul_pd(__alpha, __imx));
            __m256d __coefy = _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(__tem, __iqy), __omegag), _mm256_mul_pd(__alpha, __imy));
            __f[0] = _mm256_sub_pd(__f[0], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_mul_pd(__coefx, __ux), _mm256_mul_pd(__coefy, __uy))));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_add_pd(__f[c], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_mul_pd(__coefx, _mm256_sub_pd(P::__cx[c], __ux)), _mm256_mul_pd(__coefy, _mm256_sub_pd(P::__cy[c], __uy)))));
            }
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 3D
        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d &__imx, const __m256d &__imy, const __m256d &__imz, const __m256d &__tem, const __m256d &__iqx, const __m256d &__iqy, const __m256d &__iqz, const __m256d &__omegag, __m256d *__f, const __m256d &__alpha) {
            __m256d __coef = _mm256_div_pd(_mm256_set1_pd(3.0), _mm256_add_pd(__rho, __alpha));
            __m256d __coefx = _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(__tem, __iqx), __omegag), _mm256_mul_pd(__alpha, __imx));
            __m256d __coefy = _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(__tem, __iqy), __omegag), _mm256_mul_pd(__alpha, __imy));
            __m256d __coefz = _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(__tem, __iqz), __omegag), _mm256_mul_pd(__alpha, __imz));
            __f[0] = _mm256_sub_pd(__f[0], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__coefx, __ux), _mm256_mul_pd(__coefy, __uy)), _mm256_mul_pd(__coefz, __uz))));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_add_pd(__f[c], _mm256_mul_pd(__coef, 
                    _mm256_add_pd(
                        _mm256_add_pd(
                            _mm256_mul_pd(__coefx, _mm256_sub_pd(P::__cx[c], __ux)), 
                            _mm256_mul_pd(__coefy, _mm256_sub_pd(P::__cy[c], __uy))
                        ), _mm256_mul_pd(__coefz, _mm256_sub_pd(P::__cz[c], __uz))
                    )
                ));
            }
        }

        //  Function of applying external force with heat exchange of AAD for 2D/3D
        template<class Q>
        void ExternalForceHeatExchange(const __m256d &__item, __m256d *__g, const __m256d &__beta) {
            __m256d __coef = _mm256_mul_pd(__beta, _mm256_div_pd(_mm256_add_pd(_mm256_set1_pd(1.0), __item), _mm256_add_pd(_mm256_set1_pd(1.0), __beta)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_sub_pd(__g[c], __coef);
            }
        }

        //  Function of applying external force with natural convection of AAD for 2D
        template<class Q>
        void ExternalForceNaturalConvection(const __m256d &__imx, const __m256d &__imy, const __m256d &__gx, const __m256d &__gy, __m256d *__g) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_mul_pd(__imx, __gx), _mm256_mul_pd(__imy, __gy)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_add_pd(__g[c], __coef);
            }
        }

        //  Function of applying external force with natural convection of AAD for 3D
        template<class Q>
        void ExternalForceNaturalConvection(const __m256d &__imx, const __m256d &__imy, const __m256d &__imz, const __m256d &__gx, const __m256d &__gy, const __m256d &__gz, __m256d *__g) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__imx, __gx), _mm256_mul_pd(__imy, __gy)), _mm256_mul_pd(__imz, __gz)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_add_pd(__g[c], __coef);
            }
        }

        //  Function of applying external force with mass flow maximization of AAD for 2D
        template<class P>
        void ExternalForceMassFlow(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__directionx, const __m256d &__directiony, __m256d *__f) {
            for (int c = 0; c < P::nc; ++c) {
                __f[c]=_mm256_sub_pd(__f[c], _mm256_div_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(P::__cx[c], __ux), __directionx), _mm256_mul_pd(_mm256_sub_pd(P::__cy[c], __uy), __directiony)), __rho));
            }
        }

        //  Function of applying external force with mass flow maximization of AAD for 3D
        template<class P>
        void ExternalForceMassFlow(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d &__directionx, const __m256d &__directiony, const __m256d &__directionz, __m256d *__f) {
            for (int c = 0; c < P::nc; ++c) {
                __f[c]=_mm256_sub_pd(__f[c], _mm256_div_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(P::__cx[c], __ux), __directionx), _mm256_mul_pd(_mm256_sub_pd(P::__cy[c], __uy), __directiony)), _mm256_mul_pd(_mm256_sub_pd(P::__cz[c], __uz), __directionz)), __rho));
            }
        }

        //  Function of applying external force with heat compliance of AAD for 2D/3D
        template<class Q>
        void ExternalForceHeatCompliance(const __m256d &__heatsource, __m256d *__g) {
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_sub_pd(__g[c], __heatsource);
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange) and Collide of AAD for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_beta, double _diffusivity, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc];
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag), __geq[Q<double>::nc];
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                __m256d __item, __iqx, __iqy;
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __imx, __imy, __tem, __iqx, __iqy, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                __m256d __beta = _mm256_loadu_pd(&_beta[idx]);
                ExternalForceHeatExchange<Q<double> >(__item, __g, __beta);
                Macro<Q<double> >(__item, __iqx, __iqy, __g);  

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __ip, __iux, __iuy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __ux, __uy);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, imx, imy;
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy;
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                ExternalForceHeatExchange<double, Q>(item, _q.f0, _q.f, _beta[idx], idx);
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c); 
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange) and Collide of AAD for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_beta, double _diffusivity, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc];
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag), __geq[Q<double>::nc];
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __item, __iqx, __iqy, __iqz;
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __tem, __iqx, __iqy, __iqz, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __beta = _mm256_loadu_pd(&_beta[idx]);
                ExternalForceHeatExchange<Q<double> >(__item, __g, __beta);
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);  

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                    _mm256_storeu_pd(&_iqz[idx], __iqz);
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __iqz, __ux, __uy, __uz);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy, iqz;
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceHeatExchange<double, Q>(item, _q.f0, _q.f, _beta[idx], idx);
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of AAD for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_diffusivity, bool _issave = false, double *_g = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                __m256d __item, __iqx, __iqy;
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __imx, __imy, __tem, __iqx, __iqy, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_g[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __ip, __iux, __iuy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __ux, __uy);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, imx, imy;
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy;
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, bool _issave = false, double *_g = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __item, __iqx, __iqy, __iqz;
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __tem, __iqx, __iqy, __iqz, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                    _mm256_storeu_pd(&_iqz[idx], __iqz);

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_g[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __iqz, __ux, __uy, __uz);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy, iqz;
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_diffusivity, double _gx, double _gy, bool _issave = false, double *_g = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy);
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                __m256d __item, __iqx, __iqy;
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __imx, __imy, __tem, __iqx, __iqy, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                ExternalForceNaturalConvection<Q<double> >(__imx, __imy, __gx, __gy, __g);
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_g[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __ip, __iux, __iuy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __ux, __uy);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, imx, imy;
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy;
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<double, Q>(imx, imy, _gx, _gy, _q.f0, _q.f, idx);
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, double _gx, double _gy, double _gz, bool _issave = false, double *_g = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __gz = _mm256_set1_pd(_gz);
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __item, __iqx, __iqy, __iqz;
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __tem, __iqx, __iqy, __iqz, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                ExternalForceNaturalConvection<Q<double> >(__imx, __imy, __imz, __gx, __gy, __gz, __g);
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                    _mm256_storeu_pd(&_iqz[idx], __iqz);

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_g[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __iqz, __ux, __uy, __uz);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy, iqz;
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<double, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD mass flow for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvectionMassFlow(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_diffusivity, double _gx, double _gy, const double *_directionx, const double *_directiony, bool _issave = false, double *_ig = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy);
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                __m256d __item, __iqx, __iqy;
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  External force with Brinkman model and mass flow
                __m256d __directionx = _mm256_loadu_pd(&_directionx[idx]), __directiony = _mm256_loadu_pd(&_directiony[idx]);
                ExternalForceMassFlow<P<double> >(__rho, __ux, __uy, __directionx, __directiony, __f);
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __imx, __imy, __tem, __iqx, __iqy, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);
                ExternalForceNaturalConvection<Q<double> >(__imx, __imy, __gx, __gy, __g);
                Macro<Q<double> >(__item, __iqx, __iqy, __g);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);

                    if (_ig) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_ig[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __ip, __iux, __iuy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __ux, __uy);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, imx, imy;
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy;
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model and mass flow
                ExternalForceMassFlow<double, P>(_rho[idx], _ux[idx], _uy[idx], _directionx[idx], _directiony[idx], _p.f0, _p.f, idx);
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<double, Q>(imx, imy, _gx, _gy, _q.f0, _q.f, idx);
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;

                    if (_ig) {
                        int offsetf = Q<double>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD mass flow for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvectionMassFlow(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, double _gx, double _gy, double _gz, const double *_directionx, const double *_directiony, const double *_directionz, bool _issave = false, double *_ig = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __gz = _mm256_set1_pd(_gz);
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __item, __iqx, __iqy, __iqz;
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  External force with Brinkman model and mass flow
                __m256d __directionx = _mm256_loadu_pd(&_directionx[idx]), __directiony = _mm256_loadu_pd(&_directiony[idx]), __directionz = _mm256_loadu_pd(&_directionz[idx]);
                ExternalForceMassFlow<P<double> >(__rho, __ux, __uy, __directionx, __directiony, __directionz, __f);
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __tem, __iqx, __iqy, __iqz, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                ExternalForceNaturalConvection<Q<double> >(__imx, __imy, __imz, __gx, __gy, __gz, __g);
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                    _mm256_storeu_pd(&_iqz[idx], __iqz);

                    if (_ig) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_ig[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __iqz, __ux, __uy, __uz);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy, iqz;
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model and mass flow
                ExternalForceMassFlow<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], _directionx[idx], _directiony[idx], _directionz[idx], _p.f0, _p.f, idx);
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<double, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_ig) {
                        int offsetf = Q<double>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvectionHeatCompliance(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, const double *_heatsource, double _gx, double _gy, double _gz, bool _issave = false, double *_g = nullptr
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf), __feq[P<double>::nc], __geq[Q<double>::nc];
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __gz = _mm256_set1_pd(_gz);
            #pragma omp parallel for private(__feq, __geq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                _p.LoadF(idx, __f);
                _q.LoadF(idx, __g);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                __m256d __item, __iqx, __iqy, __iqz;
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __tem, __iqx, __iqy, __iqz, __omegag, __f, __alpha);
                ANS::Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);
                ExternalForceNaturalConvection<Q<double> >(__imx, __imy, __imz, __gx, __gy, __gz, __g);
                __m256d __heatsource = _mm256_loadu_pd(&_heatsource[idx]);
                ExternalForceHeatCompliance<Q<double> >(__heatsource, __g);
                Macro<Q<double> >(__item, __iqx, __iqy, __iqz, __g);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                    _mm256_storeu_pd(&_item[idx], __item);
                    _mm256_storeu_pd(&_iqx[idx], __iqx);
                    _mm256_storeu_pd(&_iqy[idx], __iqy);
                    _mm256_storeu_pd(&_iqz[idx], __iqz);

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        for (int c = 0; c < Q<double>::nc; ++c) {
                            _mm256_storeu_pd(&_g[offsetf + Q<double>::packsize*c], __g[c]);
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, __feq[c]));
                }
                _p.StoreF(idx, __f);
                Equilibrium<Q<double> >(__geq, __item, __iqx, __iqy, __iqz, __ux, __uy, __uz);
                for (int c = 0; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, __geq[c]));
                }
                _q.StoreF(idx, __g);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy, iqz;
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<double, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                ExternalForceHeatCompliance<double, Q>(_heatsource[idx], _q.f0, _q.f, idx);
                Macro<double, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_g) {
                        int offsetf = Q<double>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<double>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<double>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of getting sensitivity of heat exchange
        template<template<class>class Q>
        void SensitivityHeatExchange(
            Q<double>& _q, double *_dfds, 
            const double *_ux, const double *_uy, const double *_imx, const double *_imy, const double *_dads, 
            const double *_tem, const double *_item, const double *_dbds
        ) {
            const int ne = _q.nxyz/Q<double>::packsize;
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*Q<double>::packsize;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dbds = _mm256_loadu_pd(&_dbds[idx]);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                __m256d __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]), __item = _mm256_loadu_pd(&_item[idx]);
                __dfds=_mm256_add_pd(__dfds, _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __dads), _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy))), _mm256_mul_pd(_mm256_mul_pd(__dbds, _mm256_sub_pd(_mm256_set1_pd(1.0), __tem)), _mm256_add_pd(_mm256_set1_pd(1.0), __item))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*Q<double>::packsize; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]) - _dbds[idx]*(1.0 - _tem[idx])*(1.0 + _item[idx]);
            }
        }

        //  Function of getting sensitivity of heat exchange
        template<template<class>class Q>
        void SensitivityHeatExchange(
            Q<double>& _q, double *_dfds, 
            const double *_ux, const double *_uy, const double *_uz, const double *_imx, const double *_imy, const double *_imz, const double *_dads, 
            const double *_tem, const double *_item, const double *_dbds
        ) {
            const int ne = _q.nxyz/Q<double>::packsize;
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*Q<double>::packsize;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dbds = _mm256_loadu_pd(&_dbds[idx]);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]), __tem = _mm256_loadu_pd(&_tem[idx]);
                __m256d __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]), __imz = _mm256_loadu_pd(&_imz[idx]), __item = _mm256_loadu_pd(&_item[idx]);
                __dfds=_mm256_add_pd(__dfds, _mm256_sub_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __dads), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy)), _mm256_mul_pd(__uz, __imz))), _mm256_mul_pd(_mm256_mul_pd(__dbds, _mm256_sub_pd(_mm256_set1_pd(1.0), __tem)), _mm256_add_pd(_mm256_set1_pd(1.0), __item))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*Q<double>::packsize; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]) - _dbds[idx]*(1.0 - _tem[idx])*(1.0 + _item[idx]);
            }
        }

        //  Function of getting sensitivity of Brinkman model and diffusivity term
        template<template<class>class Q>
        void SensitivityBrinkmanDiffusivity(
            Q<double>& _q, double *_dfds, const double *_ux, const double *_uy, const double *_imx, const double *_imy, const double *_dads,
            const double *_tem, const double *_item, const double *_iqx, const double *_iqy, const double *_g, const double *_ig,
            const double *_diffusivity, const double *_dkds
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            //  Brinkman term and diffusivity term
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;

                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dkds = _mm256_loadu_pd(&_dkds[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]);
                __m256d __tem = _mm256_loadu_pd(&_tem[idx]), __item = _mm256_loadu_pd(&_item[idx]), __iqx = _mm256_loadu_pd(&_iqx[idx]), __iqy = _mm256_loadu_pd(&_iqy[idx]);

                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy)))));

                int offsetf = nc*idx;
                __m256d __sumg = _mm256_setzero_pd(), __taug = _mm256_add_pd(_mm256_mul_pd(__3, _mm256_loadu_pd(&_diffusivity[idx])), _mm256_set1_pd(0.5));
                for (int c = 0; c < nc; ++c) {
                    int idxf = offsetf + ps*c;
                    __m256d __g = _mm256_loadu_pd(&_g[idxf]), __ig = _mm256_loadu_pd(&_ig[idxf]);
                    __sumg = _mm256_add_pd(__sumg, _mm256_mul_pd(__g, __ig));
                }
                __dfds = _mm256_sub_pd(__dfds, _mm256_div_pd(
                    _mm256_mul_pd(__3, _mm256_mul_pd(__dkds, _mm256_sub_pd(
                        __sumg, _mm256_mul_pd(__tem, _mm256_add_pd(__item, _mm256_mul_pd(__3, _mm256_add_pd(_mm256_mul_pd(__ux, __iqx), _mm256_mul_pd(__uy, __iqy)))))
                    ))), _mm256_mul_pd(__taug, __taug)
                ));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _q.nxyz; ++idx) {
                
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                
                int offsetf = nc*idx;
                double sumg = 0.0;
                for (int c = 0; c < nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
            }
        }

        //  Function of getting sensitivity of Brinkman model and diffusivity term
        template<template<class>class Q>
        void SensitivityBrinkmanDiffusivity(
            Q<double>& _q, double *_dfds, const double *_ux, const double *_uy, const double *_uz, const double *_imx, const double *_imy, const double *_imz, const double *_dads,
            const double *_tem, const double *_item, const double *_iqx, const double *_iqy, const double *_iqz, const double *_g, const double *_ig,
            const double *_diffusivity, const double *_dkds
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            //  Brinkman term and diffusivity term
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;

                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dkds = _mm256_loadu_pd(&_dkds[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                __m256d __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]), __imz = _mm256_loadu_pd(&_imz[idx]);
                __m256d __tem = _mm256_loadu_pd(&_tem[idx]), __item = _mm256_loadu_pd(&_item[idx]), __iqx = _mm256_loadu_pd(&_iqx[idx]), __iqy = _mm256_loadu_pd(&_iqy[idx]), __iqz = _mm256_loadu_pd(&_iqz[idx]);

                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_add_pd(_mm256_mul_pd(__uy, __imy), _mm256_mul_pd(__uz, __imz))))));

                int offsetf = nc*idx;
                __m256d __sumg = _mm256_setzero_pd(), __taug = _mm256_add_pd(_mm256_mul_pd(__3, _mm256_loadu_pd(&_diffusivity[idx])), _mm256_set1_pd(0.5));
                for (int c = 0; c < nc; ++c) {
                    int idxf = offsetf + ps*c;
                    __m256d __g = _mm256_loadu_pd(&_g[idxf]), __ig = _mm256_loadu_pd(&_ig[idxf]);
                    __sumg = _mm256_add_pd(__sumg, _mm256_mul_pd(__g, __ig));
                }
                __dfds = _mm256_sub_pd(__dfds, _mm256_div_pd(
                    _mm256_mul_pd(__3, _mm256_mul_pd(__dkds, _mm256_sub_pd(
                        __sumg, _mm256_mul_pd(__tem, _mm256_add_pd(__item, _mm256_mul_pd(
                            __3, _mm256_add_pd(_mm256_mul_pd(__ux, __iqx), _mm256_add_pd(_mm256_mul_pd(__uy, __iqy), _mm256_mul_pd(__uz, __iqz)))
                        )))
                    ))), _mm256_mul_pd(__taug, __taug)
                ));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                int offsetf = nc*idx;
                double sumg = 0.0;
                for (int c = 0; c < nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D2Q9
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            Q<double>& _q, double *_dfds, const double *_ux, const double *_uy, const double *_imx, const double *_imy, const double *_dads,
            const double *_tem, const double *_item, const double *_iqx, const double *_iqy, const double *_g, const double *_ig,
            const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            //  Brinkman term and diffusivity term
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;

                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dkds = _mm256_loadu_pd(&_dkds[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]);
                __m256d __tem = _mm256_loadu_pd(&_tem[idx]), __item = _mm256_loadu_pd(&_item[idx]), __iqx = _mm256_loadu_pd(&_iqx[idx]), __iqy = _mm256_loadu_pd(&_iqy[idx]);

                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy)))));

                int offsetf = nc*idx;
                __m256d __sumg = _mm256_setzero_pd(), __taug = _mm256_add_pd(_mm256_mul_pd(__3, _mm256_loadu_pd(&_diffusivity[idx])), _mm256_set1_pd(0.5));
                for (int c = 0; c < nc; ++c) {
                    int idxf = offsetf + ps*c;
                    __m256d __g = _mm256_loadu_pd(&_g[idxf]), __ig = _mm256_loadu_pd(&_ig[idxf]);
                    __sumg = _mm256_add_pd(__sumg, _mm256_mul_pd(__g, __ig));
                }
                __dfds = _mm256_sub_pd(__dfds, _mm256_div_pd(
                    _mm256_mul_pd(__3, _mm256_mul_pd(__dkds, _mm256_sub_pd(
                        __sumg, _mm256_mul_pd(__tem, _mm256_add_pd(__item, _mm256_mul_pd(__3, _mm256_add_pd(_mm256_mul_pd(__ux, __iqx), _mm256_mul_pd(__uy, __iqy)))))
                    ))), _mm256_mul_pd(__taug, __taug)
                ));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _q.nxyz; ++idx) {
                
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                
                int offsetf = nc*idx;
                double sumg = 0.0;
                for (int c = 0; c < nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
            }

            SensitivityTemperatureAtHeatSourceAlongXEdge(_q, 0, -1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);         //  Boundary term along xmin
            SensitivityTemperatureAtHeatSourceAlongXEdge(_q, _q.lx - 1, 1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);  //  Boundary term along xmax
            SensitivityTemperatureAtHeatSourceAlongYEdge(_q, 0, -1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);         //  Boundary term along ymin
            SensitivityTemperatureAtHeatSourceAlongYEdge(_q, _q.ly - 1, 1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);  //  Boundary term along ymax
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15
        template<template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            Q<double>& _q, double *_dfds, const double *_ux, const double *_uy, const double *_uz, const double *_imx, const double *_imy, const double *_imz, const double *_dads,
            const double *_tem, const double *_item, const double *_iqx, const double *_iqy, const double *_iqz, const double *_g, const double *_ig,
            const double *_diffusivity, const double *_dkds, Fv _qnbc, Ff _bctype
        ) {
            const int ps = Q<double>::packsize, ne = _q.nxyz/ps, nc = Q<double>::nc;
            auto IndexG = [=](int _idx, int _c) {
                return _idx < ne*ps ? (_idx/ps)*ps*nc + ps*_c + _idx%ps : nc*_idx + _c; 
            };

            //  Brinkman term and diffusivity term
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;

                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __dkds = _mm256_loadu_pd(&_dkds[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                __m256d __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]), __imz = _mm256_loadu_pd(&_imz[idx]);
                __m256d __tem = _mm256_loadu_pd(&_tem[idx]), __item = _mm256_loadu_pd(&_item[idx]), __iqx = _mm256_loadu_pd(&_iqx[idx]), __iqy = _mm256_loadu_pd(&_iqy[idx]), __iqz = _mm256_loadu_pd(&_iqz[idx]);

                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_add_pd(_mm256_mul_pd(__uy, __imy), _mm256_mul_pd(__uz, __imz))))));

                int offsetf = nc*idx;
                __m256d __sumg = _mm256_setzero_pd(), __taug = _mm256_add_pd(_mm256_mul_pd(__3, _mm256_loadu_pd(&_diffusivity[idx])), _mm256_set1_pd(0.5));
                for (int c = 0; c < nc; ++c) {
                    int idxf = offsetf + ps*c;
                    __m256d __g = _mm256_loadu_pd(&_g[idxf]), __ig = _mm256_loadu_pd(&_ig[idxf]);
                    __sumg = _mm256_add_pd(__sumg, _mm256_mul_pd(__g, __ig));
                }
                __dfds = _mm256_sub_pd(__dfds, _mm256_div_pd(
                    _mm256_mul_pd(__3, _mm256_mul_pd(__dkds, _mm256_sub_pd(
                        __sumg, _mm256_mul_pd(__tem, _mm256_add_pd(__item, _mm256_mul_pd(
                            __3, _mm256_add_pd(_mm256_mul_pd(__ux, __iqx), _mm256_add_pd(_mm256_mul_pd(__uy, __iqy), _mm256_mul_pd(__uz, __iqz)))
                        )))
                    ))), _mm256_mul_pd(__taug, __taug)
                ));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                int offsetf = nc*idx;
                double sumg = 0.0;
                for (int c = 0; c < nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
            }

            SensitivityTemperatureAtHeatSourceAlongXFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along xmin
            SensitivityTemperatureAtHeatSourceAlongXFace(_q, _q.lx - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along xmax
            SensitivityTemperatureAtHeatSourceAlongYFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along ymin
            SensitivityTemperatureAtHeatSourceAlongYFace(_q, _q.ly - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along ymax
            SensitivityTemperatureAtHeatSourceAlongZFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along zmin
            SensitivityTemperatureAtHeatSourceAlongZFace(_q, _q.lz - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along zmax
        }
    }
}