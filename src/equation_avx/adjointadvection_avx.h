//*****************************************************************************
//  Title       :   src/equation_avx/adjointadvection_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <immintrin.h>
#include "adjointnavierstokes_avx.h"

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
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
        __m256d Equilibrium(__m256d __item, __m256d __iqx, __m256d __iqy, __m256d __ux, __m256d __uy) {
            return _mm256_add_pd(__item, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_mul_pd(__iqx, __ux), _mm256_mul_pd(__iqy, __uy))));
        }

        //  Function of getting equilibrium of AAD for 3D
        template<class Q>
        __m256d Equilibrium(__m256d __item, __m256d __iqx, __m256d __iqy, __m256d __iqz, __m256d __ux, __m256d __uy, __m256d __uz) {
            return _mm256_add_pd(__item, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__iqx, __ux), _mm256_mul_pd(__iqy, __uy)), _mm256_mul_pd(__iqz, __uz))));
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class P>
        void ExternalForceBrinkman(__m256d __rho, __m256d __ux, __m256d __uy, __m256d __imx, __m256d __imy, __m256d __tem, __m256d __iqx, __m256d __iqy, __m256d __omegag, __m256d *__f, __m256d __alpha) {
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
        void ExternalForceBrinkman(__m256d __rho, __m256d __ux, __m256d __uy, __m256d __uz, __m256d __imx, __m256d __imy, __m256d __imz, __m256d __tem, __m256d __iqx, __m256d __iqy, __m256d __iqz, __m256d __omegag, __m256d *__f, __m256d __alpha) {
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
        void ExternalForceHeatExchange(__m256d __item, __m256d *__g, __m256d __beta) {
            __m256d __coef = _mm256_mul_pd(__beta, _mm256_div_pd(_mm256_add_pd(_mm256_set1_pd(1.0), __item), _mm256_add_pd(_mm256_set1_pd(1.0), __beta)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_sub_pd(__g[c], __coef);
            }
        }

        //  Function of applying external force with natural convection of AAD for 2D
        template<class Q>
        void ExternalForceNaturalConvection(__m256d __imx, __m256d __imy, __m256d __gx, __m256d __gy, __m256d *__g) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_mul_pd(__imx, __gx), _mm256_mul_pd(__imy, __gy)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_add_pd(__g[c], __coef);
            }
        }

        //  Function of applying external force with natural convection of AAD for 3D
        template<class Q>
        void ExternalForceNaturalConvection(__m256d __imx, __m256d __imy, __m256d __imz, __m256d __gx, __m256d __gy, __m256d __gz, __m256d *__g) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__imx, __gx), _mm256_mul_pd(__imy, __gy)), _mm256_mul_pd(__imz, __gz)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_add_pd(__g[c], __coef);
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
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]);
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
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
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
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, double *_ip, double *_iux, double *_iuy, double _iuz, double *_imx, double *_imy, double *_imz, const double *_alpha, double _viscosity,
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_beta, double _diffusivity, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
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
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
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
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_diffusivity, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]);
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
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
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
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
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
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
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
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, const double *_diffusivity, double _gx, double _gy, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]);
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
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __ip, __iux, __iuy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double ip, iux, iuy, imx, imy;
                ANS::Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                double item, iqx, iqy;
                Macro<double, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
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
            Q<double>& _q, const double *_tem, double *_item, double *_iqx, double *_iqy, double *_iqz, const double *_diffusivity, double _gx, double _gy, double _gz, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __m256d __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __m256d __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
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
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, ANS::Equilibrium<P<double> >(__ux, __uy, __uz, __ip, __iux, __iuy, __iuz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__item, __iqx, __iqy, __iqz, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
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
    }
}