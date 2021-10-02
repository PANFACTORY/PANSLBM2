//*****************************************************************************
//  Title       :   src/equation_avx/advection_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/13
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace AD {
        template<class T, template<class>class Q>void Macro(T &, T &, T &, T, T, const T *, const T *, T, int);         //  Function of updating macroscopic values of AD for 2D
        template<class T, template<class>class Q>void Macro(T &, T &, T &, T &, T, T, T, const T *, const T *, T, int); //  Function of updating macroscopic values of AD for 3D       
        template<class T, template<class>class Q>void Equilibrium(T *, T, T, T);                                        //  Function of getting equilibrium of AD for 2D
        template<class T, template<class>class Q>void Equilibrium(T *, T, T, T, T);                                     //  Function of getting equilibrium of AD for 3D        
        template<class T, template<class>class P>void ExternalForceNaturalConvection(T, T, T, T, T *, int);             //  Function of applying external force of AD with natural convection for 2D
        template<class T, template<class>class P>void ExternalForceNaturalConvection(T, T, T, T, T, T *, int);          //  Function of applying external force of AD with natural convection for 3D
        template<class T, template<class>class Q>void ExternalForceHeatExchange(T, T, T *, T *, int);                   //  Function of applying external force of AD with heat exchange for 2D/3D

        //  Function of updating macroscopic values of AD for 2D  
        template<class Q>
        void Macro(__m256d &__tem, __m256d &__qx, __m256d &__qy, __m256d __ux, __m256d __uy, const __m256d *__g, __m256d __omegag) {
            __tem = __g[0];
            __qx = _mm256_setzero_pd();
            __qy = _mm256_setzero_pd();
            for (int c = 1; c < Q::nc; ++c) {
                __tem = _mm256_add_pd(__tem, __g[c]);
                __qx = _mm256_add_pd(__qx, _mm256_mul_pd(Q::__cx[c], __g[c]));
                __qy = _mm256_add_pd(__qy, _mm256_mul_pd(Q::__cy[c], __g[c]));
            }
            __m256d __coef = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(0.5), __omegag));
            __qx = _mm256_mul_pd(__coef, _mm256_sub_pd(__qx, _mm256_mul_pd(__tem, __ux)));
            __qy = _mm256_mul_pd(__coef, _mm256_sub_pd(__qy, _mm256_mul_pd(__tem, __uy)));
        }

        //  Function of updating macroscopic values of AD for 3D  
        template<class Q>
        void Macro(__m256d &__tem, __m256d &__qx, __m256d &__qy, __m256d &__qz, __m256d __ux, __m256d __uy, __m256d __uz, const __m256d *__g, __m256d __omegag) {
            __tem = __g[0];
            __qx = _mm256_setzero_pd();
            __qy = _mm256_setzero_pd();
            __qz = _mm256_setzero_pd();
            for (int c = 1; c < Q::nc; ++c) {
                __tem = _mm256_add_pd(__tem, __g[c]);
                __qx = _mm256_add_pd(__qx, _mm256_mul_pd(Q::__cx[c], __g[c]));
                __qy = _mm256_add_pd(__qy, _mm256_mul_pd(Q::__cy[c], __g[c]));
                __qz = _mm256_add_pd(__qz, _mm256_mul_pd(Q::__cz[c], __g[c]));
            }
            __m256d __coef = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(0.5), __omegag));
            __qx = _mm256_mul_pd(__coef, _mm256_sub_pd(__qx, _mm256_mul_pd(__tem, __ux)));
            __qy = _mm256_mul_pd(__coef, _mm256_sub_pd(__qy, _mm256_mul_pd(__tem, __uy)));
            __qz = _mm256_mul_pd(__coef, _mm256_sub_pd(__qz, _mm256_mul_pd(__tem, __uz)));
        }

        //  Function of getting equilibrium of AD for 2D
        template<class Q>
        __m256d Equilibrium(__m256d __tem, __m256d __ux, __m256d __uy, int _c) {
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0);
            __m256d __ciu = _mm256_add_pd(_mm256_mul_pd(Q::__cx[_c], __ux), _mm256_mul_pd(Q::__cy[_c], __uy));
            return _mm256_mul_pd(Q::__ei[_c], _mm256_mul_pd(__tem, _mm256_add_pd(__1, _mm256_mul_pd(__3, __ciu))));
        }

        //  Function of getting equilibrium of AD for 3D
        template<class Q>
        __m256d Equilibrium(__m256d __tem, __m256d __ux, __m256d __uy, __m256d __uz, int _c) {
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0);
            __m256d __cu = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(Q::__cx[_c], __ux), _mm256_mul_pd(Q::__cy[_c], __uy)), _mm256_mul_pd(Q::__cz[_c], __uz));
            return _mm256_mul_pd(Q::__ei[_c], _mm256_mul_pd(__tem, _mm256_add_pd(__1, _mm256_mul_pd(__3, __cu))));
        }

        //  Function of applying external force of AD with natural convection for 2D
        template<class P>
        void ExternalForceNaturalConvection(__m256d __tem, __m256d __gx, __m256d __gy, __m256d __tem0, __m256d *__f) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_mul_pd(P::__ei[0], _mm256_sub_pd(__tem, __tem0)));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_add_pd(__f[c], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_mul_pd(P::__cx[c], __gx), _mm256_mul_pd(P::__cy[c], __gy))));
            }
        }

        //  Function of applying external force of AD with natural convection for 3D
        template<class P>
        void ExternalForceNaturalConvection(__m256d __tem, __m256d __gx, __m256d __gy, __m256d __gz, __m256d __tem0, __m256d *__f) {
            __m256d __coef = _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_mul_pd(P::__ei[0], _mm256_sub_pd(__tem, __tem0)));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_add_pd(__f[c], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(P::__cx[c], __gx), _mm256_mul_pd(P::__cy[c], __gy)), _mm256_mul_pd(P::__cz[c], __gz))));
            }
        }

        //  Function of applying external force of AD with heat exchange for 2D/3D
        template<class Q>
        void ExternalForceHeatExchange(__m256d __tem, __m256d __beta, __m256d *__g) {
            __m256d __1 = _mm256_set1_pd(1.0);
            __m256d __coef = _mm256_mul_pd(__beta, _mm256_div_pd(_mm256_sub_pd(__1, __tem), _mm256_add_pd(__1, __beta)));
            for (int c = 0; c < Q::nc; ++c) {
                __g[c] = _mm256_add_pd(__g[c], _mm256_mul_pd(Q::__ei[c], __coef));
            }
        }

        //  Function of Update macro and Collide of force convection for 2D
        template<template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double _diffusivity, 
            bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __tem, __qx, __qy;
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy;
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                double tem, qx, qy;
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of force convection for 3D
        template<template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double *_qz, double _diffusivity, 
            bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __tem, __qx, __qy, __qz;
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                    _mm256_storeu_pd(&_qz[idx], __qz);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy, uz;
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                double tem, qx, qy, qz;
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of natural convection for 2D
        template<template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double _diffusivity, 
            double _gx, double _gy, double _tem0, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __tem0 = _mm256_set1_pd(_tem0);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __tem, __qx, __qy;
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  External force with natural convection
                ExternalForceNaturalConvection<P<double> >(__tem, __gx, __gy, __tem0, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy;
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                double tem, qx, qy;
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<double, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of natural convection for 3D
        template<template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double *_qz, double _diffusivity, 
            double _gx, double _gy, double _gz, double _tem0, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            double omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<double>::nc];
            __m256d __omegag = _mm256_set1_pd(omegag), __iomegag = _mm256_set1_pd(iomegag);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __gz = _mm256_set1_pd(_gz), __tem0 = _mm256_set1_pd(_tem0);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __tem, __qx, __qy, __qz;
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  External force with natural convection
                ExternalForceNaturalConvection<P<double> >(__tem, __gx, __gy, __gz, __tem0, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                    _mm256_storeu_pd(&_qz[idx], __qz);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy, uz;
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                double tem, qx, qy, qz;
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<double, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c); 
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of Brinkman and heat exchange for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<double>& _p, double *_rho, double *_ux, double *_uy, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, const double *_beta, double _diffusivity, bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __tem, __qx, __qy;
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  External force with Brinkman and heat exchange
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __beta = _mm256_loadu_pd(&_beta[idx]);
                ExternalForceHeatExchange<Q<double> >(__tem, __beta, __g);
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy;
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                double tem, qx, qy;
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<double, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and heat exchange for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double *_qz, const double *_beta, double _diffusivity, bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __tem, __qx, __qy, __qz;
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  External force with Brinkman and heat exchange
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __beta = _mm256_loadu_pd(&_beta[idx]);
                ExternalForceHeatExchange<Q<double> >(__tem, __beta, __g);
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                    _mm256_storeu_pd(&_qz[idx], __qz);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy, uz;
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                double tem, qx, qy, qz;
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<double, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {            
                    int idxf = Q<double>::IndexF(idx, c);      
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and force convection for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, const double *_diffusivity, 
            bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __tem, __qx, __qy;
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  External force with Brinkman
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double rho, ux, uy;
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                double tem, qx, qy;
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and force convection for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double *_qz, const double *_diffusivity, 
            bool _issave = false
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
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __tem, __qx, __qy, __qz;
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  External force with Brinkman and heat exchange
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                    _mm256_storeu_pd(&_qz[idx], __qz);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double rho, ux, uy, uz;
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                double tem, qx, qy, qz;
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection for 2D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, const double *_diffusivity, 
            double _gx, double _gy, double _tem0, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __tem0 = _mm256_set1_pd(_tem0);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                __m256d __tem, __qx, __qy;
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  External force with Brinkman and heat exchange
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceNaturalConvection<P<double> >(__tem, __gx, __gy, __tem0, __f);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __f);
                Macro<Q<double> >(__tem, __qx, __qy, __ux, __uy, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double rho, ux, uy;
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                double tem, qx, qy;
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<double, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<double, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection for 3D
        template<template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, const double *_alpha, double _viscosity,
            Q<double>& _q, double *_tem, double *_qx, double *_qy, double *_qz, const double *_diffusivity, 
            double _gx, double _gy, double _gz, double _tem0, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<double>::nc], geq[Q<double>::nc];
            __m256d __omegaf = _mm256_set1_pd(omegaf), __iomegaf = _mm256_set1_pd(iomegaf);
            __m256d __gx = _mm256_set1_pd(_gx), __gy = _mm256_set1_pd(_gy), __gz = _mm256_set1_pd(_gz), __tem0 = _mm256_set1_pd(_tem0);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                __m256d __diffusivity = _mm256_loadu_pd(&_diffusivity[idx]);                 
                __m256d __omegag = _mm256_div_pd(_mm256_set1_pd(1.0), _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __diffusivity), _mm256_set1_pd(0.5)));
                __m256d __iomegag = _mm256_sub_pd(_mm256_set1_pd(1.0), __omegag);

                //  Pack f0, f, g0 and g
                __m256d __f[P<double>::nc], __g[Q<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                __g[0] = _mm256_load_pd(&_q.f0[idx]);
                Q<double>::ShuffleToSoA(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                __m256d __tem, __qx, __qy, __qz;
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  External force with Brinkman and heat exchange
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceNaturalConvection<P<double> >(__tem, __gx, __gy, __gz, __tem0, __f);
                NS::ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __alpha, __f);
                NS::Macro<P<double> >(__rho, __ux, __uy, __uz, __f);
                Macro<Q<double> >(__tem, __qx, __qy, __qz, __ux, __uy, __uz, __g, __omegag);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                    _mm256_storeu_pd(&_tem[idx], __tem);
                    _mm256_storeu_pd(&_qx[idx], __qx);
                    _mm256_storeu_pd(&_qy[idx], __qy);
                    _mm256_storeu_pd(&_qz[idx], __qz);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[0]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomegaf, __f[c]), _mm256_mul_pd(__omegaf, NS::Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                _mm256_store_pd(&_q.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[0]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, 0))));
                for (int c = 1; c < Q<double>::nc; ++c) {
                    __g[c] = _mm256_add_pd(_mm256_mul_pd(__iomegag, __g[c]), _mm256_mul_pd(__omegag, Equilibrium<Q<double> >(__tem, __ux, __uy, __uz, c)));
                }
                Q<double>::ShuffleToAoS(&_q.f[Q<double>::IndexF(idx, 1)], &__g[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                double omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                double rho, ux, uy, uz;
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                double tem, qx, qy, qz;
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<double, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<double, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<double, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<double, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<double>::nc; ++c) {
                    int idxf = Q<double>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    }
}