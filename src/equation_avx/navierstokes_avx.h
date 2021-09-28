//*****************************************************************************
//  Title       :   src/equation_avx/navierstokes_avx.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/01/24
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <immintrin.h>

//  compile option for g++(MinGW) : -mavx

namespace PANSLBM2 {
    namespace NS {
        template<class T, template<class>class P>void Macro(T &, T &, T &, const T *, const T *, int);      //  Function of updating macroscopic values of NS for 2D
        template<class T, template<class>class P>void Macro(T &, T &, T &, T &, const T *, const T *, int); //  Function of updating macroscopic values of NS for 3D        
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T);                            //  Function of getting equilibrium of NS for 2D  
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T, T);                         //  Function of getting equilibrium of NS for 3D

        //  Function of updating macroscopic values of NS for 2D
        template<class P>
        void Macro(__m256d &__rho, __m256d &__ux, __m256d &__uy, __m256d __f0, const __m256d *__f) {
            __rho = __f0;
            __ux = _mm256_setzero_pd();
            __uy = _mm256_setzero_pd();
            for (int c = 1; c < P::nc; ++c) {
                __rho = _mm256_add_pd(__rho, __f[c]);
                __ux = _mm256_add_pd(__ux, _mm256_mul_pd(__f[c], P::__cx[c]));
                __uy = _mm256_add_pd(__uy, _mm256_mul_pd(__f[c], P::__cy[c]));
            }
            __m256d __invrho = _mm256_div_pd(_mm256_set1_pd(1.0), __rho);
            __ux = _mm256_mul_pd(__ux, __invrho);
            __uy = _mm256_mul_pd(__uy, __invrho);
        }

        //  Function of updating macroscopic values of NS for 3D
        template<class P>
        void Macro(__m256d &__rho, __m256d &__ux, __m256d &__uy, __m256d &__uz, __m256d __f0, const __m256d *__f) {
            __rho = __f0;
            __ux = _mm256_setzero_pd();
            __uy = _mm256_setzero_pd();
            __uz = _mm256_setzero_pd();
            for (int c = 1; c < P::nc; ++c) {
                __rho = _mm256_add_pd(__rho, __f[c]);
                __ux = _mm256_add_pd(__ux, _mm256_mul_pd(__f[c], P::__cx[c]));
                __uy = _mm256_add_pd(__uy, _mm256_mul_pd(__f[c], P::__cy[c]));
                __uz = _mm256_add_pd(__uz, _mm256_mul_pd(__f[c], P::__cz[c]));
            }
            __m256d __invrho = _mm256_div_pd(_mm256_set1_pd(1.0), __rho);
            __ux = _mm256_mul_pd(__ux, __invrho);
            __uy = _mm256_mul_pd(__uy, __invrho);
            __uz = _mm256_mul_pd(__uz, __invrho);
        }

        //  Function of getting equilibrium of NS for 2D
        template<class P>
        __m256d Equilibrium(__m256d __rho, __m256d __ux, __m256d __uy, int _c) {
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0), __45 = _mm256_set1_pd(4.5), __15 = _mm256_set1_pd(1.5);
            __m256d __1m15uu = _mm256_sub_pd(__1, _mm256_mul_pd(__15, _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy))));
            __m256d __cu = _mm256_add_pd(_mm256_mul_pd(P::__cx[_c], __ux), _mm256_mul_pd(P::__cy[_c], __uy));
            __m256d __3cup45cucu = _mm256_add_pd(_mm256_mul_pd(__3, __cu), _mm256_mul_pd(__45, _mm256_mul_pd(__cu, __cu)));
            return _mm256_mul_pd(P::__ei[_c], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, __3cup45cucu)));
        }

        //  Function of getting equilibrium of NS for 3D
        template<class P>
        __m256d Equilibrium(__m256d __rho, __m256d __ux, __m256d __uy, __m256d __uz, int _c) {
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0), __45 = _mm256_set1_pd(4.5), __15 = _mm256_set1_pd(1.5);
            __m256d __1m15uu = _mm256_sub_pd(__1, _mm256_mul_pd(__15, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy)), _mm256_mul_pd(__uz, __uz))));
            __m256d __cu = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(P::__cx[_c], __ux), _mm256_mul_pd(P::__cy[_c], __uy)), _mm256_mul_pd(P::__cz[_c], __uz));
            __m256d __3cup45cucu = _mm256_add_pd(_mm256_mul_pd(__3, __cu), _mm256_mul_pd(__45, _mm256_mul_pd(__cu, __cu)));
            return _mm256_mul_pd(P::__ei[_c], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, __3cup45cucu)));
        }

        //  Function of Update macro and Collide of NS for 2D
        template<template<class>class P>
        void MacroCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f0 = _mm256_load_pd(&_p.f0[idx]), __f[P<double>::nc];
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                Macro<P<double> >(__rho, __ux, __uy, __f0, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f0), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, c)));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy;
                Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                }

                //  Collide
                Equilibrium<double, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro and Collide of NS for 3D
        template<template<class>class P>
        void MacroCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, double _viscosity, bool _issave = false) {
            double omega = 1.0/(3.0*_viscosity + 0.5);
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(1.0 - omega);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; idx += P<double>::packsize) {
                //  Pack f0 and f
                __m256d __f0 = _mm256_load_pd(&_p.f0[idx]), __f[P<double>::nc] = { 0 };
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_load_pd(&_p.f[P<double>::IndexF(idx, c)]);
                }

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                Macro<P<double> >(__rho, __ux, __uy, __uz, __f0, __f);

                //  Save macro if need
                if (_issave) {
                    if (idx + P<double>::packsize <= _p.nxyz) {
                        _mm256_storeu_pd(&_rho[idx], __rho);
                        _mm256_storeu_pd(&_ux[idx], __ux);
                        _mm256_storeu_pd(&_uy[idx], __uy);
                        _mm256_storeu_pd(&_uz[idx], __uz);
                    } else {
                        double rho[P<double>::packsize], ux[P<double>::packsize], uy[P<double>::packsize], uz[P<double>::packsize];
                        _mm256_storeu_pd(rho, __rho);
                        _mm256_storeu_pd(ux, __ux);
                        _mm256_storeu_pd(uy, __uy);
                        _mm256_storeu_pd(uz, __uz);
                        for (int didx = 0; didx < _p.nxyz%P<double>::packsize; ++didx) {
                            _rho[idx + didx] = rho[didx];
                            _ux[idx + didx] = ux[didx];
                            _uy[idx + didx] = uy[didx];
                            _uz[idx + didx] = uz[didx];
                        }
                    }
                }

                //  Collide
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f0), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, __uz, 0))));
                for (int c = 1; c < P<double>::nc; ++c) {
                    _mm256_store_pd(&_p.f[P<double>::IndexF(idx, c)], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, Equilibrium<P<double> >(__rho, __ux, __uy, __uz, c))));
                }
            }
        }
    }
}