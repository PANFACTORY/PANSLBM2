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
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T *, int);          //  Function of applying external force of NS with Brinkman model for 2D        
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T, T *, int);       //  Function of applying external force of NS with Brinkman model for 3D

        //  Function of updating macroscopic values of NS for 2D
        template<class P>
        void Macro(__m256d &__rho, __m256d &__ux, __m256d &__uy, const __m256d *__f) {
            __rho = __f[0];
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
        void Macro(__m256d &__rho, __m256d &__ux, __m256d &__uy, __m256d &__uz, const __m256d *__f) {
            __rho = __f[0];
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
        void Equilibrium(__m256d *__feq, const __m256d &__rho, const __m256d &__ux, const __m256d &__uy) {
            __m256d __1m15uu = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(1.5), _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy))));
            for (int c = 0; c < P::nc; ++c) {
                __m256d __cu = _mm256_add_pd(_mm256_mul_pd(P::__cx[c], __ux), _mm256_mul_pd(P::__cy[c], __uy));
                __feq[c] = _mm256_mul_pd(P::__ei[_c], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __cu), _mm256_mul_pd(_mm256_set1_pd(4.5), _mm256_mul_pd(__cu, __cu))))));
            }
        }

        //  Function of getting equilibrium of NS for 3D
        template<class P>
        void Equilibrium(__m256d *__feq, const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz) {
            __m256d __1m15uu = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(1.5), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy)), _mm256_mul_pd(__uz, __uz))));
            for (int c = 0; c < P::nc; ++c) {
                __m256d __cu = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(P::__cx[_c], __ux), _mm256_mul_pd(P::__cy[_c], __uy)), _mm256_mul_pd(P::__cz[_c], __uz));
                __feq[c] = _mm256_mul_pd(P::__ei[_c], _mm256_mul_pd(__rho, _mm256_add_pd(__1m15uu, _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __cu), _mm256_mul_pd(_mm256_set1_pd(4.5), _mm256_mul_pd(__cu, __cu))))));
            } 
        }

        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__alpha, __m256d *__f) {
            __m256d __coef = _mm256_div_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __alpha), __rho), _mm256_add_pd(__rho, __alpha));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_sub_pd(__f[c], _mm256_mul_pd(_mm256_mul_pd(__coef, P::__ei[c]), _mm256_add_pd(_mm256_mul_pd(P::__cx[c], __ux), _mm256_mul_pd(P::__cy[c], __uy))));
            }
        }

        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d &__alpha, __m256d *__f) {
            __m256d __coef = _mm256_div_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __alpha), __rho), _mm256_add_pd(__rho, __alpha));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_sub_pd(__f[c], __mm256_mul_pd(__mm256_mul_pd(__coef, P::__ei[c]), __mm256_add_pd(__mm256_add_pd(__mm256_mul_pd(P::__cx[c], __ux), __mm256_mul_pd(P::__cy[c], __uy)), __mm256_mul_pd(P::__cz[c], __uz))));
            }
        }

        //  Function of Update macro and Collide of NS for 2D
        template<template<class>class P>
        void MacroCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);

                //  Update macro
                __m256d __rho, __ux, __uy;
                Macro<P<double> >(__rho, __ux, __uy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy);
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[0]), _mm256_mul_pd(__omega, __feq[0])));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
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
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                Macro<P<double> >(__rho, __ux, __uy, __uz, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy, __uz);
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[0]), _mm256_mul_pd(__omega, __feq[0])));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy, uz;
                Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                }

                //  Collide
                Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of NS for 2D
        template<template<class>class P>
        void MacroBrinkmanCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double _viscosity, const double *_alpha, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
                
                //  Update macro
                __m256d __rho, __ux, __uy;
                Macro<P<double> >(__rho, __ux, __uy, __f);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __alpha, __f);
                Macro<P<double> >(__rho, __ux, __uy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy);
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[0]), _mm256_mul_pd(__omega, __feq[0])));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
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

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                Macro<double, P>(rho, ux, uy, _p.f0, _p.f, idx);

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
        void MacroBrinkmanCollide(P<double>& _p, double *_rho, double *_ux, double *_uy, double *_uz, double _viscosity, const double *_alpha, bool _issave = false) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                __f[0] = _mm256_load_pd(&_p.f0[idx]);
                P<double>::ShuffleToSoA(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);

                //  Update macro
                __m256d __rho, __ux, __uy, __uz;
                Macro<P<double> >(__rho, __ux, __uy, __uz, __f);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __alpha, __f);
                Macro<P<double> >(__rho, __ux, __uy, __uz, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_rho[idx], __rho);
                    _mm256_storeu_pd(&_ux[idx], __ux);
                    _mm256_storeu_pd(&_uy[idx], __uy);
                    _mm256_storeu_pd(&_uz[idx], __uz);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __rho, __ux, __uy, __uz);
                _mm256_store_pd(&_p.f0[idx], _mm256_add_pd(_mm256_mul_pd(__iomega, __f[0]), _mm256_mul_pd(__omega, __feq[0])));
                for (int c = 1; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                P<double>::ShuffleToAoS(&_p.f[P<double>::IndexF(idx, 1)], &__f[1]);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double rho, ux, uy, uz;
                Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                Macro<double, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                }

                //  Collide
                Equilibrium<double, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }
    }
}