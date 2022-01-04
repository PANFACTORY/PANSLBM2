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
    namespace ANS {
        template<class T, template<class>class P>void Macro(T &, T &, T &, T &, T &, T, T, T, const T *, const T *, int);               //  Function of updating macroscopic values of ANS for 2D
        template<class T, template<class>class P>void Macro(T &, T &, T &, T &, T &, T &, T &, T, T, T, T, const T *, const T *, int);  //  Function of updating macroscopic values of ANS for 3D
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T, T, T);                                                  //  Function of getting equilibrium of ANS for 2D
        template<class T, template<class>class P>void Equilibrium(T *, T, T, T, T, T, T, T);                                            //  Function of getting equilibrium of ANS for 3D
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T, T *, T *, T, int);                           //  Function of applying external force with Brinkman model of ANS for 2D
        template<class T, template<class>class P>void ExternalForceBrinkman(T, T, T, T, T, T, T, T *, T *, T, int);                     //  Function of applying external force with Brinkman model of ANS for 3D

        //  Function of updating macroscopic values of ANS for 2D
        template<class P>
        void Macro(__m256d &__ip, __m256d &__iux, __m256d &__iuy, __m256d &__imx, __m256d &__imy, const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d *__f) {           
            __ip = _mm256_setzero_pd();
            __iux = _mm256_setzero_pd();
            __iuy = _mm256_setzero_pd();
            __imx = _mm256_setzero_pd();
            __imy = _mm256_setzero_pd();

            __m256d __1uu = _mm256_sub_pd(_mm256_set1_pd(1.0), _mm256_mul_pd(_mm256_set1_pd(1.5), _mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy))));
            
            for (int c = 0; c < P::nc; ++c) {
                __m256d __fei = _mm256_mul_pd(__f[c], P::__ei[c]);
                __m256d __cu = _mm256_add_pd(_mm256_mul_pd(P::__cx[c], __ux), _mm256_mul_pd(P::__cy[c], __uy));
                __ip = _mm256_add_pd(__ip, _mm256_mul_pd(__fei, _mm256_add_pd(__1uu, _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), __cu), _mm256_mul_pd(_mm256_set1_pd(4.5), _mm256_mul_pd(__cu, __cu))))));               
                __iux = _mm256_add_pd(__iux, _mm256_mul_pd(__fei, _mm256_add_pd(P::__cx[c], _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_mul_pd(__cu, P::__cx[c])), __ux))));
                __iuy = _mm256_add_pd(__iuy, _mm256_mul_pd(__fei, _mm256_add_pd(P::__cy[c], _mm256_sub_pd(_mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_mul_pd(__cu, P::__cy[c])), __uy))));
                __imx = _mm256_add_pd(__imx, _mm256_mul_pd(__fei, P::__cx[c]));
                __imy = _mm256_add_pd(__imy, _mm256_mul_pd(__fei, P::__cy[c]));
            }
        }

        //  Function of updating macroscopic values of ANS for 3D
        template<class P>
        void Macro(__m256d &__ip, __m256d &__iux, __m256d &__iuy, __m256d &__iuz, __m256d &__imx, __m256d &__imy, __m256d &__imz, const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d *__f) {
            __ip = _mm256_setzero_pd();
            __iux = _mm256_setzero_pd();
            __iuy = _mm256_setzero_pd();
            __iuz = _mm256_setzero_pd();
            __imx = _mm256_setzero_pd();
            __imy = _mm256_setzero_pd();
            __imz = _mm256_setzero_pd();

            __m256d __uu = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__ux, __ux), _mm256_mul_pd(__uy, __uy)), _mm256_mul_pd(__uz, __uz));
            __m256d __1 = _mm256_set1_pd(1.0), __3 = _mm256_set1_pd(3.0), __45 = _mm256_set1_pd(4.5), __15 = _mm256_set1_pd(1.5);

            for (int c = 0; c < P::nc; ++c) {
                __m256d __fei = _mm256_mul_pd(__f[c], P::__ei[c]);
                __m256d __cu = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(P::__cx[c], __ux), _mm256_mul_pd(P::__cy[c], __uy)), _mm256_mul_pd(P::__cz[c], __uz));
                __ip = _mm256_add_pd(__ip, _mm256_mul_pd(__fei, _mm256_add_pd(__1, _mm256_add_pd(_mm256_mul_pd(__3, __cu), _mm256_sub_pd(_mm256_mul_pd(__45, _mm256_mul_pd(__cu, __cu)), _mm256_mul_pd(__15, __uu))))));
                __iux = _mm256_add_pd(__iux, _mm256_mul_pd(__fei, _mm256_add_pd(P::__cx[c], _mm256_sub_pd(_mm256_mul_pd(__3, _mm256_mul_pd(__cu, P::__cx[c])), __ux))));
                __iuy = _mm256_add_pd(__iuy, _mm256_mul_pd(__fei, _mm256_add_pd(P::__cy[c], _mm256_sub_pd(_mm256_mul_pd(__3, _mm256_mul_pd(__cu, P::__cy[c])), __uy))));
                __iuz = _mm256_add_pd(__iuz, _mm256_mul_pd(__fei, _mm256_add_pd(P::__cz[c], _mm256_sub_pd(_mm256_mul_pd(__3, _mm256_mul_pd(__cu, P::__cz[c])), __uz))));
                __imx = _mm256_add_pd(__imx, _mm256_mul_pd(__fei, P::__cx[c]));
                __imy = _mm256_add_pd(__imy, _mm256_mul_pd(__fei, P::__cy[c]));
                __imz = _mm256_add_pd(__imz, _mm256_mul_pd(__fei, P::__cz[c]));
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class P>
        void Equilibrium(__m256d *__feq, const __m256d &__ux, const __m256d &__uy, const __m256d &__ip, const __m256d &__iux, const __m256d &__iuy) {
            for (int c = 0; c < P::nc; ++c) {
                __feq[c] = _mm256_add_pd(__ip, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_mul_pd(__iux, _mm256_sub_pd(P::__cx[c], __ux)), _mm256_mul_pd(__iuy, _mm256_sub_pd(P::__cy[c], __uy))))); 
            }
        }

        //  Function of getting equilibrium of ANS for 3D
        template<class P>
        void Equilibrium(__m256d *__feq, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d &__ip, const __m256d &__iux, const __m256d &__iuy, const __m256d &__iuz) {
            for (int c = 0; c < P::nc; ++c) {
                __feq[c] = _mm256_add_pd(__ip, _mm256_mul_pd(_mm256_set1_pd(3.0), _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__iux, _mm256_sub_pd(P::__cx[c], __ux)), _mm256_mul_pd(__iuy, _mm256_sub_pd(P::__cy[c], __uy))), _mm256_mul_pd(__iuz, _mm256_sub_pd(P::__cz[c], __uz)))));
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 2D
        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__imx, const __m256d &__imy, __m256d *__f, const __m256d &__alpha) {
            __m256d __3 = _mm256_set1_pd(3.0);
            __m256d __coef = _mm256_mul_pd(__3, _mm256_div_pd(__alpha, _mm256_add_pd(__rho, __alpha)));
            __f[0] = _mm256_add_pd(__f[0], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy))));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_sub_pd(__f[c], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(P::__cx[c], __ux), __imx), _mm256_mul_pd(_mm256_sub_pd(P::__cy[c], __uy), __imy))));
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 3D
        template<class P>
        void ExternalForceBrinkman(const __m256d &__rho, const __m256d &__ux, const __m256d &__uy, const __m256d &__uz, const __m256d &__imx, const __m256d &__imy, const __m256d &__imz, __m256d *__f, const __m256d &__alpha) {
            __m256d __3 = _mm256_set1_pd(3.0);
            __m256d __coef = _mm256_mul_pd(__3, _mm256_div_pd(__alpha, _mm256_add_pd(__rho, __alpha)));
            __f[0] = _mm256_add_pd(__f[0], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy)), _mm256_mul_pd(__uz, __imz))));
            for (int c = 1; c < P::nc; ++c) {
                __f[c] = _mm256_sub_pd(__f[c], _mm256_mul_pd(__coef, _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_sub_pd(P::__cx[c], __ux), __imx), _mm256_mul_pd(_mm256_sub_pd(P::__cy[c], __uy), __imy)), _mm256_mul_pd(_mm256_sub_pd(P::__cz[c], __uz), __imz))));
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of ANS for 2D
        template<template<class>class P>
        void MacroBrinkmanCollide(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, 
            double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, 
            double _viscosity, const double *_alpha, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
             
                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]);
                Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __imx, __imy, __f, __alpha);
                Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, __ux, __uy, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __ux, __uy, __ip, __iux, __iuy);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, imx, imy;
                Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _p.f0, _p.f, _alpha[idx], idx);
                Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                }

                //  Collide
                Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of ANS for 3D
        template<template<class>class P>
        void MacroBrinkmanCollide(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, 
            double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, 
            double _viscosity, const double *_alpha, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
                
                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);

                //  External force with Brinkman model
                __m256d __alpha = _mm256_loadu_pd(&_alpha[idx]);
                ExternalForceBrinkman<P<double> >(__rho, __ux, __uy, __uz, __imx, __imy, __imz, __f, __alpha);
                Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, __ux, __uy, __uz, __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, __iux, __iuy, __iuz);
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<double, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _p.f0, _p.f, _alpha[idx], idx);
                Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                }

                //  Collide and stream
                Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        } 
    
        //  Function of Update macro and Collide of ANS with Level-Set-Method for 2D
        template<template<class>class P>
        void MacroCollideLSM(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, 
            double *_ip, double *_iux, double *_iuy, double *_imx, double *_imy, 
            double _viscosity, const double *_chi, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
                __m256d __chi = _mm256_loadu_pd(&_chi[idx]);
             
                //  Update macro
                __m256d __ip, __iux, __iuy, __imx, __imy;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]);
                Macro<P<double> >(__ip, __iux, __iuy, __imx, __imy, __rho, _mm256_mul_pd(__ux, __chi), _mm256_mul_pd(__uy, __chi), __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __ux, __uy, __ip, _mm256_mul_pd(__iux, __chi), _mm256_mul_pd(__iuy, __chi));
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, imx, imy;
                Macro<double, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx]*_chi[idx], _uy[idx]*_chi[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                }

                //  Collide
                Equilibrium<double, P>(feq, _ux[idx], _uy[idx], ip, iux*_chi[idx], iuy*_chi[idx]);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro and Collide of ANS with Level-Set-Method for 3D
        template<template<class>class P>
        void MacroCollideLSM(
            P<double>& _p, const double *_rho, const double *_ux, const double *_uy, const double *_uz, 
            double *_ip, double *_iux, double *_iuy, double *_iuz, double *_imx, double *_imy, double *_imz, 
            double _viscosity, const double *_chi, bool _issave = false
        ) {
            const int ne = _p.nxyz/P<double>::packsize;
            double omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<double>::nc];
            __m256d __omega = _mm256_set1_pd(omega), __iomega = _mm256_set1_pd(iomega), __feq[P<double>::nc];
            #pragma omp parallel for private(__feq)
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*P<double>::packsize;

                //  Pack f0 and f
                __m256d __f[P<double>::nc];
                _p.LoadF(idx, __f);
                __m256d __chi = _mm256_loadu_pd(&_chi[idx]);
                
                //  Update macro
                __m256d __ip, __iux, __iuy, __iuz, __imx, __imy, __imz;
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                Macro<P<double> >(__ip, __iux, __iuy, __iuz, __imx, __imy, __imz, __rho, _mm256_mul_pd(__ux, __chi), _mm256_mul_pd(__uy, __chi), _mm256_mul_pd(__uz, __chi), __f);

                //  Save macro if need
                if (_issave) {
                    _mm256_storeu_pd(&_ip[idx], __ip);
                    _mm256_storeu_pd(&_iux[idx], __iux);
                    _mm256_storeu_pd(&_iuy[idx], __iuy);
                    _mm256_storeu_pd(&_iuz[idx], __iuz);
                    _mm256_storeu_pd(&_imx[idx], __imx);
                    _mm256_storeu_pd(&_imy[idx], __imy);
                    _mm256_storeu_pd(&_imz[idx], __imz);
                }

                //  Collide
                Equilibrium<P<double> >(__feq, __ux, __uy, __uz, __ip, _mm256_mul_pd(__iux, __chi), _mm256_mul_pd(__iuy, __chi), _mm256_mul_pd(__iuz, __chi));
                for (int c = 0; c < P<double>::nc; ++c) {
                    __f[c] = _mm256_add_pd(_mm256_mul_pd(__iomega, __f[c]), _mm256_mul_pd(__omega, __feq[c]));
                }
                _p.StoreF(idx, __f);
            }
            for (int idx = ne*P<double>::packsize; idx < _p.nxyz; ++idx) {
                //  Update macro
                double ip, iux, iuy, iuz, imx, imy, imz;
                Macro<double, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx]*_chi[idx], _uy[idx]*_chi[idx], _uz[idx]*_chi[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                }

                //  Collide and stream
                Equilibrium<double, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux*_chi[idx], iuy*_chi[idx], iuz*_chi[idx]);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<double>::nc; ++c) {
                    int idxf = P<double>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        } 

        //  Function of getting sensitivity of Brinkman model for 2D
        template<template<class>class P>
        void SensitivityBrinkman(P<double>& _p, double *_dfds, const double *_ux, const double *_uy, const double *_imx, const double *_imy, const double *_dads) {
            const int ps = P<double>::packsize, ne = _p.nxyz/ps, nc = P<double>::nc;
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]);
                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_mul_pd(__uy, __imy)))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _p.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
            }
        }

        //  Function of getting sensitivity of Brinkman model for 3D
        template<template<class>class P>
        void SensitivityBrinkman(P<double>& _p, double *_dfds, const double *_ux, const double *_uy, const double *_uz, const double *_imx, const double *_imy, const double *_imz, const double *_dads) {
            const int ps = P<double>::packsize, ne = _p.nxyz/ps, nc = P<double>::nc;
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __dads = _mm256_loadu_pd(&_dads[idx]), __3 = _mm256_set1_pd(3.0);
                __m256d __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                __m256d __imx = _mm256_loadu_pd(&_imx[idx]), __imy = _mm256_loadu_pd(&_imy[idx]), __imz = _mm256_loadu_pd(&_imz[idx]);
                __dfds = _mm256_add_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__dads, _mm256_add_pd(_mm256_mul_pd(__ux, __imx), _mm256_add_pd(_mm256_mul_pd(__uy, __imy), _mm256_mul_pd(__uz, __imz))))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _p.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
            }
        }
    
        //  Function of getting sensitivity with Level-Set-Method for 2D
        template<template<class>class P>
        void SensitivityLSM(P<double>& _p, double *_dfds, const double *_rho, const double *_ux, const double *_uy, const double *_iux, const double *_iuy, double _viscosity) {
            const int ps = P<double>::packsize, ne = _p.nxyz/ps, nc = P<double>::nc;
            const double omega = 1.0/(3.0*_viscosity + 0.5);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __omega = _mm256_set1_pd(omega), __3 = _mm256_set1_pd(3.0);
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __iux = _mm256_loadu_pd(&_iux[idx]), __iuy = _mm256_loadu_pd(&_iuy[idx]);
                __dfds = _mm256_sub_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__omega, _mm256_mul_pd(__rho, _mm256_add_pd(_mm256_mul_pd(__ux, __iux), _mm256_mul_pd(__uy, __iuy))))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _p.nxyz; ++idx) {
                _dfds[idx] -= 3.0*omega*_rho[idx]*(_ux[idx]*_iux[idx] + _uy[idx]*_iuy[idx]);
            }
        }

        //  Function of getting sensitivity with Level-Set-Method for 3D
        template<template<class>class P>
        void SensitivityLSM(P<double>& _p, double *_dfds, const double *_rho, const double *_ux, const double *_uy, const double *_uz, const double *_iux, const double *_iuy, const double *_iuz, double _viscosity) {
            const int ps = P<double>::packsize, ne = _p.nxyz/ps, nc = P<double>::nc;
            const double omega = 1.0/(3.0*_viscosity + 0.5);
            #pragma omp parallel for
            for (int pidx = 0; pidx < ne; ++pidx) {
                int idx = pidx*ps;
                __m256d __dfds = _mm256_loadu_pd(&_dfds[idx]), __omega = _mm256_set1_pd(omega), __3 = _mm256_set1_pd(3.0);
                __m256d __rho = _mm256_loadu_pd(&_rho[idx]), __ux = _mm256_loadu_pd(&_ux[idx]), __uy = _mm256_loadu_pd(&_uy[idx]), __uz = _mm256_loadu_pd(&_uz[idx]);
                __m256d __iux = _mm256_loadu_pd(&_iux[idx]), __iuy = _mm256_loadu_pd(&_iuy[idx]), __iuz = _mm256_loadu_pd(&_iuz[idx]);
                __dfds = _mm256_sub_pd(__dfds, _mm256_mul_pd(__3, _mm256_mul_pd(__omega, _mm256_mul_pd(__rho, _mm256_add_pd(_mm256_mul_pd(__ux, __iux), _mm256_add_pd(_mm256_mul_pd(__uy, __iuy), _mm256_mul_pd(__uz, __iuz)))))));
                _mm256_storeu_pd(&_dfds[idx], __dfds);
            }
            for (int idx = ne*ps; idx < _p.nxyz; ++idx) {
                _dfds[idx] -= 3.0*omega*_rho[idx]*(_ux[idx]*_iux[idx] + _uy[idx]*_iuy[idx] + _uz[idx]*_iuz[idx]);
            }
        }
    }  
}