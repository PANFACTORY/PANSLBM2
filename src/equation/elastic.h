//*****************************************************************************
//  Title       :   src/equation/elastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/navierstokes_avx.h"
#endif

namespace PANSLBM2 {
    namespace EL {
        //  Function of updating macroscopic values of EL for 2D
        template<class T, template<class>class P>
        void Macro(T _rho, T &_ux, T &_uy, T &_sxx, T &_sxy, T &_syx, T &_syy, const T *_f, int _idx) {
            _ux = T();
            _uy = T();
            _sxx = T();
            _sxy = T();
            _syx = T();
            _syy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T f = _f[P<T>::IndexF(_idx, c)];
                _ux += P<T>::cx[c]*f;
                _uy += P<T>::cy[c]*f;
                _sxx -= P<T>::cx[c]*P<T>::cx[c]*f;
                _sxy -= P<T>::cx[c]*P<T>::cy[c]*f;
                _syx -= P<T>::cy[c]*P<T>::cx[c]*f;
                _syy -= P<T>::cy[c]*P<T>::cy[c]*f;
            }
            T invrho = 1.0/_rho;
            _ux *= invrho;
            _uy *= invrho;
        }

        //  Function of updating macroscopic values of EL with topology optimization for 2D
        template<class T, template<class>class P>
        void Macro(T _rho, T &_ux, T &_uy, T &_sxx, T &_sxy, T &_syx, T &_syy, const T *_f, T _gamma, int _idx) {
            Macro<T, P>(_rho, _ux, _uy, _sxx, _sxy, _syx, _syy, _f, _idx);
            _sxx *= _gamma;
            _sxy *= _gamma;
            _syx *= _gamma;
            _syy *= _gamma;
        }
    
        //  Function of getting equilibrium of EL for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _rho, T _ux, T _uy, T _sxx, T _sxy, T _syx, T _syy) {
            for (int c = 0; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy;
                T cisci = P<T>::cx[c]*_sxx*P<T>::cx[c] + P<T>::cx[c]*_sxy*P<T>::cy[c] + P<T>::cy[c]*_syx*P<T>::cx[c] + P<T>::cy[c]*_syy*P<T>::cy[c];
                T trss = _sxx + _syy;
                _feq[c] = P<T>::ei[c]*(3.0*_rho*ciu - 4.5*cisci + 1.5*trss);
            }
        }

        //  Function of setting boundary condition of EL set Stress for 2D along x edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetStressAlongXEdge(P<T>& _p, int _i, int _directionx, Fv0 _txbc, Fv1 _tybc, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T tx = _txbc(i + _p.offsetx, j + _p.offsety), ty = _tybc(i + _p.offsetx, j + _p.offsety);
                        if (_directionx == -1) {
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - 4.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + 2.0*tx/3.0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 6)] - (_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + (tx + 3.0*ty)/6.0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 7)] - (_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + (tx - 3.0*ty)/6.0;
                        } else if (_directionx == 1) {
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - 2.0*tx/3.0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 5)] - (_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - (tx - 3.0*ty)/6.0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 8)] - (_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - (tx + 3.0*ty)/6.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of EL set Stress for 2D along y edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetStressAlongYEdge(P<T>& _p, int _j, int _directiony, Fv0 _txbc, Fv1 _tybc, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T tx = _txbc(i + _p.offsetx, j + _p.offsety), ty = _tybc(i + _p.offsetx, j + _p.offsety);
                        if (_directiony == -1) {
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - 4.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + 2.0*ty/3.0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 8)] - (_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + (ty + 3.0*tx)/6.0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 7)] - (_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + (ty - 3.0*tx)/6.0;
                        } else if (_directiony == 1) {
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - 2.0*ty/3.0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 6)] - (_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - (ty + 3.0*tx)/6.0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 5)] - (_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - (ty - 3.0*tx)/6.0;
                        }
                    }
                }
            }
        }
    }

    namespace EL {
        //  Function of Update macro and Collide of EL for 2D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, bool _issave = false) {
            T omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ux, uy, sxx, sxy, syx, syy;
                Macro<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _sxx[idx] = sxx;
                    _sxy[idx] = sxy;
                    _syx[idx] = syx;
                    _syy[idx] = syy;
                }

                //  Collide
                Equilibrium<T, P>(feq, _rho[idx], ux, uy, sxx, sxy, syx, syy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro and Collide of EL with topology optimization for 2D
        template<class T, template<class>class P>
        void MacroExtendedCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, const T *_gamma, bool _issave = false) {
            T omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ux, uy, sxx, sxy, syx, syy;
                Macro<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, _p.f, _gamma[idx], idx);

                //  Save macro if need
                if (_issave) {
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _sxx[idx] = sxx;
                    _sxy[idx] = sxy;
                    _syx[idx] = syx;
                    _syy[idx] = syy;
                }

                //  Collide
                Equilibrium<T, P>(feq, _rho[idx], ux, uy, sxx, sxy, syx, syy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }
    
        //  Function of setting initial condition of EL for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_sxx, const T *_sxy, const T *_syx, const T *_syy) {
            T feq[P<T>::nc];
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _rho[idx], _ux[idx], _uy[idx], _sxx[idx], _sxy[idx], _syx[idx], _syy[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting boundary condition of EL set Stress for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetStress(P<T>& _p, Fv0 _txbc, Fv1 _tybc, Ff _bctype) {
            BoundaryConditionSetStressAlongXEdge(_p, 0, -1, _txbc, _tybc, _bctype);         //  On xmin
            BoundaryConditionSetStressAlongXEdge(_p, _p.lx - 1, 1, _txbc, _tybc, _bctype);  //  On xmax
            BoundaryConditionSetStressAlongYEdge(_p, 0, -1, _txbc, _tybc, _bctype);         //  On ymin
            BoundaryConditionSetStressAlongYEdge(_p, _p.ly - 1, 1, _txbc, _tybc, _bctype);  //  On ymax
        }
    }
}