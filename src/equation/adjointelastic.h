//*****************************************************************************
//  Title       :   src/equation/adjointelastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace AEL {
        //  Function of updating macroscopic values of AEL for 2D
        template<class T, class P>
        void Macro(T &_irho, T &_imx, T &_imy, T &_isxx, T &_isxy, T &_isyx, T &_isyy, const T *_f, int _idx) {
            _irho = T();
            _imx = T();
            _imy = T();
            _isxx = T();
            _isxy = T();
            _isyx = T();
            _isyy = T();
            for (int c = 0; c < P::nc; ++c) {
                _irho += P::ei[c]*_f[P::IndexF(_idx, c)];
                _imx += P::ei[c]*P::cx[c]*_f[P::IndexF(_idx, c)];
                _imy += P::ei[c]*P::cy[c]*_f[P::IndexF(_idx, c)];
                _isxx += P::ei[c]*P::cx[c]*P::cx[c]*_f[P::IndexF(_idx, c)];
                _isxy += P::ei[c]*P::cx[c]*P::cy[c]*_f[P::IndexF(_idx, c)];
                _isyx += P::ei[c]*P::cy[c]*P::cx[c]*_f[P::IndexF(_idx, c)];
                _isyy += P::ei[c]*P::cy[c]*P::cy[c]*_f[P::IndexF(_idx, c)];
            }
        }

        //  Function of getting equilibrium of AEL for 2D
        template<class T, class P>
        T Equilibrium(T _irho, T _imx, T _imy, T _isxx, T _isxy, T _isyx, T _isyy, T _gamma, int _c) {
            T imc = _imx*P::cx[_c] + _imy*P::cy[_c];
            T cisc = P::cx[_c]*_isxx*P::cx[_c] + P::cx[_c]*_isxy*P::cy[_c] + P::cy[_c]*_isyx*P::cx[_c] + P::cy[_c]*_isyy*P::cy[_c];
            T irhocc = _irho*(P::cx[_c]*P::cx[_c] + P::cy[_c]*P::cy[_c]);
            return 3.0*imc + 4.5*_gamma*cisc - 1.5*_gamma*irhocc;
        }

        //  Function of Update macro, Collide and Stream of AEL for 2D
        template<class T, class P>
        void Macro_Collide_Stream(
            P& _p, T *_irho, T *_imx, T *_imy, T *_isxx, T *_isxy, T *_isyx, T *_isyy, T _tau, const T *_gamma, bool _issave = false
        ) {
            T omega = 1/_tau;
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T irho, imx, imy, isxx, isxy, isyx, isyy;
                    Macro<T, P>(irho, imx, imy, isxx, isxy, isyx, isyy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _irho[idx] = irho;
                        _imx[idx] = imx;
                        _imy[idx] = imy;
                        _isxx[idx] = isxx;
                        _isxy[idx] = isxy;
                        _isyx[idx] = isyx;
                        _isyy[idx] = isyy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(irho, imx, imy, isxx, isxy, isyx, isyy, _gamma[idx], c);
                    }
                }
            }
        }

        //  Function of setting initial condition of AEL for 2D
        template<class T, class P>
        void InitialCondition(P& _p, const T *_irho, const T *_imx, const T *_imy, const T *_isxx, const T *_isxy, const T *_isyx, const T *_isyy, const T *_gamma) {
            for (int idx = 0; idx < _p.nxy; ++idx) {
                for (int c = 0; c < P::nc; ++c) {
                    _p.f[P::IndexF(idx, c)] = Equilibrium<T, P>(_irho[idx], _imx[idx], _imy[idx], _isxx[idx], _isxy[idx], _isyx[idx], _isyy[idx], _gamma[idx], c);
                }
            }
        }

        //  Function of setting boundary condition of AEL set iStress for 2D
        template<class T, class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetiStress(P& _p, Fv0 _txbc, Fv1 _tybc, const T *_rho, Ff _bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(0, j);
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] - (4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 + 2.0*_txbc(0 + _p.offsetx, j + _p.offsety)/_rho[idx];
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 5)] - (4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 + 2.0*(_txbc(0 + _p.offsetx, j + _p.offsety) - _tybc(0 + _p.offsetx, j + _p.offsety))/_rho[idx];
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 8)] - (4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 + 2.0*(_txbc(0 + _p.offsetx, j + _p.offsety) + _tybc(0 + _p.offsetx, j + _p.offsety))/_rho[idx];
                }

                //  On xmax
                if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(_p.nx - 1, j);
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] - (4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 - 2.0*_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)/_rho[idx];
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 6)] - (4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 - 2.0*(_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety) + _tybc((_p.nx - 1) + _p.offsetx, j + _p.offsety))/_rho[idx];
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 7)] - (4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 - 2.0*(_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety) - _tybc((_p.nx - 1) + _p.offsetx, j + _p.offsety))/_rho[idx];
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                    int idx = _p.Index(i, 0);
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] - (4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 + 2.0*_tybc(i + _p.offsetx, 0 + _p.offsety)/_rho[idx];
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 6)] - (4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 + 2.0*(_tybc(i + _p.offsetx, 0 + _p.offsety) + _txbc(i + _p.offsetx, 0 + _p.offsety))/_rho[idx];
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 5)] - (4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 + 2.0*(_tybc(i + _p.offsetx, 0 + _p.offsety) - _txbc(i + _p.offsetx, 0 + _p.offsety))/_rho[idx];
                }

                //  On ymax
                if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                    int idx = _p.Index(i, _p.ny - 1);
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] - (4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 - 2.0*_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)/_rho[idx];
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 8)] - (4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 - 2.0*(_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety) + _txbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety))/_rho[idx];
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 7)] - (4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 - 2.0*(_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety) - _txbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety))/_rho[idx];
                }
            }
        } 
    }
}