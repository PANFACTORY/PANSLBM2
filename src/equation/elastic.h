//*****************************************************************************
//  Title       :   src/equation/elastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

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
            for (int c = 0; c < P<T>::nc; ++c) {
                _ux += P<T>::cx[c]*_f[P<T>::IndexF(_idx, c)];
                _uy += P<T>::cy[c]*_f[P<T>::IndexF(_idx, c)];
                _sxx -= P<T>::cx[c]*P<T>::cx[c]*_f[P<T>::IndexF(_idx, c)];
                _sxy -= P<T>::cx[c]*P<T>::cy[c]*_f[P<T>::IndexF(_idx, c)];
                _syx -= P<T>::cy[c]*P<T>::cx[c]*_f[P<T>::IndexF(_idx, c)];
                _syy -= P<T>::cy[c]*P<T>::cy[c]*_f[P<T>::IndexF(_idx, c)];
            }
            _ux /= _rho;
            _uy /= _rho;
        }

        //  Function of updating macroscopic values of EL with topology optimization for 2D
        template<class T, template<class>class P>
        void Macro(T _rho, T &_ux, T &_uy, T &_sxx, T &_sxy, T &_syx, T &_syy, const T *_f, const T *_gamma, int _idx) {
            Macro<T, P>(_rho, _ux, _uy, _sxx, _sxy, _syx, _syy, _f, _idx);
            _sxx *= _gamma[_idx];
            _sxy *= _gamma[_idx];
            _syx *= _gamma[_idx];
            _syy *= _gamma[_idx];
        }
    
        //  Function of getting equilibrium of EL for 2D
        template<class T, template<class>class P>
        T Equilibrium(T _rho, T _ux, T _uy, T _sxx, T _sxy, T _syx, T _syy, int _c) {
            T ciu = P<T>::cx[_c]*_ux + P<T>::cy[_c]*_uy;
            T cisci = P<T>::cx[_c]*_sxx*P<T>::cx[_c] + P<T>::cx[_c]*_sxy*P<T>::cy[_c] + P<T>::cy[_c]*_syx*P<T>::cx[_c] + P<T>::cy[_c]*_syy*P<T>::cy[_c];
            T trss = _sxx + _syy;
            return P<T>::ei[_c]*(3.0*_rho*ciu - 4.5*cisci + 1.5*trss);
        }

        //  Function of Update macro, Collide and Stream of EL for 2D
        template<class T, template<class>class P>
        void MacroCollideStream(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, bool _issave = false) {
            T omega = 1/_tau;
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

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

                    //  Collide and stream
                    for (int c = 0; c < P<T>::nc; ++c) {
                        int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c]);
                        _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of EL with topology optimization for 2D
        template<class T, template<class>class P>
        void MacroExtendedCollideStream(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, const T *_gamma, bool _issave = false) {
            T omega = 1/_tau;
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ux, uy, sxx, sxy, syx, syy;
                    Macro<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, _p.f, _gamma, idx);

                    //  Save macro if need
                    if (_issave) {
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _sxx[idx] = sxx;
                        _sxy[idx] = sxy;
                        _syx[idx] = syx;
                        _syy[idx] = syy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P<T>::nc; ++c) {
                        int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c]);
                        _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, c);
                    }
                }
            }
        }
    
        //  Function of setting initial condition of EL for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_sxx, const T *_sxy, const T *_syx, const T *_syy) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    for (int c = 0; c < P<T>::nc; ++c) {
                        _p.f[P<T>::IndexF(idx, c)] = Equilibrium<T, P>(_rho[idx], _ux[idx], _uy[idx], _sxx[idx], _sxy[idx], _syx[idx], _syy[idx], c);
                    }
                }
            }
        }

        //  Function of setting boundary condition of EL set Stress for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetStress(P<T>& _p, Fv0 _txbc, Fv1 _tybc, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - 4.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + 2.0*_txbc(0 + _p.offsetx, j + _p.offsety)/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 6)] - (_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + (_txbc(0 + _p.offsetx, j + _p.offsety) + 3.0*_tybc(0 + _p.offsetx, j + _p.offsety))/6.0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 7)] - (_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 + (_txbc(0 + _p.offsetx, j + _p.offsety) - 3.0*_tybc(0 + _p.offsetx, j + _p.offsety))/6.0;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - 2.0*_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)/3.0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 5)] - (_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - (_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety) - 3.0*_tybc((_p.nx - 1) + _p.offsetx, j + _p.offsety))/6.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 8)] - (_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - (_txbc((_p.nx - 1) + _p.offsetx, j + _p.offsety) + 3.0*_tybc((_p.nx - 1) + _p.offsetx, j + _p.offsety))/6.0;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - 4.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + 2.0*_tybc(i + _p.offsetx, 0 + _p.offsety)/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 8)] - (_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + (_tybc(i + _p.offsetx, 0 + _p.offsety) + 3.0*_txbc(i + _p.offsetx, 0 + _p.offsety))/6.0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 7)] - (_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + (_tybc(i + _p.offsetx, 0 + _p.offsety) - 3.0*_txbc(i + _p.offsetx, 0 + _p.offsety))/6.0;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - 2.0*_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 6)] - (_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - (_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety) + 3.0*_txbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety))/6.0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 5)] - (_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 - (_tybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety) - 3.0*_txbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety))/6.0;
                    }
                }
            }
        }
    }
}