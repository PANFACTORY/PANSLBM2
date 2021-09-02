//*****************************************************************************
//  Title       :   src/equation/elastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace {
    const int SetS = 1;
}

namespace PANSLBM2 {
    namespace EL {
        //  Function of updating macroscopic values of EL for 2D
        template<class T, class P>
        void Macro(T _rho, T &_ux, T &_uy, T &_sxx, T &_sxy, T &_syx, T &_syy, const T *_f, int _idx) {
            _ux = T();
            _uy = T();
            _sxx = T();
            _sxy = T();
            _syx = T();
            _syy = T();
            for (int c = 0; c < P::nc; ++c) {
                _ux += P::cx[c]*_f[P::IndexF(_idx, c)];
                _uy += P::cy[c]*_f[P::IndexF(_idx, c)];
                _sxx -= P::cx[c]*P::cx[c]*_f[P::IndexF(_idx, c)];
                _sxy -= P::cx[c]*P::cy[c]*_f[P::IndexF(_idx, c)];
                _syx -= P::cy[c]*P::cx[c]*_f[P::IndexF(_idx, c)];
                _syy -= P::cy[c]*P::cy[c]*_f[P::IndexF(_idx, c)];
            }
            _ux /= _rho;
            _uy /= _rho;
        }

        //  Function of updating macroscopic values of EL with topology optimization for 2D
        template<class T, class P>
        void Macro(T _rho, T &_ux, T &_uy, T &_sxx, T &_sxy, T &_syx, T &_syy, const T *_f, const T *_gamma, int _idx) {
            Macro<T, P>(_rho, _ux, _uy, _sxx, _sxy, _syx, _syy, _f, _idx);
            _sxx *= _gamma[_idx];
            _sxy *= _gamma[_idx];
            _syx *= _gamma[_idx];
            _syy *= _gamma[_idx];
        }
    
        //  Function of getting equilibrium of EL for 2D
        template<class T, class P>
        T Equilibrium(T _rho, T _ux, T _uy, T _sxx, T _sxy, T _syx, T _syy, int _c) {
            T ciu = P::cx[_c]*_ux + P::cy[_c]*_uy;
            T cisci = P::cx[_c]*_sxx*P::cx[_c] + P::cx[_c]*_sxy*P::cy[_c] + P::cy[_c]*_syx*P::cx[_c] + P::cy[_c]*_syy*P::cy[_c];
            T trss = _sxx + _syy;
            return P::ei[_c]*(3.0*_rho*ciu - 4.5*cisci + 1.5*trss);
        }

        //  Function of Update macro, Collide and Stream of EL for 2D
        template<class T, class P>
        void Macro_Collide_Stream(P& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, bool _issave = false) {
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
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of EL with topology optimization for 2D
        template<class T, class P>
        void MacroExtended_Collide_Stream(P& _p, T *_rho, T *_ux, T *_uy, T *_sxx, T *_sxy, T *_syx, T *_syy, T _tau, const T *_gamma, bool _issave = false) {
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
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(_rho[idx], ux, uy, sxx, sxy, syx, syy, c);
                    }
                }
            }
        }
    
        //  Function of setting initial condition of EL for 2D
        template<class T, class P>
        void InitialCondition(P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_sxx, const T *_sxy, const T *_syx, const T *_syy) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    for (int c = 0; c < P::nc; ++c) {
                        _p.f[P::IndexF(idx, c)] = Equilibrium<T, P>(_rho[idx], _ux[idx], _uy[idx], _sxx[idx], _sxy[idx], _syx[idx], _syy[idx], c);
                    }
                }
            }
        }

        //  Function of setting boundary condition of EL set Stress for 2D
        template<class T, class P>
        void BoundaryConditionSetStress(P& _p, const T *_txbc, const T *_tybc, const int *_bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype[j + _p.offsetxmin] == SetS) {
                    int idx = _p.Index(0, j), idxbc = j + _p.offsetxmin;
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] - 4.0*(_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 + 2.0*_txbc[idxbc]/3.0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 6)] - (_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 + (_txbc[idxbc] + 3.0*_tybc[idxbc])/6.0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 7)] - (_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0 + (_txbc[idxbc] - 3.0*_tybc[idxbc])/6.0;
                }

                //  On xmax
                if (_bctype[j + _p.offsetxmax] == SetS) {
                    int idx = _p.Index(_p.nx - 1, j), idxbc = j + _p.offsetxmax;
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] - 4.0*(_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 - 2.0*_txbc[idxbc]/3.0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 5)] - (_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 - (_txbc[idxbc] - 3.0*_tybc[idxbc])/6.0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 8)] - (_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0 - (_txbc[idxbc] + 3.0*_tybc[idxbc])/6.0;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype[i + _p.offsetymin] == SetS) {
                    int idx = _p.Index(i, 0), idxbc = i + _p.offsetymin;
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] - 4.0*(_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 + 2.0*_tybc[idxbc]/3.0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 8)] - (_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 + (_tybc[idxbc] + 3.0*_txbc[idxbc])/6.0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 7)] - (_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0 + (_tybc[idxbc] - 3.0*_txbc[idxbc])/6.0;
                }

                //  On ymax
                if (_bctype[i + _p.offsetymax] == SetS) {
                    int idx = _p.Index(i, _p.ny - 1), idxbc = i + _p.offsetymax;
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] - 4.0*(_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 - 2.0*_tybc[idxbc]/3.0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 6)] - (_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 - (_tybc[idxbc] + 3.0*_txbc[idxbc])/6.0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 5)] - (_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0 - (_tybc[idxbc] - 3.0*_txbc[idxbc])/6.0;
                }
            }
        }
    }
}