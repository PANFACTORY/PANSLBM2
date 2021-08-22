//*****************************************************************************
//  Title       :   src/equation/adjointnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>

namespace PANSLBM2 {
    namespace {
        const int SetiU = 1;
        const int SetiRho = 2;
    }

    namespace ANS {
        //  Function of updating macroscopic values of ANS for 2D
        template<class T, class P>
        void Macro(
            T &_ip, T &_iux, T &_iuy, T &_imx, T &_imy, 
            const T *_rho, const T *_ux, const T *_uy, const T *_f, int _idx
        ) {
            _ip = T();
            _iux = T();
            _iuy = T();
            _imx = T();
            _imy = T();
            T uu = _ux[_idx]*_ux[_idx] + _uy[_idx]*_uy[_idx];
            for (int c = 0; c < P::nc; ++c) {
                T ciu = P::cx[c]*_ux[_idx] + P::cy[c]*_uy[_idx];
                T fei = _f[P::IndexF(_idx, c)]*P::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P::cx[c] + 3.0*ciu*P::cx[c] - _ux[_idx]);
                _iuy += fei*(P::cy[c] + 3.0*ciu*P::cy[c] - _uy[_idx]);
                _imx += fei*P::cx[c];
                _imy += fei*P::cy[c];
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class T, class P>
        T Equilibrium(T _ux, T _uy, T _ip, T _iux, T _iuy, int _c) {
            return _ip + 3.0*(_iux*(P::cx[_c] - _ux) + _iuy*(P::cy[_c] - _uy));
        }

        //  Function of applying external force with Brinkman model of ANS for 2D
        template<class T, class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, 
            T _imx, T _imy, T *_f, const T *_alpha, int _idx
        ) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] -= 3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*((P::cx[c] - _ux[_idx])*_imx + (P::cy[c] - _uy[_idx])*_imy);
            }
        }

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of ANS for 2D
        template<class T, class P>
        void Macro_Brinkman_Collide_Stream(
            P& _p, const T *_rho, const T *_ux, const T *_uy, 
            T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _p.f, _alpha, idx);
                    Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _ip[idx] = ip;
                        _iux[idx] = iux;
                        _iuy[idx] = iuy;
                        _imx[idx] = imx;
                        _imy[idx] = imy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                }
            }
        } 
    
        //  Function of setting initial condition of ANS for 2D
        template<class T, class P>
        void InitialCondition(P& _p, const T *_ux, const T *_uy, const T *_ip, const T *_iux, const T *_iuy) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    for (int c = 0; c < P::nc; ++c) {
                        _p.f[P::IndexF(idx, c)] = Equilibrium<T, P>(_ux[idx], _uy[idx], _ip[idx], _iux[idx], _iuy[idx], c);
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for 2D
        template<class T, class P>
        void BoundaryConditionSetiU(P& _p, const T *_uxbc, const T *_uybc, const int *_bctype, T _eps = 1.0) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype[j + _p.offsetxmin] == SetiU) {
                    int idx = _p.Index(0, j), idxbc = j + _p.offsetxmin;
                    T rho0 = (-2.0*_eps + _uxbc[idxbc]*(4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)]) + 3.0*_uybc[idxbc]*(_p.f[P::IndexF(idx, 5)] - _p.f[P::IndexF(idx, 8)]))/(3.0*(1.0 - _uxbc[idxbc]));
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] + rho0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0;
                }

                //  On xmax
                if (_bctype[j + _p.offsetxmax] == SetiU) {
                    int idx = _p.Index(_p.nx - 1, j), idxbc = j + _p.offsetxmax;
                    T rho0 = (-2.0*_eps - _uxbc[idxbc]*(4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)]) + 3.0*_uybc[idxbc]*(_p.f[P::IndexF(idx, 6)] - _p.f[P::IndexF(idx, 7)]))/(3.0*(1.0 + _uxbc[idxbc]));
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] + rho0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype[i + _p.offsetymin] == SetiU) {
                    int idx = _p.Index(i, 0), idxbc = i + _p.offsetymin;
                    T rho0 = (-2.0*_eps + _uxbc[idxbc]*(4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)]) + 3.0*_uxbc[idxbc]*(_p.f[P::IndexF(idx, 5)] - _p.f[P::IndexF(idx, 6)]))/(3.0*(1.0 - _uxbc[idxbc]));
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] + rho0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0;
                }

                //  On ymax
                if (_bctype[i + _p.offsetymax] == SetiU) {
                    int idx = _p.Index(i, _p.ny - 1), idxbc = i + _p.offsetymax;
                    T rho0 = (-2.0*_eps - _uxbc[idxbc]*(4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)]) + 3.0*_uxbc[idxbc]*(_p.f[P::IndexF(idx, 8)] - _p.f[P::IndexF(idx, 7)]))/(3.0*(1.0 + _uxbc[idxbc]));
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] + rho0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0;
                }
            }
        }
    
        //  Function of setting boundary condition of ANS set iRho for 2D
        template<class T, class P>
        void BoundaryConditionSetiRho(P& _p, const int *_bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype[j + _p.offsetxmin] == SetiRho) {
                    int idx = _p.Index(0, j), idxbc = j + _p.offsetxmin;
                    T rho0 = (4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] - rho0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] - rho0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] - rho0;
                }

                //  On xmax
                if (_bctype[j + _p.offsetxmax] == SetiRho) {
                    int idx = _p.Index(_p.nx - 1, j), idxbc = j + _p.offsetxmax;
                    T rho0 = (4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0;
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] - rho0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] - rho0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] - rho0;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype[i + _p.offsetymin] == SetiRho) {
                    int idx = _p.Index(i, 0), idxbc = i + _p.offsetymin;
                    T rho0 = (4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0;
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] - rho0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] - rho0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] - rho0;
                }

                //  On ymax
                if (_bctype[i + _p.offsetymax] == SetiRho) {
                    int idx = _p.Index(i, _p.ny - 1), idxbc = i + _p.offsetymax;
                    T rho0 = (4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] - rho0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] - rho0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] - rho0;
                }
            }
        }
    }
}