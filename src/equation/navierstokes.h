//*****************************************************************************
//  Title       :   src/equation/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace {
        const int INLET = 1;
        const int OUTLET = 2;
    }

    namespace NS {
        //  Function of updating macroscopic values of NS for 2D
        template<class T, class P>
        void Macro(T &_rho, T &_ux, T &_uy, const T *_f, int _idx) {
            _rho = T();
            _ux = T();
            _uy = T();
            for (int c = 0; c < P::nc; ++c) {
                _rho += _f[P::IndexF(_idx, c)];
                _ux += P::cx[c]*_f[P::IndexF(_idx, c)];
                _uy += P::cy[c]*_f[P::IndexF(_idx, c)];
            }
            _ux /= _rho;
            _uy /= _rho;
        }

        //  Function of getting equilibrium of NS for 2D
        template<class T, class P>
        T Equilibrium(T _rho, T _ux, T _uy, int _c) {
            T uu = _ux*_ux + _uy*_uy;
            T ciu = P::cx[_c]*_ux + P::cy[_c]*_uy;
            return P::ei[_c]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
        }

        //  Function of Update macro, Collide and Stream of NS for 2D
        template<class T, class P>
        void Macro_Collide_Stream(P& _p, T *_rho, T *_ux, T *_uy, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    Macro<T, P>(rho, ux, uy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.IndexStream(i, j, c);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(rho, ux, uy, c);
                    }
                }
            }
        }

        //  Function of setting initial condition of NS for 2D
        template<class T, class P>
        void InitialCondition(P& _p, const T *_rho, const T *_ux, const T *_uy) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    for (int c = 0; c < P::nc; ++c) {
                        _p.f[P::IndexF(idx, c)] = Equilibrium<T, P>(_rho[idx], _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for 2D
        template<class T, class P>
        void BoundaryConditionSetU(P& _p, const T *_uxbc, const T *_uybc, const int *_bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype[j + _p.offsetxmin] == INLET) {
                    int idx = _p.Index(0, j), idxbc = j + _p.offsetxmin;
                    T rho0 = (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 4)] + 2.0*(_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)]))/(1.0 - _uxbc[idxbc]);
                    T mx = rho0*_uxbc[idxbc]/6.0;
                    T my = 0.5*(_p.f[P::IndexF(idx, 2)] - _p.f[P::IndexF(idx, 4)] - rho0*_uybc[idxbc]);
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] + 4.0*mx;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + mx - my;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + mx + my;
                }

                //  On xmax
                if (_bctype[j + _p.offsetxmax] == INLET) {
                    int idx = _p.Index(_p.nx - 1, j), idxbc = j + _p.offsetxmax;
                    T rho0 = (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 4)] + 2.0*(_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)]))/(1.0 + _uxbc[idxbc]);
                    T mx = rho0*_uxbc[idxbc]/6.0;
                    T my = 0.5*(_p.f[P::IndexF(idx, 2)] - _p.f[P::IndexF(idx, 4)] - rho0*_uybc[idxbc]);
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] - 4.0*mx;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] - mx + my;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] - mx - my;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype[i + _p.offsetymin] == INLET) {
                    int idx = _p.Index(i, 0), idxbc = i + _p.offsetymin;
                    T rho0 = (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 3)] + 2.0*(_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)]))/(1.0 - _uybc[idxbc]);
                    T mx = 0.5*(_p.f[P::IndexF(idx, 1)] - _p.f[P::IndexF(idx, 3)] - rho0*_uxbc[idxbc]);
                    T my = rho0*_uybc[idxbc]/6.0;
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] + 4.0*my;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] - mx + my;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + mx + my;
                }

                //  On ymax
                if (_bctype[i + _p.offsetymax] == INLET) {
                    int idx = _p.Index(i, _p.ny - 1), idxbc = i + _p.offsetymax;
                    T rho0 = (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 3)] + 2.0*(_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)]))/(1.0 + _uybc[idxbc]);
                    T mx = 0.5*(_p.f[P::IndexF(idx, 1)] - _p.f[P::IndexF(idx, 3)] - rho0*_uxbc[idxbc]);
                    T my = rho0*_uybc[idxbc]/6.0;
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] - 4.0*my;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + mx - my;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] - mx - my;
                }
            }
        }
    
        //  Function of setting boundary condition of NS set rho for 2D
        template<class T, class P>
        void BoundaryConditionSetRho(P& _p, const T *_rhobc, const T *_usbc, const int *_bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype[j + _p.offsetxmin] == OUTLET) {
                    int idx = _p.Index(0, j), idxbc = j + _p.offsetxmin;
                    T ux0 = 1.0 - (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 4)] + 2.0*(_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)]))/_rhobc[idxbc];
                    T mx = _rhobc[idxbc]*ux0/6.0;
                    T my = 0.5*(_p.f[P::IndexF(idx, 2)] - _p.f[P::IndexF(idx, 4)] - _rhobc[idxbc]*_usbc[idxbc]);
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] + 4.0*mx;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + mx - my;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + mx + my;
                }

                //  On xmax
                if (_bctype[j + _p.offsetxmax] == OUTLET) {
                    int idx = _p.Index(_p.nx - 1, j), idxbc = j + _p.offsetxmax;
                    T ux0 = -1.0 + (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 4)] + 2.0*(_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)]))/_rhobc[idxbc];
                    T mx = _rhobc[idxbc]*ux0/6.0;
                    T my = 0.5*(_p.f[P::IndexF(idx, 2)] - _p.f[P::IndexF(idx, 4)] - _rhobc[idxbc]*_usbc[idxbc]);
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] - 4.0*mx;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] - mx - my;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] - mx + my;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype[i + _p.offsetymin] == OUTLET) {
                    int idx = _p.Index(i, 0), idxbc = i + _p.offsetymin;
                    T uy0 = 1.0 - (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 3)] + 2.0*(_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)]))/_rhobc[idxbc];
                    T mx = 0.5*(_p.f[P::IndexF(idx, 1)] - _p.f[P::IndexF(idx, 3)] - _rhobc[idxbc]*_usbc[idxbc]);
                    T my = _rhobc[idxbc]*uy0/6.0;
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] + 4.0*my;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] - mx + my;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + mx + my;
                }

                //  On ymax
                if (_bctype[i + _p.offsetymax] == OUTLET) {
                    int idx = _p.Index(i, _p.ny - 1), idxbc = i + _p.offsetymax;
                    T uy0 = -1.0 + (_p.f[P::IndexF(idx, 0)] + _p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 3)] + 2.0*(_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)]))/_rhobc[idxbc];
                    T mx = 0.5*(_p.f[P::IndexF(idx, 1)] - _p.f[P::IndexF(idx, 3)] - _rhobc[idxbc]*_usbc[idxbc]);
                    T my = _rhobc[idxbc]*uy0/6.0;
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] - 4.0*my;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + mx - my;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] - mx - my;
                }
            }
        }
    }
}