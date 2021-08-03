//*****************************************************************************
//  Title       :   src/equation/adjointadvection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>
#include "adjointnavierstokes.h"

namespace PANSLBM2 {
    namespace {
        const int SetiT = 1;
        const int SetiQ = 2;
    }

    namespace AAD {
        //  Function of updating macroscopic values of AAD for 2D
        template<class T, class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, const T *_g, int _idx) {
            _item = T();
            _iqx = T();
            _iqy = T();
            for (int c = 0; c <Q::nc; ++c) {
                _item += Q::ei[c]*_g[Q::IndexF(_idx, c)];
                _iqx += Q::ei[c]*Q::cx[c]*_g[Q::IndexF(_idx, c)];
                _iqy += Q::ei[c]*Q::cy[c]*_g[Q::IndexF(_idx, c)];
            }
        }

        //  Function of getting equilibrium of AAD for 2D
        template<class T, class Q>
        T Equilibrium(T _item, T _iqx, T _iqy, T _ux, T _uy, int _c) {
            return _item + 3.0*(_ux*_iqx + _uy*_iqy);
        }

        //  Function of applying external force from advection of AAD for 2D
        template<class T, class P>
        void ExternalForceFromAdvection(T _iqx, T _iqy, const T *_rho, const T *_ux, const T *_uy, T *_f, const T *_tem, T _omegag, int _idx) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] += 3.0*_tem[_idx]*_omegag*((P::cx[c] - _ux[_idx])*_iqx + (P::cy[c] - _uy[_idx])*_iqy)/_rho[_idx];
            }
        }

        //  Function of applying external force with heat exchange of AAD for 2D
        template<class T, class Q>
        void ExternalForceHeatExchange(T _item, T *_g, const T *_beta, int _idx) {
            for (int c = 0; c < Q::nc; ++c) {
                _g[Q::IndexF(_idx, c)] -= _beta[_idx]*(1.0 + _item)/(1.0 + _beta[_idx]);
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange), Collide and Stream of AAD for 2D
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_HeatExchange(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceFromAdvection<T, P>(iqx, iqy, _rho, _ux, _uy, _p.f, _tem, omegag, idx);
                    ANS::ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    ExternalForceHeatExchange<T, Q>(item, _q.f, _beta, idx);
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _ip[idx] = ip;
                        _iux[idx] = iux;
                        _iuy[idx] = iuy;
                        _imx[idx] = imx;
                        _imy[idx] = imy;
                        _item[idx] = item;
                        _iqx[idx] = iqx;
                        _iqy[idx] = iqy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.IndexiStream(i, j, c);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.IndexiStream(i, j, c);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }
    
        //  Function of setting initial condition of AAD for 2D
        template<class T, class Q>
        void InitialCondition(
            Q& _q, const T *_ux, const T *_uy, const T *_item, const T *_iqx, const T *_iqy
        ) {
            for (int idx = 0; idx < _q.nxy; ++idx) {
                for (int c = 0; c < Q::nc; ++c) {
                    _q.f[Q::IndexF(idx, c)] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for 2D
        template<class T, class Q>
        void BoundaryConditionSetiT(Q& _q, const T *_ux, const T *_uy, const int *_bctype) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == SetiT) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T rho0 = -(4.0*(1.0 + 3.0*_ux[idx])*_q.f[Q::IndexF(idx, 1)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 5)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 8)])/(6.0*(1.0 + 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 3)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == SetiT) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T rho0 = -(4.0*(1.0 - 3.0*_ux[idx])*_q.f[Q::IndexF(idx, 3)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 6)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 7)])/(6.0*(1.0 - 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 1)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == SetiT) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T rho0 = -(4.0*(1.0 + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 2)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 5)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 6)])/(6.0*(1.0 + 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 4)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == SetiT) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T rho0 = -(4.0*(1.0 - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 4)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 7)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 8)])/(6.0*(1.0 - 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 2)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                }
            }
        }
    
        //  Function of setting boundary condition set iQ of AAD for 2D
        template<class T, class Q>
        void BoundaryConditionSetiQ(Q& _q, const T *_ux, const T *_uy, const int *_bctype, T _eps = T()) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == SetiQ) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T rho0 = (
                        (1.0 + 3.0*_ux[idx])*(4.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])
                        + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])
                        - 12.0*_eps
                    )/(6.0*(1.0 - 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 3)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == SetiQ) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T rho0 = (
                        (1.0 - 3.0*_ux[idx])*(4.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])
                        + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])
                        - 12.0*_eps
                    )/(6.0*(1.0 + 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 1)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == SetiQ) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T rho0 = (
                        (1.0 + 3.0*_uy[idx])*(4.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])
                        + 3.0*_ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])
                        - 12.0*_eps
                    )/(6.0*(1.0 - 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 4)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == SetiQ) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T rho0 = (
                        (1.0 - 3.0*_uy[idx])*(4.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])
                        + 3.0*_ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])
                        - 12.0*_eps
                    )/(6.0*(1.0 + 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 2)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                }
            }
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for 2D
        template<class T, class P, class Q>
        void BoundaryConditionSetiRho(P& _p, Q& _q, const T *_rho, const T *_ux, const T *_uy, const T *_tem, const int *_bctypef, const int *_bctypeg, T _eps = T()) {
            int idx, idxbc; 
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                idx = _q.Index(0, j);
                idxbc = j + _q.offsetxmin;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])/(2.0*(1.0 + 3.0*_ux[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])/3.0
                            + _uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])/2.0
                        )/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                }

                //  On xmax
                idx = _q.Index(_q.nx - 1, j);
                idxbc = j + _q.offsetxmax;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_ux[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])/3.0
                             + _uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])/2.0
                        )/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                idx = _q.Index(i, 0);
                idxbc = i + _q.offsetymin;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])/(2.0*(1.0 + 3.0*_uy[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])/3.0
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])/2.0
                        )/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                }

                //  On ymax
                idx = _q.Index(i, _q.ny - 1);
                idxbc = i + _q.offsetymax;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_uy[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])/3.0 
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])/2.0
                        )/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                }
            }
        }
    }
}