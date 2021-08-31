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

        //  Function of updating macroscopic values of ANS for 3D
        template<class T, class P>
        void Macro(
            T &_ip, T &_iux, T &_iuy, T &_iuz, T &_imx, T &_imy, T &_imz, 
            const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_f, int _idx
        ) {
            _ip = T();
            _iux = T();
            _iuy = T();
            _iuz = T();
            _imx = T();
            _imy = T();
            _imz = T();
            T uu = _ux[_idx]*_ux[_idx] + _uy[_idx]*_uy[_idx] + _uz[_idx]*_uz[_idx];
            for (int c = 0; c < P::nc; ++c) {
                T ciu = P::cx[c]*_ux[_idx] + P::cy[c]*_uy[_idx] + P::cz[c]*_uz[_idx];
                T fei = _f[P::IndexF(_idx, c)]*P::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P::cx[c] + 3.0*ciu*P::cx[c] - _ux[_idx]);
                _iuy += fei*(P::cy[c] + 3.0*ciu*P::cy[c] - _uy[_idx]);
                _iuz += fei*(P::cz[c] + 3.0*ciu*P::cz[c] - _uz[_idx]);
                _imx += fei*P::cx[c];
                _imy += fei*P::cy[c];
                _imz += fei*P::cz[c];
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class T, class P>
        T Equilibrium(T _ux, T _uy, T _ip, T _iux, T _iuy, int _c) {
            return _ip + 3.0*(_iux*(P::cx[_c] - _ux) + _iuy*(P::cy[_c] - _uy));
        }

        //  Function of getting equilibrium of ANS for 3D
        template<class T, class P>
        T Equilibrium(T _ux, T _uy, T _uz, T _ip, T _iux, T _iuy, T _iuz, int _c) {
            return _ip + 3.0*(_iux*(P::cx[_c] - _ux) + _iuy*(P::cy[_c] - _uy) + _iuz*(P::cz[_c] - _uz));
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

        //  Function of applying external force with Brinkman model of ANS for 3D
        template<class T, class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, const T *_uz, 
            T _imx, T _imy, T _imz, T *_f, const T *_alpha, int _idx
        ) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] -= 3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*((P::cx[c] - _ux[_idx])*_imx + (P::cy[c] - _uy[_idx])*_imy + (P::cz[c] - _uz[_idx])*_imz);
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

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of ANS for 3D
        template<class T, class P>
        void Macro_Brinkman_Collide_Stream(
            P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, 
            T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T ip, iux, iuy, iuz, imx, imy, imz;
                        Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _p.f, _alpha, idx);
                        Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);

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
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i - P::cx[c], j - P::cy[c], k - P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P::IndexF(idx, c)] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                        }
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

        //  Function of setting initial condition of ANS for 3D
        template<class T, class P>
        void InitialCondition(P& _p, const T *_ux, const T *_uy, const T *_uz, const T *_ip, const T *_iux, const T *_iuy, const T *_iuz) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        for (int c = 0; c < P::nc; ++c) {
                            _p.f[P::IndexF(idx, c)] = Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], _ip[idx], _iux[idx], _iuy[idx], _iuz[idx], c);
                        }
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
    
        //  Function of setting boundary condition of ANS set iU for 3D
        template<class T, class P>
        void BoundaryConditionSetiU(P& _p, const T *_uxbc, const T *_uybc, const T *_uzbc, const int *_bctype, T _eps = 1.0) {
            int idx, idxbc;

            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    //  On xmin
                    idx = _p.Index(0, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmin;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 + _d_uxbc[idxbc]*(8.0*_d_f[P::IndexF(idx, 1)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 12)])
                            + 3.0*_d_uybc[idxbc]*(_d_f[P::IndexF(idx, 7)] - _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 12)])
                            + 3.0*_d_uzbc[idxbc]*(_d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 9)] - _d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 12)])
                        )/(6.0*(1.0 - _d_uxbc[idxbc]));
                        _d_f[P::IndexF(idx, 4)] = _d_f[P::IndexF(idx, 1)] + rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] + rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] + rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] + rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] + rho0;
                    }

                    //  On xmax
                    idx = _p.Index(_p.nx - 1, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmax;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 - _d_uxbc[idxbc]*(8.0*_d_f[P::IndexF(idx, 4)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 13)] + _d_f[P::IndexF(idx, 14)])
                            + 3.0*_d_uybc[idxbc]*(_d_f[P::IndexF(idx, 8)] - _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 13)] - _d_f[P::IndexF(idx, 14)])
                            + 3.0*_d_uzbc[idxbc]*(_d_f[P::IndexF(idx, 8)] - _d_f[P::IndexF(idx, 11)] - _d_f[P::IndexF(idx, 13)] + _d_f[P::IndexF(idx, 14)])
                        )/(6.0*(1.0 + _d_uxbc[idxbc]));
                        _d_f[P::IndexF(idx, 1)] = _d_f[P::IndexF(idx, 4)] + rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] + rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] + rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] + rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] + rho0;
                    }
                }
            }
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    //  On ymin
                    idx = _p.Index(i, 0, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymin;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 + 3.0*_d_uxbc[idxbc]*(_d_f[P::IndexF(idx, 7)] - _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 13)]) 
                            + _d_uybc[idxbc]*(8.0*_d_f[P::IndexF(idx, 2)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 13)])
                            + 3.0*_d_uzbc[idxbc]*(_d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] - _d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 13)])
                        )/(6.0*(1.0 - _d_uybc[idxbc]));
                        _d_f[P::IndexF(idx, 5)] = _d_f[P::IndexF(idx, 2)] + rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] + rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] + rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] + rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] + rho0;
                    }

                    //  On ymax
                    idx = _p.Index(i, _p.ny - 1, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymax;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 + 3.0*_d_uxbc[idxbc]*(_d_f[P::IndexF(idx, 9)] - _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] - _d_f[P::IndexF(idx, 14)]) 
                            - _d_uybc[idxbc]*(8.0*_d_f[P::IndexF(idx, 5)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 14)])
                            + 3.0*_d_uzbc[idxbc]*(_d_f[P::IndexF(idx, 9)] - _d_f[P::IndexF(idx, 11)] - _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 14)])
                        )/(6.0*(1.0 + _d_uybc[idxbc]));
                        _d_f[P::IndexF(idx, 2)] = _d_f[P::IndexF(idx, 5)] + rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] + rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] + rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] + rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] + rho0;
                    }
                }
            }
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    //  On zmin
                    idx = _p.Index(i, j, 0);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmin;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 + 3.0*_d_uxbc[idxbc]*(_d_f[P::IndexF(idx, 7)] - _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 9)] - _d_f[P::IndexF(idx, 14)])
                            + 3.0*_d_uybc[idxbc]*(_d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] - _d_f[P::IndexF(idx, 9)] - _d_f[P::IndexF(idx, 14)])
                            + _d_uzbc[idxbc]*(8.0*_d_f[P::IndexF(idx, 3)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 14)])
                        )/(6.0*(1.0 - _d_uzbc[idxbc]));
                        _d_f[P::IndexF(idx, 6)] = _d_f[P::IndexF(idx, 3)] + rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] + rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] + rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] + rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] + rho0;
                    }

                    //  On zmax
                    idx = _p.Index(i, j, _p.nz - 1);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmax;
                    if (_d_bctype[idxbc] == INLET) {
                        T rho0 = (-4.0 + 3.0*_d_uxbc[idxbc]*(_d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] - _d_f[P::IndexF(idx, 13)])
                            + 3.0*_d_uybc[idxbc]*(_d_f[P::IndexF(idx, 10)] - _d_f[P::IndexF(idx, 11)] - _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 13)])
                            - _d_uzbc[idxbc]*(8.0*_d_f[P::IndexF(idx, 6)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 13)])
                        )/(6.0*(1.0 + _d_uzbc[idxbc]));
                        _d_f[P::IndexF(idx, 3)] = _d_f[P::IndexF(idx, 6)] + rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] + rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] + rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] + rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] + rho0;
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iRho for 2D
        template<class T, class P>
        void BoundaryConditionSetiRho2D(P& _p, const int *_bctype) {
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
    
        //  Function of setting boundary condition of ANS set iRho for 3D
        template<class T, class P>
        void BoundaryConditionSetiRho3D(P& _p, const int *_bctype) {
            int idx, idxbc;

            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    //  On xmin
                    idx = _p.Index(0, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmin;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 1)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 12)])/6.0;
                        _d_f[P::IndexF(idx, 4)] = _d_f[P::IndexF(idx, 1)] - rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] - rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] - rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] - rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] - rho0;
                    }

                    //  On xmax
                    idx = _p.Index(_p.nx - 1, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmax;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 4)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 13)] + _d_f[P::IndexF(idx, 14)])/6.0;
                        _d_f[P::IndexF(idx, 1)] = _d_f[P::IndexF(idx, 4)] - rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] - rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] - rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] - rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] - rho0;
                    }
                }
            }
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    //  On ymin
                    idx = _p.Index(i, 0, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymin;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 2)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 13)])/6.0;
                        _d_f[P::IndexF(idx, 5)] = _d_f[P::IndexF(idx, 2)] - rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] - rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] - rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] - rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] - rho0;
                    }

                    //  On ymax
                    idx = _p.Index(i, _p.ny - 1, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymax;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 5)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 14)])/6.0;
                        _d_f[P::IndexF(idx, 2)] = _d_f[P::IndexF(idx, 5)] - rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] - rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] - rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] - rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] - rho0;
                    }
                }
            }
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    //  On zmin
                    idx = _p.Index(i, j, 0);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmin;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 3)] + _d_f[P::IndexF(idx, 7)] + _d_f[P::IndexF(idx, 8)] + _d_f[P::IndexF(idx, 9)] + _d_f[P::IndexF(idx, 14)])/6.0;
                        _d_f[P::IndexF(idx, 6)] = _d_f[P::IndexF(idx, 3)] - rho0;
                        _d_f[P::IndexF(idx, 10)] = _d_f[P::IndexF(idx, 14)] - rho0;
                        _d_f[P::IndexF(idx, 11)] = _d_f[P::IndexF(idx, 7)] - rho0;
                        _d_f[P::IndexF(idx, 12)] = _d_f[P::IndexF(idx, 8)] - rho0;
                        _d_f[P::IndexF(idx, 13)] = _d_f[P::IndexF(idx, 9)] - rho0;
                    }

                    //  On zmax
                    idx = _p.Index(i, j, _p.nz - 1);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmax;
                    if (_d_bctype[idxbc] == OUTLET) {
                        T rho0 = (8.0*_d_f[P::IndexF(idx, 6)] + _d_f[P::IndexF(idx, 10)] + _d_f[P::IndexF(idx, 11)] + _d_f[P::IndexF(idx, 12)] + _d_f[P::IndexF(idx, 13)])/6.0;
                        _d_f[P::IndexF(idx, 3)] = _d_f[P::IndexF(idx, 6)] - rho0;
                        _d_f[P::IndexF(idx, 7)] = _d_f[P::IndexF(idx, 11)] - rho0;
                        _d_f[P::IndexF(idx, 8)] = _d_f[P::IndexF(idx, 12)] - rho0;
                        _d_f[P::IndexF(idx, 9)] = _d_f[P::IndexF(idx, 13)] - rho0;
                        _d_f[P::IndexF(idx, 14)] = _d_f[P::IndexF(idx, 10)] - rho0;
                    }
                }
            }
        }
    }
}