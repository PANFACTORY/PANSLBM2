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
    namespace ANS {
        //  Function of updating macroscopic values of ANS for 2D
        template<class T, template<class>class P>
        void Macro(
            T &_ip, T &_iux, T &_iuy, T &_imx, T &_imy, 
            const T *_rho, const T *_ux, const T *_uy, const T *_f0, const T *_f, int _idx
        ) {
            T uu = _ux[_idx]*_ux[_idx] + _uy[_idx]*_uy[_idx];
            _ip = _f0[_idx]*P<T>::ei[0]*(1.0 - 1.5*uu);
            _iux = -_f0[_idx]*P<T>::ei[0]*_ux[_idx];
            _iuy = -_f0[_idx]*P<T>::ei[0]*_uy[_idx];;
            _imx = T();
            _imy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux[_idx] + P<T>::cy[c]*_uy[_idx];
                T fei = _f[P<T>::IndexF(_idx, c)]*P<T>::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P<T>::cx[c] + 3.0*ciu*P<T>::cx[c] - _ux[_idx]);
                _iuy += fei*(P<T>::cy[c] + 3.0*ciu*P<T>::cy[c] - _uy[_idx]);
                _imx += fei*P<T>::cx[c];
                _imy += fei*P<T>::cy[c];
            }
        }

        //  Function of updating macroscopic values of ANS for 3D
        template<class T, template<class>class P>
        void Macro(
            T &_ip, T &_iux, T &_iuy, T &_iuz, T &_imx, T &_imy, T &_imz, 
            const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_f0, const T *_f, int _idx
        ) {
            T uu = _ux[_idx]*_ux[_idx] + _uy[_idx]*_uy[_idx] + _uz[_idx]*_uz[_idx];
            _ip = _f0[_idx]*P<T>::ei[0]*(1.0 - 1.5*uu);;
            _iux = -_f0[_idx]*P<T>::ei[0]*_ux[_idx];
            _iuy = -_f0[_idx]*P<T>::ei[0]*_uy[_idx];;
            _iuz = -_f0[_idx]*P<T>::ei[0]*_uz[_idx];;
            _imx = T();
            _imy = T();
            _imz = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux[_idx] + P<T>::cy[c]*_uy[_idx] + P<T>::cz[c]*_uz[_idx];
                T fei = _f[P<T>::IndexF(_idx, c)]*P<T>::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P<T>::cx[c] + 3.0*ciu*P<T>::cx[c] - _ux[_idx]);
                _iuy += fei*(P<T>::cy[c] + 3.0*ciu*P<T>::cy[c] - _uy[_idx]);
                _iuz += fei*(P<T>::cz[c] + 3.0*ciu*P<T>::cz[c] - _uz[_idx]);
                _imx += fei*P<T>::cx[c];
                _imy += fei*P<T>::cy[c];
                _imz += fei*P<T>::cz[c];
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class T, template<class>class P>
        T Equilibrium(T _ux, T _uy, T _ip, T _iux, T _iuy, int _c) {
            return _ip + 3.0*(_iux*(P<T>::cx[_c] - _ux) + _iuy*(P<T>::cy[_c] - _uy));
        }

        //  Function of getting equilibrium of ANS for 3D
        template<class T, template<class>class P>
        T Equilibrium(T _ux, T _uy, T _uz, T _ip, T _iux, T _iuy, T _iuz, int _c) {
            return _ip + 3.0*(_iux*(P<T>::cx[_c] - _ux) + _iuy*(P<T>::cy[_c] - _uy) + _iuz*(P<T>::cz[_c] - _uz));
        }

        //  Function of applying external force with Brinkman model of ANS for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, 
            T _imx, T _imy, T *_f0, T *_f, const T *_alpha, int _idx
        ) {
            _f0[_idx] -= -3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*(_ux[_idx]*_imx + _uy[_idx]*_imy);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= 3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*((P<T>::cx[c] - _ux[_idx])*_imx + (P<T>::cy[c] - _uy[_idx])*_imy);
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, const T *_uz, 
            T _imx, T _imy, T _imz, T *_f0, T *_f, const T *_alpha, int _idx
        ) {
            _f0[_idx] -= -3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*(_ux[_idx]*_imx + _uy[_idx]*_imy + _uz[_idx]*_imz);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= 3.0*_alpha[_idx]/(_rho[_idx] + _alpha[_idx])*((P<T>::cx[c] - _ux[_idx])*_imx + (P<T>::cy[c] - _uy[_idx])*_imy + (P<T>::cz[c] - _uz[_idx])*_imz);
            }
        }

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of ANS for 2D
        template<class T, template<class>class P>
        void MacroBrinkmanCollideStream(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, 
            T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T ip, iux, iuy, imx, imy;
                Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _p.f0, _p.f, _alpha, idx);
                Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omega)*_p.f0[idx] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of ANS for 3D
        template<class T, template<class>class P>
        void MacroBrinkmanCollideStream(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, 
            T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                    //  Update macro
                    T ip, iux, iuy, iuz, imx, imy, imz;
                    Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _p.f0, _p.f, _alpha, idx);
                    Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);

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
                    _p.f0[idx] = (1.0 - omega)*_p.f0[idx] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, 0);
                    for (int c = 1; c < P<T>::nc; ++c) {
                        _p.f[P<T>::IndexF(idx, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                    }
                }
            }
        } 
    
        //  Function of setting initial condition of ANS for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_ux, const T *_uy, const T *_ip, const T *_iux, const T *_iuy) {
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _p.f0[idx] = Equilibrium<T, P>(_ux[idx], _uy[idx], _ip[idx], _iux[idx], _iuy[idx], 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = Equilibrium<T, P>(_ux[idx], _uy[idx], _ip[idx], _iux[idx], _iuy[idx], c);
                }
            }
        }

        //  Function of setting initial condition of ANS for 3D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_ux, const T *_uy, const T *_uz, const T *_ip, const T *_iux, const T *_iuy, const T *_iuz) {
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _p.f0[idx] = Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], _ip[idx], _iux[idx], _iuy[idx], _iuz[idx], 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], _ip[idx], _iux[idx], _iuy[idx], _iuz[idx], c);
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D2Q9
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Ff _bctype, T _eps = T()) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        T rho0 = (-2.0*_eps + _uxbc(0 + _p.offsetx, j + _p.offsety)*(4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]) + 3.0*_uybc(0 + _p.offsetx, j + _p.offsety)*(_p.f[P<T>::IndexF(idx, 5)] - _p.f[P<T>::IndexF(idx, 8)]))/(3.0*(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety)));
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] + rho0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T rho0 = (-2.0*_eps - _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*(4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]) + 3.0*_uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*(_p.f[P<T>::IndexF(idx, 6)] - _p.f[P<T>::IndexF(idx, 7)]))/(3.0*(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)));
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + rho0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T rho0 = (-2.0*_eps + _uxbc(i + _p.offsetx, 0 + _p.offsety)*(4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]) + 3.0*_uxbc(i + _p.offsetx, 0 + _p.offsety)*(_p.f[P<T>::IndexF(idx, 5)] - _p.f[P<T>::IndexF(idx, 6)]))/(3.0*(1.0 - _uxbc(i + _p.offsetx, 0 + _p.offsety)));
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] + rho0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        T rho0 = (-2.0*_eps - _uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*(4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]) + 3.0*_uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 7)]))/(3.0*(1.0 + _uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)));
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + rho0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition of ANS set iU for D3Q15
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void iBoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype, T _eps = T()) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(0, j, k);
                            T rho0 = (-4.0*_eps + _uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])
                                + 3.0*_uybc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 12)])
                                + 3.0*_uzbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 12)])
                            )/(6.0*(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] + rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0;
                        }
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(_p.nx - 1, j, k);
                            T rho0 = (-4.0*_eps - _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])
                                + 3.0*_uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] - _p.f[P<T>::IndexF(idx, 14)])
                                + 3.0*_uzbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])
                            )/(6.0*(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                        }
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, 0, k);
                            T rho0 = (-4.0*_eps + 3.0*_uxbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 13)]) 
                                + _uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])
                                + 3.0*_uzbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 13)])
                            )/(6.0*(1.0 - _uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] + rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0;
                        }
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, _p.ny - 1, k);
                            T rho0 = (-4.0*_eps + 3.0*_uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] - _p.f[P<T>::IndexF(idx, 14)]) 
                                - _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])
                                + 3.0*_uzbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])
                            )/(6.0*(1.0 + _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                        }
                    }
                }
            }
            //  On zmin
            if (_p.PEz == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)) {
                            int idx = _p.Index(i, j, 0);
                            T rho0 = (-4.0*_eps + 3.0*_uxbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 14)])
                                + 3.0*_uybc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 14)])
                                + _uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])
                            )/(6.0*(1.0 - _uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] + rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                        }
                    }
                }
            }
            //  On zmax
            if (_p.PEz == _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)) {
                            int idx = _p.Index(i, j, _p.nz - 1);
                            T rho0 = (-4.0*_eps + 3.0*_uxbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] - _p.f[P<T>::IndexF(idx, 13)])
                                + 3.0*_uybc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*(_p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])
                                - _uzbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*(8.0*_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])
                            )/(6.0*(1.0 + _uzbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)));
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iRho for D2Q9
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRho2D(P<T>& _p, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - rho0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - rho0;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0;
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - rho0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - rho0;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0;
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - rho0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - rho0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - rho0;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - rho0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition of ANS set iRho for D3Q15
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRho3D(P<T>& _p, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(0, j, k);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])/6.0;
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - rho0;
                        }
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(_p.nx - 1, j, k);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] - rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                        }
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, 0, k);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - rho0;
                        }
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, _p.ny - 1, k);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] - rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - rho0;
                        }
                    }
                }
            }
            //  On zmin
            if (_p.PEz == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)) {
                            int idx = _p.Index(i, j, 0);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] - rho0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - rho0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - rho0;
                        }
                    }
                }
            }
            //  On zmax
            if (_p.PEz == _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)) {
                            int idx = _p.Index(i, j, _p.nz - 1);
                            T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] - rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - rho0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - rho0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - rho0;
                        }
                    }
                }
            }
        }
    }
}