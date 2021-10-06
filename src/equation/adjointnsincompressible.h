//*****************************************************************************
//  Title       :   src/equation/adjointnsincompressible.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/10/07
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once

namespace PANSLBM2 {
    namespace ANSin {
        //  Function of updating macroscopic values of ANS for 2D
        template<class T, template<class>class P>
        void Macro(T &_ip, T &_iux, T &_iuy, T &_imx, T &_imy, T _rho, T _ux, T _uy, const T *_f0, const T *_f, int _idx) {
            _ip = _f0[_idx]*P<T>::ei[0];
            _iux = -_f0[_idx]*P<T>::ei[0]*_ux;
            _iuy = -_f0[_idx]*P<T>::ei[0]*_uy;
            _imx = T();
            _imy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy;
                T fei = _f[P<T>::IndexF(_idx, c)]*P<T>::ei[c];
                _ip += fei;
                _iux += fei*(P<T>::cx[c] + 3.0*ciu*P<T>::cx[c] - _ux);
                _iuy += fei*(P<T>::cy[c] + 3.0*ciu*P<T>::cy[c] - _uy);
                _imx += fei*P<T>::cx[c];
                _imy += fei*P<T>::cy[c];
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _ip, T _iux, T _iuy) {
            for (int c = 0; c < P<T>::nc; ++c) {
                _feq[c] = _ip + 3.0*(_iux*P<T>::cx[c] + _iuy*P<T>::cy[c]);
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _imx, T _imy, T *_f, T _alpha, int _idx) {
            T coef = 3.0*_alpha/(1.0 + _alpha);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= coef*(P<T>::cx[c]*_imx + P<T>::cy[c]*_imy);
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of ANS for 2D
        template<class T, template<class>class P>
        void MacroBrinkmanCollide(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, 
            T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(imx, imy, _p.f, _alpha[idx], idx);
                Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                }

                //  Collide
                Equilibrium<T, P>(feq, ip, iux, iuy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of setting initial condition of ANS for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_ux, const T *_uy, const T *_ip, const T *_iux, const T *_iuy) {
            T feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _ip[idx], _iux[idx], _iuy[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D2Q9
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetU(P<T>& _p, Ff _bctype, T _eps = T()) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - 2.0*_eps/3.0;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - 2.0*_eps/3.0;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - 2.0*_eps/3.0;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - 2.0*_eps/3.0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - 2.0*_eps/3.0;
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
    }
}