//*****************************************************************************
//  Title       :   src/equation/nsincompressible.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/10/07
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once

namespace PANSLBM2 {
    namespace NSin {
        //  Function of updating macroscopic values of NS for 2D
        template<class T, template<class>class P>
        void Macro(T &_rho, T &_ux, T &_uy, const T *_f0, const T *_f, int _idx) {
            _rho = _f0[_idx];
            _ux = T();
            _uy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T f = _f[P<T>::IndexF(_idx, c)];
                _rho += f;
                _ux += P<T>::cx[c]*f;
                _uy += P<T>::cy[c]*f;
            }
        }

        //  Function of getting equilibrium of NS for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _rho, T _ux, T _uy) {
            T rhouu = _rho - 1.5*(_ux*_ux + _uy*_uy);
            for (int c = 0; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy;
                _feq[c] = P<T>::ei[c]*(3.0*ciu + 4.5*ciu*ciu + rhouu);
            }
        }

        //  Function of applying external force of NS with Brinkman model for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _alpha, T *_f, int _idx) {
            T coef = 3.0*_alpha*_rho/(_rho + _alpha);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= coef*P<T>::ei[c]*(P<T>::cx[c]*_ux + P<T>::cy[c]*_uy);
            }
        }

        //  Function of Update macro and Collide of NS for 2D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                }

                //  Collide
                Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of NS for 2D
        template<class T, template<class>class P>
        void MacroBrinkmanCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, const T *_alpha, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                }

                //  Collide
                Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of setting initial condition of NS for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy) {
            T feq[P<T>::nc];
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _rho[idx], _ux[idx], _uy[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting boundary condition of NS set U for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        T ux = _uxbc(0 + _p.offsetx, j + _p.offsety), uy = _uybc(0 + _p.offsetx, j + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 2.0*ux/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + ux/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + ux/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T ux = _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety), uy = _uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 2.0*ux/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - ux/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - ux/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T ux = _uxbc(i + _p.offsetx, 0 + _p.offsety), uy = _uybc(i + _p.offsetx, 0 + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 2.0*uy/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + uy/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + uy/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/(1.0 + _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                        T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*_uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                        T my = rho0*_uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)/6.0;

                        T ux = _uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety), uy = _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 2.0*uy/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - uy/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - uy/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set rho for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(0, j);
                        T ux = _rhobc(0 + _p.offsetx, j + _p.offsety) - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]));
                        T uy = _usbc(0 + _p.offsetx, j + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 2.0*ux/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + ux/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + ux/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T ux = -_rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety) + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]));
                        T uy = _usbc((_p.nx - 1) + _p.offsetx, j + _p.offsety);
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 2.0*ux/3.0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - ux/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - ux/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - uy);
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T ux = _usbc(i + _p.offsetx, 0 + _p.offsety);
                        T uy = _rhobc(i + _p.offsetx, 0 + _p.offsety) - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]));
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 2.0*uy/3.0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + uy/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + uy/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        T ux = _usbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety);
                        T uy = -_rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety) + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]));
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 2.0*uy/3.0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - uy/6.0 + 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - uy/6.0 - 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - ux);
                    }
                }
            }
            
        }
    }
}