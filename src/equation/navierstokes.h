//*****************************************************************************
//  Title       :   src/equation/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/navierstokes_avx.h"
#endif

namespace PANSLBM2 {
    namespace NS {
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
            T invrho = 1.0/_rho;
            _ux *= invrho;
            _uy *= invrho;
        }

        //  Function of updating macroscopic values of NS for 3D
        template<class T, template<class>class P>
        void Macro(T &_rho, T &_ux, T &_uy, T &_uz, const T *_f0, const T *_f, int _idx) {
            _rho = _f0[_idx];
            _ux = T();
            _uy = T();
            _uz = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T f = _f[P<T>::IndexF(_idx, c)];
                _rho += f;
                _ux += P<T>::cx[c]*f;
                _uy += P<T>::cy[c]*f;
                _uz += P<T>::cz[c]*f;
            }
            T invrho = 1.0/_rho;
            _ux *= invrho;
            _uy *= invrho;
            _uz *= invrho;
        }

        //  Function of getting equilibrium of NS for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _rho, T _ux, T _uy) {
            T uu = 1.0 - 1.5*(_ux*_ux + _uy*_uy);
            for (int c = 0; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy;
                _feq[c] = P<T>::ei[c]*_rho*(3.0*ciu + 4.5*ciu*ciu + uu);
            }
        }

        //  Function of getting equilibrium of NS for 3D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _rho, T _ux, T _uy, T _uz) {
            T uu = 1.0 - 1.5*(_ux*_ux + _uy*_uy + _uz*_uz);
            for (int c = 0; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy + P<T>::cz[c]*_uz;
                _feq[c] = P<T>::ei[c]*_rho*(3.0*ciu + 4.5*ciu*ciu + uu);
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

        //  Function of applying external force of NS with Brinkman model for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _uz, T _alpha, T *_f, int _idx) {
            T coef = 3.0*_alpha*_rho/(_rho + _alpha);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= coef*P<T>::ei[c]*(P<T>::cx[c]*_ux + P<T>::cy[c]*_uy + P<T>::cz[c]*_uz);
            }
        }
    }

#ifndef _USE_AVX_DEFINES
    namespace NS {
        //  Function of Update macro and Collide of NS for 2D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
#pragma omp parallel for
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

        //  Function of Update macro and Collide of NS for 3D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                }

                //  Collide
                Equilibrium<T, P>(feq, rho, ux, uy, uz);
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
#pragma omp parallel for
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

        //  Function of Update macro, External force(Brinkman model) and Collide of NS for 3D
        template<class T, template<class>class P>
        void MacroBrinkmanCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, const T *_alpha, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                }

                //  Collide
                Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }
    }   
#endif 
    
    //  Initial conditions and boundary conditions
    namespace NS {
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

        //  Function of setting initial condition of NS for 3D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz) {
            T feq[P<T>::nc];
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _rho[idx], _ux[idx], _uy[idx], _uz[idx]);
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
                        T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety));
                        T mx = rho0*_uxbc(0 + _p.offsetx, j + _p.offsety)/6.0;
                        T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uybc(0 + _p.offsetx, j + _p.offsety));
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                        T mx = rho0*_uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)/6.0;
                        T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 - _uybc(i + _p.offsetx, 0 + _p.offsety));
                        T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*_uxbc(i + _p.offsetx, 0 + _p.offsety));
                        T my = rho0*_uybc(i + _p.offsetx, 0 + _p.offsety)/6.0;
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
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
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition of NS set U for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(0, j, k);
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mx = rho0*_uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)/12.0;
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
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
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mx = rho0*_uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)/12.0;
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - 8.0*mx;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - mx + my - mz;
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
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - _uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                            T my = rho0*_uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)/12.0;
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
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
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 + _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
                            T my = rho0*_uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)/12.0;
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - 8.0*my;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx - my - mz;
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
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 - _uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                            T mz = rho0*_uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)/12.0;
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
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
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 + _uzbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
                            T mz = rho0*_uzbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)/12.0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] - 8.0*mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx - my - mz;
                        }
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
                        T ux0 = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/_rhobc(0 + _p.offsetx, j + _p.offsety);
                        T mx = _rhobc(0 + _p.offsetx, j + _p.offsety)*ux0/6.0;
                        T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(0 + _p.offsetx, j + _p.offsety)*_usbc(0 + _p.offsetx, j + _p.offsety));
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(_p.nx - 1, j);
                        T ux0 = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/_rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety);
                        T mx = _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*ux0/6.0;
                        T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*_usbc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                        int idx = _p.Index(i, 0);
                        T uy0 = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/_rhobc(i + _p.offsetx, 0 + _p.offsety);
                        T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - _rhobc(i + _p.offsetx, 0 + _p.offsety)*_usbc(i + _p.offsetx, 0 + _p.offsety));
                        T my = _rhobc(i + _p.offsetx, 0 + _p.offsety)*uy0/6.0;
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                        int idx = _p.Index(i, _p.ny - 1);
                        T uy0 = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/_rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety);
                        T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*_usbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                        T my = _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*uy0/6.0;
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                    }
                }
            }
            
        }
    
        //  Function of setting boundary condition of NS set rho for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(0, j, k);
                            T ux0 = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            T mx = _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*ux0/12.0;
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_usbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_utbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
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
                            T ux0 = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/_rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            T mx = _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*ux0/12.0;
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_usbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_utbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - 8.0*mx;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - mx + my - mz;
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
                            T uy0 = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz);
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*_utbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                            T my = _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*uy0/12.0;
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*_usbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
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
                            T uy0 = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/_rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz);
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*_utbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
                            T my = _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*uy0/12.0;
                            T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)*_usbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - 8.0*my;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx - my - mz;
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
                            T uz0 = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/_rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz);
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*_usbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*_utbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                            T mz = _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*uz0/12.0;
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
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
                            T uz0 = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz);
                            T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*_usbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
                            T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*_utbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
                            T mz = _rhobc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)*uz0/12.0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] - 8.0*mz;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx - my - mz;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx + my - mz;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx - my - mz;
                        }
                    }
                }
            }
        }
    }
}