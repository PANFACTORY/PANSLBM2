//*****************************************************************************
//  Title       :   src/equation/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace NS {
        //  Function of updating macroscopic values of NS for 2D
        template<class T, template<class>class P>
        void Macro(T &_rho, T &_ux, T &_uy, const T *_f, int _idx) {
            _rho = T();
            _ux = T();
            _uy = T();
            for (int c = 0; c < P<T>::nc; ++c) {
                _rho += _f[P<T>::IndexF(_idx, c)];
                _ux += P<T>::cx[c]*_f[P<T>::IndexF(_idx, c)];
                _uy += P<T>::cy[c]*_f[P<T>::IndexF(_idx, c)];
            }
            _ux /= _rho;
            _uy /= _rho;
        }

        //  Function of updating macroscopic values of NS for 3D
        template<class T, template<class>class P>
        void Macro(T &_rho, T &_ux, T &_uy, T &_uz, const T *_f, int _idx) {
            _rho = T();
            _ux = T();
            _uy = T();
            _uz = T();
            for (int c = 0; c < P<T>::nc; ++c) {
                _rho += _f[P<T>::IndexF(_idx, c)];
                _ux += P<T>::cx[c]*_f[P<T>::IndexF(_idx, c)];
                _uy += P<T>::cy[c]*_f[P<T>::IndexF(_idx, c)];
                _uz += P<T>::cz[c]*_f[P<T>::IndexF(_idx, c)];
            }
            _ux /= _rho;
            _uy /= _rho;
            _uz /= _rho;
        }

        //  Function of getting equilibrium of NS for 2D
        template<class T, template<class>class P>
        T Equilibrium(T _rho, T _ux, T _uy, int _c) {
            T uu = _ux*_ux + _uy*_uy;
            T ciu = P<T>::cx[_c]*_ux + P<T>::cy[_c]*_uy;
            return P<T>::ei[_c]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
        }

        //  Function of getting equilibrium of NS for 3D
        template<class T, template<class>class P>
        T Equilibrium(T _rho, T _ux, T _uy, T _uz, int _c) {
            T uu = _ux*_ux + _uy*_uy + _uz*_uz;
            T ciu = P<T>::cx[_c]*_ux + P<T>::cy[_c]*_uy + P<T>::cz[_c]*_uz;
            return P<T>::ei[_c]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
        }

        //  Function of applying external force of NS with Brinkman model for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _alpha, T *_f, int _idx) {
            for (int c = 0; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= 3.0*_alpha/(1.0 + _alpha/_rho)*P<T>::ei[c]*(P<T>::cx[c]*_ux + P<T>::cy[c]*_uy);
            }
        }

        //  Function of applying external force of NS with Brinkman model for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _uz, T _alpha, T *_f, int _idx) {
            for (int c = 0; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= 3.0*_alpha/(1.0 + _alpha/_rho)*P<T>::ei[c]*(P<T>::cx[c]*_ux + P<T>::cy[c]*_uy + P<T>::cz[c]*_uz);
            }
        }

        //  Function of Update macro, Collide and Stream of NS for 2D
        template<class T, template<class>class P>
        void Macro_Collide_Stream(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    Macro<T, P<T> >(rho, ux, uy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P<T>::nc; ++c) {
                        int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c]);
                        _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P<T> >(rho, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of NS for 3D
        template<class T, template<class>class P>
        void Macro_Collide_Stream(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        Macro<T, P<T> >(rho, ux, uy, uz, _p.f, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P<T>::nc; ++c) {
                            int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c], k + P<T>::cz[c]);
                            _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P<T> >(rho, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of NS for 2D
        template<class T, template<class>class P>
        void Macro_Brinkman_Collide_Stream(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, const T *_alpha, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    Macro<T, P<T> >(rho, ux, uy, _p.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P<T> >(rho, ux, uy, _alpha[idx], _p.f, idx);
                    Macro<T, P<T> >(rho, ux, uy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P<T>::nc; ++c) {
                        int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c]);
                        _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P<T> >(rho, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model), Collide and Stream of NS for 3D
        template<class T, template<class>class P>
        void Macro_Brinkman_Collide_Stream(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, const T *_alpha, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        Macro<T, P<T> >(rho, ux, uy, uz, _p.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P<T> >(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        Macro<T, P<T> >(rho, ux, uy, uz, _p.f, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P<T>::nc; ++c) {
                            int idxstream = _p.Index(i + P<T>::cx[c], j + P<T>::cy[c], k + P<T>::cz[c]);
                            _p.fnext[P<T>::IndexF(idxstream, c)] = (1.0 - omega)*_p.f[P<T>::IndexF(idx, c)] + omega*Equilibrium<T, P<T> >(rho, ux, uy, uz, c);
                        }
                    }
                    
                }
            }
        }

        //  Function of setting initial condition of NS for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    for (int c = 0; c < P<T>::nc; ++c) {
                        _p.f[P<T>::IndexF(idx, c)] = Equilibrium<T, P<T> >(_rho[idx], _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of setting initial condition of NS for 3D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        for (int c = 0; c < P<T>::nc; ++c) {
                            _p.f[P<T>::IndexF(idx, c)] = Equilibrium<T, P<T> >(_rho[idx], _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Ff _bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(0, j);
                    T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety));
                    T mx = rho0*_uxbc(0 + _p.offsetx, j + _p.offsety)/6.0;
                    T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uybc(0 + _p.offsetx, j + _p.offsety));
                    _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                    _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                    _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                }

                //  On xmax
                if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(_p.nx - 1, j);
                    T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                    T mx = rho0*_uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety)/6.0;
                    T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uybc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                    _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                    _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                    _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                    int idx = _p.Index(i, 0);
                    T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 - _uybc(i + _p.offsetx, 0 + _p.offsety));
                    T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*_uxbc(i + _p.offsetx, 0 + _p.offsety));
                    T my = rho0*_uybc(i + _p.offsetx, 0 + _p.offsety)/6.0;
                    _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                    _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                    _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
                }

                //  On ymax
                if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                    int idx = _p.Index(i, _p.ny - 1);
                    T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/(1.0 + _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                    T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*_uxbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                    T my = rho0*_uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)/6.0;
                    _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                    _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                    _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                }
            }
        }
    
        //  Function of setting boundary condition of NS set U for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    //  On xmin
                    if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(0, j, k);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - _uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                        T mx = rho0*_uxbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)/12.0;
                        T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                        T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                        _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
                    }

                    //  On xmax
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(_p.nx - 1, j, k);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/(1.0 + _uxbc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz));
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
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    //  On ymin
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(i, 0, k);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - _uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                        T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                        T my = rho0*_uybc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)/12.0;
                        T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*_uzbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                        _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
                    }

                    //  On ymax
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(i, _p.ny - 1, k);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 + _uybc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz));
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
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    //  On zmin
                    if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)) {
                        int idx = _p.Index(i, j, 0);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 - _uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                        T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*_uxbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                        T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*_uybc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                        T mz = rho0*_uzbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)/12.0;
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                        _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
                    }

                    //  On zmax
                    if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)) {
                        int idx = _p.Index(i, j, _p.nz - 1);
                        T rho0 = (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 + _uzbc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz));
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

        //  Function of setting boundary condition of NS set rho for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Ff _bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                //  On xmin
                if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(0, j);
                    T ux0 = 1.0 - (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/_rhobc(0 + _p.offsetx, j + _p.offsety);
                    T mx = _rhobc(0 + _p.offsetx, j + _p.offsety)*ux0/6.0;
                    T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(0 + _p.offsetx, j + _p.offsety)*_usbc(0 + _p.offsetx, j + _p.offsety));
                    _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                    _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                    _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                }

                //  On xmax
                if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) {
                    int idx = _p.Index(_p.nx - 1, j);
                    T ux0 = -1.0 + (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/_rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety);
                    T mx = _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*ux0/6.0;
                    T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety)*_usbc((_p.nx - 1) + _p.offsetx, j + _p.offsety));
                    _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                    _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                    _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                }
            }

            for (int i = 0; i < _p.nx; ++i) {
                //  On ymin
                if (_bctype(i + _p.offsetx, 0 + _p.offsety)) {
                    int idx = _p.Index(i, 0);
                    T uy0 = 1.0 - (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/_rhobc(i + _p.offsetx, 0 + _p.offsety);
                    T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - _rhobc(i + _p.offsetx, 0 + _p.offsety)*_usbc(i + _p.offsetx, 0 + _p.offsety));
                    T my = _rhobc(i + _p.offsetx, 0 + _p.offsety)*uy0/6.0;
                    _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                    _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                    _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
                }

                //  On ymax
                if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {
                    int idx = _p.Index(i, _p.ny - 1);
                    T uy0 = -1.0 + (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/_rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety);
                    T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*_usbc(i + _p.offsetx, (_p.ny - 1) + _p.offsety));
                    T my = _rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety)*uy0/6.0;
                    _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                    _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                    _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                }
            }
        }
    
        //  Function of setting boundary condition of NS set rho for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    //  On xmin
                    if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(0, j, k);
                        T ux0 = 1.0 - (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                        T mx = _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*ux0/12.0;
                        T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_usbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                        T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)*_utbc(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz));
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                        _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
                    }

                    //  On xmax
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(_p.nx - 1, j, k);
                        T ux0 = -1.0 + (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/_rhobc((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz);
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
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    //  On ymin
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(i, 0, k);
                        T uy0 = 1.0 - (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz);
                        T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*_utbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                        T my = _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*uy0/12.0;
                        T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - _rhobc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)*_usbc(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz));
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                        _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                        _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
                    }

                    //  On ymax
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)) {
                        int idx = _p.Index(i, _p.ny - 1, k);
                        T uy0 = -1.0 + (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/_rhobc(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz);
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
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    //  On zmin
                    if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)) {
                        int idx = _p.Index(i, j, 0);
                        T uz0 = 1.0 - (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/_rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz);
                        T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*_usbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                        T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*_utbc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz));
                        T mz = _rhobc(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)*uz0/12.0;
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                        _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                        _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
                    }

                    //  On zmax
                    if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)) {
                        int idx = _p.Index(i, j, _p.nz - 1);
                        T uz0 = -1.0 + (_p.f[P<T>::IndexF(idx, 0)] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/_rhobc(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz);
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