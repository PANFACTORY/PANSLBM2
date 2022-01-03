//*****************************************************************************
//  Title       :   src/equation/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
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

        //  Function of setting boundary condition of NS set U for D2Q9 along x edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetUAlongXEdge(P<T>& _p, int _i, int _directionx, Fv0 _uxbc, Fv1 _uybc, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T ux = _uxbc(i + _p.offsetx, j + _p.offsety), uy = _uybc(i + _p.offsetx, j + _p.offsety);
                        if (_directionx == -1) {
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/(1.0 - ux);
                            T mx = rho0*ux/6.0;
                            T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*uy);
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                        } else if (_directionx == 1) {
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 + ux);
                            T mx = rho0*ux/6.0;
                            T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*uy);
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for D2Q9 along y edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetUAlongYEdge(P<T>& _p, int _j, int _directiony, Fv0 _uxbc, Fv1 _uybc, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T ux = _uxbc(i + _p.offsetx, j + _p.offsety), uy = _uybc(i + _p.offsetx, j + _p.offsety);
                        if (_directiony == -1) {
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/(1.0 - uy);
                            T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*ux);
                            T my = rho0*uy/6.0;
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
                        } else if (_directiony == 1) {
                            T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/(1.0 + uy);
                            T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho0*ux);
                            T my = rho0*uy/6.0;
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for D3Q15 along x face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetUAlongXFace(P<T>& _p, int _i, int _directionx, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionx == -1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - ux);
                                T mx = rho0*ux/12.0;
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*uy);
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*uz);
                                _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
                            } else if (_directionx == 1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/(1.0 + ux);
                                T mx = rho0*ux/12.0;
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*uy);
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*uz);
                                _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - 8.0*mx;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - mx + my - mz;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for D3Q15 along y face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetUAlongYFace(P<T>& _p, int _j, int _directiony, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directiony == -1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 - uy);
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*ux);
                                T my = rho0*uy/12.0;
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*uz);
                                _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
                            } else if (_directiony == 1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 + uy);
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*ux);
                                T my = rho0*uy/12.0;
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho0*uz);
                                _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - 8.0*my;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx - my - mz;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set U for D3Q15 along z face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetUAlongZFace(P<T>& _p, int _k, int _directionz, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            int k = _k - _p.offsetz;
            if (0 <= k && k < _p.nz) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionz == -1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/(1.0 - uz);
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*ux);
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*uy);
                                T mz = rho0*uz/12.0;
                                _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
                            } else if (_directionz == 1) {
                                T rho0 = (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/(1.0 + uz);
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho0*ux);
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho0*uy);
                                T mz = rho0*uz/12.0;
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

        //  Function of setting boundary condition of NS set rho for D2Q9 along x edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetRhoAlongXEdge(P<T>& _p, int _i, int _directionx, Fv0 _rhobc, Fv1 _usbc, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T rho = _rhobc(i + _p.offsetx, j + _p.offsety), uy = _usbc(i + _p.offsetx, j + _p.offsety);
                        if (_directionx == -1) {
                            T ux = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]))/rho;
                            T mx = rho*ux/6.0;
                            T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho*uy);
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + 4.0*mx;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + mx + my;
                        } else if (_directionx == 1) {
                            T ux = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]))/rho;
                            T mx = rho*ux/6.0;
                            T my = 0.5*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 4)] - rho*uy);
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - 4.0*mx;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - mx + my;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set rho for D2Q9 along y edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetRhoAlongYEdge(P<T>& _p, int _j, int _directiony, Fv0 _rhobc, Fv1 _usbc, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T rho = _rhobc(i + _p.offsetx, j + _p.offsety), ux = _usbc(i + _p.offsetx, j + _p.offsety);
                        if (_directiony == -1) {
                            T uy = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]))/rho;
                            T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho*ux);
                            T my = rho*uy/6.0;
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + 4.0*my;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my;
                        } else if (_directiony == 1) {
                            T uy = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]))/rho;
                            T mx = 0.5*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 3)] - rho*ux);
                            T my = rho*uy/6.0;
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - 4.0*my;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + mx - my;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - mx - my;
                        }
                    }
                }
            }
        }
    
        //  Function of setting boundary condition of NS set rho for D3Q15 along x face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRhoAlongXFace(P<T>& _p, int _i, int _directionx, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T rho = _rhobc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _usbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _utbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionx == -1) {
                                T ux = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)]))/rho;
                                T mx = rho*ux/12.0;
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho*uy);
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho*uz);
                                _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + 8.0*mx;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + mx + my + mz;
                            } else if (_directionx == 1) {
                                T ux = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)]))/rho;
                                T mx = rho*ux/12.0;
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho*uy);
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho*uz);
                                _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - 8.0*mx;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - mx + my - mz;
                            }
                        } 
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set rho for D3Q15 along y face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRhoAlongYFace(P<T>& _p, int _j, int _directiony, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T rho = _rhobc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _usbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), ux = _utbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directiony == -1) {
                                T uy = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)]))/rho;
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho*ux);
                                T my = rho*uy/12.0;
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho*uz);
                                _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + 8.0*my;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx + my - mz;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + mx + my + mz;
                            } else if (_directiony == 1) {
                                T uy = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 6)] + 2.0*(_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)]))/rho;
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho*ux);
                                T my = rho*uy/12.0;
                                T mz = 0.25*(_p.f[P<T>::IndexF(idx, 3)] - _p.f[P<T>::IndexF(idx, 6)] - rho*uz);
                                _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - 8.0*my;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx - my - mz;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx - my - mz;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of NS set rho for D3Q15 along z face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRhoAlongZFace(P<T>& _p, int _k, int _directionz, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            int k = _k - _p.offsetz;
            if (0 <= k && k < _p.nz) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T rho = _rhobc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), ux = _usbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _utbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionz == -1) {
                                T uz = 1.0 - (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)]))/rho;
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho*ux);
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho*uy);
                                T mz = rho*uz/12.0;
                                _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + 8.0*mz;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] - mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + mx - my + mz;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - mx + my + mz;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + mx + my + mz;
                            } else if (_directionz == 1) {
                                T uz = -1.0 + (_p.f0[idx] + _p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 5)] + 2.0*(_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)]))/rho;
                                T mx = 0.25*(_p.f[P<T>::IndexF(idx, 1)] - _p.f[P<T>::IndexF(idx, 4)] - rho*ux);
                                T my = 0.25*(_p.f[P<T>::IndexF(idx, 2)] - _p.f[P<T>::IndexF(idx, 5)] - rho*uy);
                                T mz = rho*uz/12.0;
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

    namespace NS {
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

        //  Function of Update macro and Collide of NS for 3D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
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

        //  Function of Update macro, External force(Brinkman model) and Collide of NS for 3D
        template<class T, template<class>class P>
        void MacroBrinkmanCollide(P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity, const T *_alpha, bool _issave = false) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
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
    
        //  Function of Update macro and Collide of NS with Level-Set-Method for 2D
        template<class T, template<class>class P>
        void MacroCollideLSM(P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity, const T *_chi, bool _issave = false) {
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
                Equilibrium<T, P>(feq, rho, ux*_chi[idx], uy*_chi[idx]);
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
            BoundaryConditionSetUAlongXEdge(_p, 0, -1, _uxbc, _uybc, _bctype);          //  On xmin
            BoundaryConditionSetUAlongXEdge(_p, _p.lx - 1, 1, _uxbc, _uybc, _bctype);   //  On xmax
            BoundaryConditionSetUAlongYEdge(_p, 0, -1, _uxbc, _uybc, _bctype);          //  On ymin
            BoundaryConditionSetUAlongYEdge(_p, _p.ly - 1, 1, _uxbc, _uybc, _bctype);   //  On ymax
        }
    
        //  Function of setting boundary condition of NS set U for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype) {
            BoundaryConditionSetUAlongXFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype);          //  On xmin
            BoundaryConditionSetUAlongXFace(_p, _p.lx - 1, 1, _uxbc, _uybc, _uzbc, _bctype);   //  On xmax
            BoundaryConditionSetUAlongYFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype);          //  On ymin
            BoundaryConditionSetUAlongYFace(_p, _p.ly - 1, 1, _uxbc, _uybc, _uzbc, _bctype);   //  On ymax
            BoundaryConditionSetUAlongZFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype);          //  On zmin
            BoundaryConditionSetUAlongZFace(_p, _p.lz - 1, 1, _uxbc, _uybc, _uzbc, _bctype);   //  On zmax
        }

        //  Function of setting boundary condition of NS set rho for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Ff _bctype) {
            BoundaryConditionSetRhoAlongXEdge(_p, 0, -1, _rhobc, _usbc, _bctype);          //  On xmin
            BoundaryConditionSetRhoAlongXEdge(_p, _p.lx - 1, 1, _rhobc, _usbc, _bctype);   //  On xmax
            BoundaryConditionSetRhoAlongYEdge(_p, 0, -1, _rhobc, _usbc, _bctype);          //  On ymin
            BoundaryConditionSetRhoAlongYEdge(_p, _p.ly - 1, 1, _rhobc, _usbc, _bctype);   //  On ymax
        }
    
        //  Function of setting boundary condition of NS set rho for 3D
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void BoundaryConditionSetRho(P<T>& _p, Fv0 _rhobc, Fv1 _usbc, Fv2 _utbc, Ff _bctype) {
            BoundaryConditionSetRhoAlongXFace(_p, 0, -1, _rhobc, _usbc, _utbc, _bctype);          //  On xmin
            BoundaryConditionSetRhoAlongXFace(_p, _p.lx - 1, 1, _rhobc, _usbc, _utbc, _bctype);   //  On xmax
            BoundaryConditionSetRhoAlongYFace(_p, 0, -1, _rhobc, _usbc, _utbc, _bctype);          //  On ymin
            BoundaryConditionSetRhoAlongYFace(_p, _p.ly - 1, 1, _rhobc, _usbc, _utbc, _bctype);   //  On ymax
            BoundaryConditionSetRhoAlongZFace(_p, 0, -1, _rhobc, _usbc, _utbc, _bctype);          //  On zmin
            BoundaryConditionSetRhoAlongZFace(_p, _p.lz - 1, 1, _rhobc, _usbc, _utbc, _bctype);   //  On zmax
        }
    }
}