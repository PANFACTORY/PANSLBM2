//*****************************************************************************
//  Title       :   src/equation/adjointnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/adjointnavierstokes_avx.h"
#endif

namespace PANSLBM2 {
    namespace ANS {
        //  Function of updating macroscopic values of ANS for 2D
        template<class T, template<class>class P>
        void Macro(T &_ip, T &_iux, T &_iuy, T &_imx, T &_imy, T _rho, T _ux, T _uy, const T *_f0, const T *_f, int _idx) {
            T uu = _ux*_ux + _uy*_uy;
            _ip = _f0[_idx]*P<T>::ei[0]*(1.0 - 1.5*uu);
            _iux = -_f0[_idx]*P<T>::ei[0]*_ux;
            _iuy = -_f0[_idx]*P<T>::ei[0]*_uy;
            _imx = T();
            _imy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy;
                T fei = _f[P<T>::IndexF(_idx, c)]*P<T>::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P<T>::cx[c] + 3.0*ciu*P<T>::cx[c] - _ux);
                _iuy += fei*(P<T>::cy[c] + 3.0*ciu*P<T>::cy[c] - _uy);
                _imx += fei*P<T>::cx[c];
                _imy += fei*P<T>::cy[c];
            }
        }

        //  Function of updating macroscopic values of ANS for 3D
        template<class T, template<class>class P>
        void Macro(T &_ip, T &_iux, T &_iuy, T &_iuz, T &_imx, T &_imy, T &_imz, T _rho, T _ux, T _uy, T _uz, const T *_f0, const T *_f, int _idx) {
            T uu = _ux*_ux + _uy*_uy + _uz*_uz;
            _ip = _f0[_idx]*P<T>::ei[0]*(1.0 - 1.5*uu);
            _iux = -_f0[_idx]*P<T>::ei[0]*_ux;
            _iuy = -_f0[_idx]*P<T>::ei[0]*_uy;
            _iuz = -_f0[_idx]*P<T>::ei[0]*_uz;
            _imx = T();
            _imy = T();
            _imz = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T ciu = P<T>::cx[c]*_ux + P<T>::cy[c]*_uy + P<T>::cz[c]*_uz;
                T fei = _f[P<T>::IndexF(_idx, c)]*P<T>::ei[c];
                _ip += fei*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                _iux += fei*(P<T>::cx[c] + 3.0*ciu*P<T>::cx[c] - _ux);
                _iuy += fei*(P<T>::cy[c] + 3.0*ciu*P<T>::cy[c] - _uy);
                _iuz += fei*(P<T>::cz[c] + 3.0*ciu*P<T>::cz[c] - _uz);
                _imx += fei*P<T>::cx[c];
                _imy += fei*P<T>::cy[c];
                _imz += fei*P<T>::cz[c];
            }
        }

        //  Function of getting equilibrium of ANS for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _ux, T _uy, T _ip, T _iux, T _iuy) {
            for (int c = 0; c < P<T>::nc; ++c) {
                _feq[c] = _ip + 3.0*(_iux*(P<T>::cx[c] - _ux) + _iuy*(P<T>::cy[c] - _uy));
            }
        }

        //  Function of getting equilibrium of ANS for 3D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _ux, T _uy, T _uz, T _ip, T _iux, T _iuy, T _iuz) {
            for (int c = 0; c < P<T>::nc; ++c) {
                _feq[c] = _ip + 3.0*(_iux*(P<T>::cx[c] - _ux) + _iuy*(P<T>::cy[c] - _uy) + _iuz*(P<T>::cz[c] - _uz));
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _imx, T _imy, T *_f0, T *_f, T _alpha, int _idx) {
            T coef = 3.0*_alpha/(_rho + _alpha);
            _f0[_idx] -= -coef*(_ux*_imx + _uy*_imy);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= coef*((P<T>::cx[c] - _ux)*_imx + (P<T>::cy[c] - _uy)*_imy);
            }
        }

        //  Function of applying external force with Brinkman model of ANS for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _uz, T _imx, T _imy, T _imz, T *_f0, T *_f, T _alpha, int _idx) {
            T coef = 3.0*_alpha/(_rho + _alpha);
            _f0[_idx] -= -coef*(_ux*_imx + _uy*_imy + _uz*_imz);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] -= coef*((P<T>::cx[c] - _ux)*_imx + (P<T>::cy[c] - _uy)*_imy + (P<T>::cz[c] - _uz)*_imz);
            }
        }

        //  Function of setting boundary condition of ANS set iU for D2Q9 along x edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetUAlongXEdge(P<T>& _p, int _i, int _directionx, Fv0 _uxbc, Fv1 _uybc, Ff _bctype, T _eps = T()) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T ux = _uxbc(i + _p.offsetx, j + _p.offsety), uy = _uybc(i + _p.offsetx, j + _p.offsety);
                        if (_directionx == -1) {
                            T rho0 = (-2.0*_eps + ux*(4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)]) + 3.0*uy*(_p.f[P<T>::IndexF(idx, 5)] - _p.f[P<T>::IndexF(idx, 8)]))/(3.0*(1.0 - ux));
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] + rho0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                        } else if (_directionx == 1) {
                            T rho0 = (-2.0*_eps - ux*(4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)]) + 3.0*uy*(_p.f[P<T>::IndexF(idx, 6)] - _p.f[P<T>::IndexF(idx, 7)]))/(3.0*(1.0 + ux));
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + rho0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D2Q9 along y edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetUAlongYEdge(P<T>& _p, int _j, int _directiony, Fv0 _uxbc, Fv1 _uybc, Ff _bctype, T _eps = T()) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T ux = _uxbc(i + _p.offsetx, j + _p.offsety), uy = _uybc(i + _p.offsetx, j + _p.offsety);
                        if (_directiony == -1) {
                            T rho0 = (-2.0*_eps + ux*(4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)]) + 3.0*ux*(_p.f[P<T>::IndexF(idx, 5)] - _p.f[P<T>::IndexF(idx, 6)]))/(3.0*(1.0 - ux));
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] + rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0;
                        } else if (_directiony == 1) {
                            T rho0 = (-2.0*_eps - ux*(4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)]) + 3.0*ux*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 7)]))/(3.0*(1.0 + ux));
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + rho0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                        }
                    }
                }
            }
        }
        
        //  Function of setting boundary condition of ANS set iU for D3Q15 along x face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void iBoundaryConditionSetUAlongXFace(P<T>& _p, int _i, int _directionx, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype, T _eps = T()) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionx == -1) {
                                T rho0 = (-4.0*_eps + ux*(8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])
                                    + 3.0*uy*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 12)])
                                    + 3.0*uz*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 12)])
                                )/(6.0*(1.0 - ux));
                                _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] + rho0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0;
                            } else if (_directionx == 1) {
                                T rho0 = (-4.0*_eps - ux*(8.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])
                                    + 3.0*uy*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] - _p.f[P<T>::IndexF(idx, 14)])
                                    + 3.0*uz*(_p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])
                                )/(6.0*(1.0 + ux));
                                _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + rho0;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D3Q15 along y face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void iBoundaryConditionSetUAlongYFace(P<T>& _p, int _j, int _directiony, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype, T _eps = T()) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directiony == -1) {
                                T rho0 = (-4.0*_eps + 3.0*ux*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 13)]) 
                                    + uy*(8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])
                                    + 3.0*uz*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 13)])
                                )/(6.0*(1.0 - uy));
                                _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] + rho0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0;
                            } else if (_directiony == 1) {
                                T rho0 = (-4.0*_eps + 3.0*ux*(_p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] - _p.f[P<T>::IndexF(idx, 14)]) 
                                    - uy*(8.0*_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])
                                    + 3.0*uz*(_p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])
                                )/(6.0*(1.0 + uy));
                                _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + rho0;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D3Q15 along z face
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void iBoundaryConditionSetUAlongZFace(P<T>& _p, int _k, int _directionz, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype, T _eps = T()) {
            int k = _k - _p.offsetz;
            if (0 <= k && k < _p.nz) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            T ux = _uxbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uy = _uybc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz), uz = _uzbc(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                            if (_directionz == -1) {
                                T rho0 = (-4.0*_eps + 3.0*ux*(_p.f[P<T>::IndexF(idx, 7)] - _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 14)])
                                    + 3.0*uy*(_p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] - _p.f[P<T>::IndexF(idx, 9)] - _p.f[P<T>::IndexF(idx, 14)])
                                    + uz*(8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])
                                )/(6.0*(1.0 - uz));
                                _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] + rho0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0;
                            } else if (_directionz == 1) {
                                T rho0 = (-4.0*_eps + 3.0*ux*(_p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] - _p.f[P<T>::IndexF(idx, 13)])
                                    + 3.0*uy*(_p.f[P<T>::IndexF(idx, 10)] - _p.f[P<T>::IndexF(idx, 11)] - _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])
                                    - uz*(8.0*_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])
                                )/(6.0*(1.0 + uz));
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
        }

        //  Function of setting boundary condition of ANS set iRho for D2Q9 along x edge
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRhoAlongXEdge(P<T>& _p, int _i, int _directionx, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        if (_directionx == -1) {
                            T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - rho0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - rho0;
                        } else if (_directionx == 1) {
                            T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0;
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - rho0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - rho0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iRho for D2Q9 along y edge
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRhoAlongYEdge(P<T>& _p, int _j, int _directiony, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        if (_directiony == -1) {
                            T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0;
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - rho0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] - rho0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] - rho0;
                        } else if (_directiony == 1) {
                            T rho0 = (4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - rho0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of ANS set iRho for D3Q15 along x face
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRhoAlongXFace(P<T>& _p, int _i, int _directionx, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            if (_directionx == -1) {
                                T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])/6.0;
                                _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] - rho0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] - rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - rho0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - rho0;
                            } else if (_directionx == 1) {
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
            }
        }

        //  Function of setting boundary condition of ANS set iRho for D3Q15 along y face
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRhoAlongYFace(P<T>& _p, int _j, int _directiony, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            if (_directiony == -1) {
                                T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                                _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] - rho0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] - rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] - rho0;
                            } else if (_directiony == 1) {
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
            }
        }

        //  Function of setting boundary condition of ANS set iRho for D3Q15 along z face
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRhoAlongZFace(P<T>& _p, int _k, int _directionz, Ff _bctype) {
            int k = _k - _p.offsetz;
            if (0 <= k && k < _p.nz) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, j, k);
                            if (_directionz == -1) {
                                T rho0 = (8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                                _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] - rho0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] - rho0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] - rho0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] - rho0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] - rho0;
                            } else if (_directionz == 1) {
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

    namespace ANS {
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
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _p.f0, _p.f, _alpha[idx], idx);
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
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman model) and Collide of ANS for 3D
        template<class T, template<class>class P>
        void MacroBrinkmanCollide(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, 
            T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, 
            T _viscosity, const T *_alpha, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _p.f0, _p.f, _alpha[idx], idx);
                Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

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
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        } 
    
        //  Function of Update macro and collide of ANS with LSM for 2D
        template<class T, template<class>class P>
        void MacroCollideLSM(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, 
            T _viscosity, const T *_chi, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx]*_chi[idx], _uy[idx]*_chi[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                }

                //  Collide
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], ip, iux*_chi[idx], iuy*_chi[idx]);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of Update macro and collide of ANS with LSM for 3D
        template<class T, template<class>class P>
        void MacroCollideLSM(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, 
            T _viscosity, const T *_chi, bool _issave = false
        ) {
            T omega = 1.0/(3.0*_viscosity + 0.5), iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx]*_chi[idx], _uy[idx]*_chi[idx], _uz[idx]*_chi[idx], _p.f0, _p.f, idx);

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

                //  Collide
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux*_chi[idx], iuy*_chi[idx], iuz*_chi[idx]);
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
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _ip[idx], _iux[idx], _iuy[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting initial condition of ANS for 3D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_ux, const T *_uy, const T *_uz, const T *_ip, const T *_iux, const T *_iuy, const T *_iuz) {
            T feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], _ip[idx], _iux[idx], _iuy[idx], _iuz[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting boundary condition of ANS set iU for D2Q9
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Ff _bctype, T _eps = T()) {
            iBoundaryConditionSetUAlongXEdge(_p, 0, -1, _uxbc, _uybc, _bctype, _eps);           //  On xmin
            iBoundaryConditionSetUAlongXEdge(_p, _p.lx - 1, 1, _uxbc, _uybc, _bctype, _eps);    //  On xmax
            iBoundaryConditionSetUAlongYEdge(_p, 0, -1, _uxbc, _uybc, _bctype, _eps);           //  On ymin
            iBoundaryConditionSetUAlongYEdge(_p, _p.ly - 1, 1, _uxbc, _uybc, _bctype, _eps);    //  On ymax
        }
    
        //  Function of setting boundary condition of ANS set iU for D3Q15
        template<class T, template<class>class P, class Fv0, class Fv1, class Fv2, class Ff>
        void iBoundaryConditionSetU(P<T>& _p, Fv0 _uxbc, Fv1 _uybc, Fv2 _uzbc, Ff _bctype, T _eps = T()) {
            iBoundaryConditionSetUAlongXFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype, _eps);        //  On xmin
            iBoundaryConditionSetUAlongXFace(_p, _p.lx - 1, 1, _uxbc, _uybc, _uzbc, _bctype, _eps); //  On xmax
            iBoundaryConditionSetUAlongYFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype, _eps);        //  On ymin
            iBoundaryConditionSetUAlongYFace(_p, _p.ly - 1, 1, _uxbc, _uybc, _uzbc, _bctype, _eps); //  On ymax
            iBoundaryConditionSetUAlongZFace(_p, 0, -1, _uxbc, _uybc, _uzbc, _bctype, _eps);        //  On zmin
            iBoundaryConditionSetUAlongZFace(_p, _p.lz - 1, 1, _uxbc, _uybc, _uzbc, _bctype, _eps); //  On zmax
        }

        //  Function of setting boundary condition of ANS set iRho for D2Q9
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRho2D(P<T>& _p, Ff _bctype) {
            iBoundaryConditionSetRhoAlongXEdge(_p, 0, -1, _bctype);         //  On xmin
            iBoundaryConditionSetRhoAlongXEdge(_p, _p.lx - 1, 1, _bctype);  //  On xmax
            iBoundaryConditionSetRhoAlongYEdge(_p, 0, -1, _bctype);         //  On ymin
            iBoundaryConditionSetRhoAlongYEdge(_p, _p.ly - 1, 1, _bctype);  //  On ymax
        }

        //  Function of setting boundary condition of ANS set iRho for D3Q15
        template<class T, template<class>class P, class Ff>
        void iBoundaryConditionSetRho3D(P<T>& _p, Ff _bctype) {
            iBoundaryConditionSetRhoAlongXFace(_p, 0, -1, _bctype);         //  On xmin
            iBoundaryConditionSetRhoAlongXFace(_p, _p.lx - 1, 1, _bctype);  //  On xmax
            iBoundaryConditionSetRhoAlongYFace(_p, 0, -1, _bctype);         //  On ymin
            iBoundaryConditionSetRhoAlongYFace(_p, _p.ly - 1, 1, _bctype);  //  On ymax
            iBoundaryConditionSetRhoAlongZFace(_p, 0, -1, _bctype);         //  On zmin
            iBoundaryConditionSetRhoAlongZFace(_p, _p.lz - 1, 1, _bctype);  //  On zmax
        }

        //  Function of getting sensitivity of Brinkman model
        template<class T, template<class>class P>
        void SensitivityBrinkman(P<T>& _p, T *_dfds, const T *_ux, const T *_uy, const T *_imx, const T *_imy, const T *_dads) {
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
            }
        }

        //  Function of getting sensitivity of Brinkman model
        template<class T, template<class>class P>
        void SensitivityBrinkman(P<T>& _p, T *_dfds, const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz, const T *_dads) {
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
            }
        }
    
        //  Function of getting sensitivity of Level-Set-Method
        template<class T, template<class>class P>
        void SensitivityLSM(P<T>& _p, T *_dfds, const T *_rho, const T *_ux, const T *_uy, const T *_iux, const T *_iuy, T _viscosity) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _dfds[idx] -= 3.0*omega*_rho[idx]*(_ux[idx]*_iux[idx] + _uy[idx]*_iuy[idx]);
            }
        }

        //  Function of getting sensitivity of Level-Set-Method
        template<class T, template<class>class P>
        void SensitivityLSM(P<T>& _p, T *_dfds, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_iux, const T *_iuy, const T *_iuz, T _viscosity) {
            T omega = 1.0/(3.0*_viscosity + 0.5);
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                _dfds[idx] -= 3.0*omega*_rho[idx]*(_ux[idx]*_iux[idx] + _uy[idx]*_iuy[idx] + _uz[idx]*_iuz[idx]);
            }
        }
    }
}