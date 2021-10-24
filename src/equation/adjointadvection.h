//*****************************************************************************
//  Title       :   src/equation/adjointadvection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include "adjointnavierstokes.h"
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/adjointadvection_avx.h"
#endif

namespace {
    const int SetT = 1;
    const int SetQ = 2;
}

namespace PANSLBM2 {
    namespace AAD {
        //  Function of updating macroscopic values of AAD for 2D
        template<class T, template<class>class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, const T *_g0, const T *_g, int _idx) {
            _item = Q<T>::ei[0]*_g0[_idx];
            _iqx = T();
            _iqy = T();
            for (int c = 1; c <Q<T>::nc; ++c) {
                T gei = Q<T>::ei[c]*_g[Q<T>::IndexF(_idx, c)];
                _item += gei;
                _iqx += gei*Q<T>::cx[c];
                _iqy += gei*Q<T>::cy[c];
            }
        }

        //  Function of updating macroscopic values of AAD for 3D
        template<class T, template<class>class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, T &_iqz, const T *_g0, const T *_g, int _idx) {
            _item = Q<T>::ei[0]*_g0[_idx];
            _iqx = T();
            _iqy = T();
            _iqz = T();
            for (int c = 1; c <Q<T>::nc; ++c) {
                T gei = Q<T>::ei[c]*_g[Q<T>::IndexF(_idx, c)];
                _item += gei;
                _iqx += gei*Q<T>::cx[c];
                _iqy += gei*Q<T>::cy[c];
                _iqz += gei*Q<T>::cz[c];
            }
        }

        //  Function of getting equilibrium of AAD for 2D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _item, T _iqx, T _iqy, T _ux, T _uy) {
            T coef = _item + 3.0*(_ux*_iqx + _uy*_iqy);
            for (int c = 0; c < Q<T>::nc; ++c) {
                _geq[c] = coef;
            }
        }

        //  Function of getting equilibrium of AAD for 3D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _item, T _iqx, T _iqy, T _iqz, T _ux, T _uy, T _uz) {
            T coef = _item + 3.0*(_ux*_iqx + _uy*_iqy + _uz*_iqz);
            for (int c = 0; c < Q<T>::nc; ++c) {
                _geq[c] = coef;
            }
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _imx, T _imy, T _tem, T _iqx, T _iqy, T _omegag, T *_f0, T *_f, T _alpha, int _idx) {
            T coef = 3.0/(_rho + _alpha);
            _f0[_idx] += coef*(
                -_ux*(_tem*_iqx*_omegag - _alpha*_imx) + 
                -_uy*(_tem*_iqy*_omegag - _alpha*_imy)
            );
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += coef*(
                    (P<T>::cx[c] - _ux)*(_tem*_iqx*_omegag - _alpha*_imx) + 
                    (P<T>::cy[c] - _uy)*(_tem*_iqy*_omegag - _alpha*_imy)
                );
            }
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(T _rho, T _ux, T _uy, T _uz, T _imx, T _imy, T _imz, T _tem, T _iqx, T _iqy, T _iqz, T _omegag, T *_f0, T *_f, T _alpha, int _idx) {
            T coef = 3.0/(_rho + _alpha);
            _f0[_idx] += coef*(
                -_ux*(_tem*_iqx*_omegag - _alpha*_imx) + 
                -_uy*(_tem*_iqy*_omegag - _alpha*_imy) +
                -_uz*(_tem*_iqz*_omegag - _alpha*_imz)
            );
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += coef*(
                    (P<T>::cx[c] - _ux)*(_tem*_iqx*_omegag - _alpha*_imx) + 
                    (P<T>::cy[c] - _uy)*(_tem*_iqy*_omegag - _alpha*_imy) +
                    (P<T>::cz[c] - _uz)*(_tem*_iqz*_omegag - _alpha*_imz)
                );
            }
        }

        //  Function of applying external force with heat exchange of AAD for 2D/3D
        template<class T, template<class>class Q>
        void ExternalForceHeatExchange(T _item, T *_g0, T *_g, T _beta, int _idx) {
            T coef = _beta*(1.0 + _item)/(1.0 + _beta);
            _g0[_idx] -= coef;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] -= coef;
            }
        }

        //  Function of applying external force with natural convection of AAD for 2D
        template<class T, template<class>class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _gx, T _gy, T *_g0, T *_g, int _idx) {
            T coef = 3.0*(_imx*_gx + _imy*_gy);
            _g0[_idx] += coef;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += coef;
            }
        }

        //  Function of applying external force with natural convection of AAD for 3D
        template<class T, template<class>class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _imz, T _gx, T _gy, T _gz, T *_g0, T *_g, int _idx) {
            T coef = 3.0*(_imx*_gx + _imy*_gy + _imz*_gz);
            _g0[_idx] += coef;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += coef;
            }
        }

        //  Function of setting boundary condition set iT of AAD for D2Q9 along x edge
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetTAlongXEdge(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        if (_directionx == -1) {
                            T rho0 = -(4.0*(1.0 + 3.0*_ux[idx])*_q.f[Q<T>::IndexF(idx, 1)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 5)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 8)])/(6.0*(1.0 + 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                        } else if (_directionx == 1) {
                            T rho0 = -(4.0*(1.0 - 3.0*_ux[idx])*_q.f[Q<T>::IndexF(idx, 3)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 6)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 7)])/(6.0*(1.0 - 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D2Q9 along y edge
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetTAlongYEdge(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        if (_directiony == -1) {
                            T rho0 = -(4.0*(1.0 + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 2)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 5)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 6)])/(6.0*(1.0 + 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                        } else if (_directiony == 1) {
                            T rho0 = -(4.0*(1.0 - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 4)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 7)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 8)])/(6.0*(1.0 - 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                        }
                    }
                }
            }
        }    
        
        //  Function of setting boundary condition set iT of AAD for D3Q15 along x face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetTAlongXFace(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directionx == -1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/12.0
                                    -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]))
                                    -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]));
                                _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            } else if (_directionx == 1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                    -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]))
                                    -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]));
                                _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D3Q15 along y face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetTAlongYFace(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directiony == -1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/12.0
                                    -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]))
                                    -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]));
                                _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            } else if (_directiony == 1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                    -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]))
                                    -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]));
                                _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D3Q15 along z face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetTAlongZFace(Q<T>& _q, int _k, int _directionz, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directionz == -1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                    -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]))
                                    -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]));
                                _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            } else if (_directionz == 1) {
                                T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/12.0
                                    -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]))
                                    -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]));
                                _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iQ of AAD for D2Q9 along x edge
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQAlongXEdge(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, Ff _bctype, T _eps = T()) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        if (_directionx == -1) {
                            T rho0 = (
                                (1.0 + 3.0*_ux[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])
                                + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])
                                - 12.0*_eps
                            )/(6.0*(1.0 - 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                        } else if (_directionx == 1) {
                            T rho0 = (
                                (1.0 - 3.0*_ux[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])
                                + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])
                                - 12.0*_eps
                            )/(6.0*(1.0 + 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                        }
                    }
                }
            }
        }
            
        //  Function of setting boundary condition set iQ of AAD for D2Q9 along y edge
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQAlongYEdge(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, Ff _bctype, T _eps = T()) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        if (_directiony == -1) {
                            T rho0 = (
                                (1.0 + 3.0*_uy[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])
                                + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])
                                - 12.0*_eps
                            )/(6.0*(1.0 - 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                        } else if (_directiony == 1) {
                            T rho0 = (
                                (1.0 - 3.0*_uy[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])
                                + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])
                                - 12.0*_eps
                            )/(6.0*(1.0 + 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                        }
                    }
                }
            }
        }    
        
        //  Function of setting boundary condition set iQ of AAD for D3Q15 along x face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQAlongXFace(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, const T *_uz, Ff _bctype, T _eps = T()) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directionx == -1) {
                                T rho0 = (
                                    (1.0 + 3.0*_ux[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])
                                    + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])
                                    + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 - 3.0*_ux[idx]));
                                _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            } else if (_directionx == 1) {
                                T rho0 = (
                                    (1.0 - 3.0*_ux[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])
                                    + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])
                                    + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 + 3.0*_ux[idx]));
                                _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iQ of AAD for D3Q15 along y face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQAlongYFace(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, const T *_uz, Ff _bctype, T _eps = T()) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directiony == -1) {
                                T rho0 = (
                                    (1.0 + 3.0*_uy[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])
                                    + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])
                                    + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 - 3.0*_uy[idx]));
                                _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            } else if (_directiony == 1) {
                                T rho0 = (
                                    (1.0 - 3.0*_uy[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 + 3.0*_uy[idx]));
                                _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iQ of AAD for D3Q15 along z face
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQAlongZFace(Q<T>& _q, int _k, int _directionz, const T *_ux, const T *_uy, const T *_uz, Ff _bctype, T _eps = T()) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            if (_directionz == -1) {
                                T rho0 = (
                                    (1.0 + 3.0*_uz[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 - 3.0*_uz[idx]));
                                _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            } else if (_directionz == 1) {
                                T rho0 = (
                                    (1.0 - 3.0*_uz[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])
                                    - 24.0*_eps
                                )/(12.0*(1.0 + 3.0*_uz[idx]));
                                _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                                _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D2Q9 along x edge
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRhoAlongXEdge(P<T>& _p, Q<T>& _q, int _i, int _directionx, const T *_rho, const T *_ux, const T *_uy, const T *_tem, Ff _bctype, T _eps = T()) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    int bctype = _bctype(i + _p.offsetx, j + _p.offsety);
                    if (bctype) {
                        int idx = _q.Index(i, j);
                        if (_directionx == -1) {
                            T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                            T flux0 = T();
                            if (bctype == SetT) {
                                flux0 = _tem[idx]*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/(2.0*(1.0 + 3.0*_ux[idx])*_rho[idx]);
                            } else if (bctype == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (4.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/3.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/2.0
                                )/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                        } else if (_directionx == 1) {
                            T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0;
                            T flux0 = T();
                            if (bctype == SetT) {
                                flux0 = _tem[idx]*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_ux[idx])*_rho[idx]);
                            } else if (bctype == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (4.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/3.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/2.0
                                )/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D2Q9 along y edge
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRhoAlongYEdge(P<T>& _p, Q<T>& _q, int _j, int _directiony, const T *_rho, const T *_ux, const T *_uy, const T *_tem, Ff _bctype, T _eps = T()) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    int bctype = _bctype(i + _p.offsetx, j + _p.offsety);
                    if (bctype) { 
                        int idx = _q.Index(i, j);
                        if (_directiony == -1) {
                            T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0;
                            T flux0 = T();
                            if (bctype == SetT) {
                                flux0 = _tem[idx]*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/(2.0*(1.0 + 3.0*_uy[idx])*_rho[idx]);
                            } else if (bctype == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (4.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/3.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/2.0
                                )/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                        } else if (_directiony == 1) {
                            T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                            T flux0 = T();
                            if (bctype == SetT) {
                                flux0 = _tem[idx]*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_uy[idx])*_rho[idx]);
                            } else if (bctype == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (4.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/3.0 
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])/2.0
                                )/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
        }
        
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15 along x face
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRhoAlongXFace(P<T>& _p, Q<T>& _q, int _i, int _directionx, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, Ff _bctypeg, T _eps = T()) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int bctype = _bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                        if (bctype) {
                            int idx = _p.Index(i, j, k);
                            if (_directionx == -1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/6.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                            } else if (_directionx == 1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                            }
                        }
                    }
                }
            }
        }
            
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15 along y face
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRhoAlongYFace(P<T>& _p, Q<T>& _q, int _j, int _directiony, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, Ff _bctypeg, T _eps = T()) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        int bctype = _bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                        if (bctype) {
                            int idx = _p.Index(i, j, k);
                            if (_directiony == -1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/6.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                            } else if (_directiony == 1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                        + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15 along z face
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRhoAlongZFace(P<T>& _p, Q<T>& _q, int _k, int _directionz, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, Ff _bctypeg, T _eps = T()) {
            int k = _k - _p.offsetz;
            if (0 <= k && k < _p.nz) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        int bctype = _bctype(i + _p.offsetx, j + _p.offsety, k + _p.offsetz);
                        if (bctype) {
                            int idx = _p.Index(i, j, k);
                            if (_directionz == -1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uz[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                            } else if (_directionz == 1) {
                                T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                                T flux0 = T();
                                if (bctype == SetT) {
                                    flux0 = _tem[idx]*(
                                        _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                                } else if (bctype == SetQ) {
                                    flux0 = -_tem[idx]*(
                                        (8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/6.0
                                        + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                        + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                                }
                                T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uz[idx])*_rho[idx]);
                                _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                                _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                            }
                        }
                    }
                }
            }
        }  
    
        //  Function of getting sensitivity of temperature at heat source for D2Q9 along x edge
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongXEdge(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, const T *_ig, T *_dfds, const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j), offsetf = Q<T>::nc*idx;
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directionx == -1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 + 3.0*_ux[idx])*(-6.0 + 4.0*_ig[offsetf + 1] + _ig[offsetf + 5] + _ig[offsetf + 8])
                                + 3.0*_uy[idx]*(_ig[offsetf + 5] - _ig[offsetf + 8])
                            )/(36.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        } else if (_directionx == 1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 - 3.0*_ux[idx])*(-6.0 + 4.0*_ig[offsetf + 3] + _ig[offsetf + 6] + _ig[offsetf + 7])
                                + 3.0*_uy[idx]*(_ig[offsetf + 6] - _ig[offsetf + 7])
                            )/(36.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D2Q9 along y edge
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongYEdge(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, const T *_ig, T *_dfds, const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j), offsetf = Q<T>::nc*idx;
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directiony == -1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 + 3.0*_uy[idx])*(-6.0 + 4.0*_ig[offsetf + 2] + _ig[offsetf + 5] + _ig[offsetf + 6])
                                + 3.0*_ux[idx]*(_ig[offsetf + 5] - _ig[offsetf + 6])
                            )/(36.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        } else if (_directiony == 1) {
                            _dfds[idx] += qn*_dkds[idx]*(
                                (1.0 - 3.0*_uy[idx])*(-6.0 + 4.0*_ig[offsetf + 4] + _ig[offsetf + 7] + _ig[offsetf + 8])
                                + 3.0*_ux[idx]*(_ig[offsetf + 8] - _ig[offsetf + 7])
                            )/(36.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
        }
        
        //  Function of getting sensitivity of temperature at heat source for D3Q15 along x face
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongXFace(Q<T>& _q, int _i, int _directionx, const T *_ux, const T *_uy, const T *_ig, const T *_uz, T *_dfds, const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k), offsetf = Q<T>::nc*idx;
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionx == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_ux[idx])*(-12.0 + 8.0*_ig[offsetf + 1] + _ig[offsetf + 7] + _ig[offsetf + 9] + _ig[offsetf + 10] + _ig[offsetf + 12])
                                    + 3.0*_uy[idx]*(_ig[offsetf + 7] - _ig[offsetf + 9] + _ig[offsetf + 10] - _ig[offsetf + 12])
                                    + 3.0*_uz[idx]*(_ig[offsetf + 7] + _ig[offsetf + 9] - _ig[offsetf + 10] - _ig[offsetf + 12])
                                )/(72.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directionx == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_ux[idx])*(-12.0 + 8.0*_ig[offsetf + 4] + _ig[offsetf + 8] + _ig[offsetf + 11] + _ig[offsetf + 13] + _ig[offsetf + 14])
                                    + 3.0*_uy[idx]*(_ig[offsetf + 8] - _ig[offsetf + 11] + _ig[offsetf + 13] - _ig[offsetf + 14])
                                    + 3.0*_uz[idx]*(_ig[offsetf + 8] - _ig[offsetf + 11] - _ig[offsetf + 13] + _ig[offsetf + 14])
                                )/(72.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15 along y face
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongYFace(Q<T>& _q, int _j, int _directiony, const T *_ux, const T *_uy, const T *_ig, const T *_uz, T *_dfds, const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k), offsetf = Q<T>::nc*idx;
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directiony == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_uy[idx])*(-12.0 + 8.0*_ig[offsetf + 2] + _ig[offsetf + 7] + _ig[offsetf + 8] + _ig[offsetf + 10] + _ig[offsetf + 13])
                                    + 3.0*_uz[idx]*(_ig[offsetf + 7] + _ig[offsetf + 8] - _ig[offsetf + 10] - _ig[offsetf + 13])
                                    + 3.0*_ux[idx]*(_ig[offsetf + 7] - _ig[offsetf + 8] + _ig[offsetf + 10] - _ig[offsetf + 13])
                                )/(72.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directiony == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_uy[idx])*(-12.0 + 8.0*_ig[offsetf + 5] + _ig[offsetf + 9] + _ig[offsetf + 11] + _ig[offsetf + 12] + _ig[offsetf + 14])
                                    + _uz[idx]*(_ig[offsetf + 9] - _ig[offsetf + 11] - _ig[offsetf + 12] + _ig[offsetf + 14])
                                    + _ux[idx]*(_ig[offsetf + 9] - _ig[offsetf + 11] + _ig[offsetf + 12] - _ig[offsetf + 14])
                                )/(72.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }
            
        //  Function of getting sensitivity of temperature at heat source for D3Q15 along z face
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSourceAlongZFace(Q<T>& _q, int _k, int _directionz, const T *_ux, const T *_uy, const T *_ig, const T *_uz, T *_dfds, const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k), offsetf = Q<T>::nc*idx;
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionz == -1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 + 3.0*_uz[idx])*(-12.0 + 8.0*_ig[offsetf + 3] + _ig[offsetf + 7] + _ig[offsetf + 8] + _ig[offsetf + 9] + _ig[offsetf + 14])
                                    + _ux[idx]*(_ig[offsetf + 7] - _ig[offsetf + 8] + _ig[offsetf + 9] - _ig[offsetf + 14])
                                    + _uy[idx]*(_ig[offsetf + 7] + _ig[offsetf + 8] - _ig[offsetf + 9] - _ig[offsetf + 14])
                                )/(72.0*(1.0 - 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                            } else if (_directionz == 1) {
                                _dfds[idx] += qn*_dkds[idx]*(
                                    (1.0 - 3.0*_uz[idx])*(-12.0 + 8.0*_ig[offsetf + 6] + _ig[offsetf + 10] + _ig[offsetf + 11] + _ig[offsetf + 12] + _ig[offsetf + 13])
                                    + _ux[idx]*(_ig[offsetf + 10] - _ig[offsetf + 11] + _ig[offsetf + 12] - _ig[offsetf + 13])
                                    + _uy[idx]*(_ig[offsetf + 10] - _ig[offsetf + 11] - _ig[offsetf + 12] + _ig[offsetf + 13])
                                )/(72.0*(1.0 + 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                            }
                        }
                    }
                }
            }
        }   
    }

    namespace AAD {
        //  Function of Update macro, External force(Brinkman, Heat exchange) and Collide of AAD for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc]; 
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(item, _q.f0, _q.f, _beta[idx], idx);
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

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

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c); 
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange) and Collide of AAD for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(item, _q.f0, _q.f, _beta[idx], idx);
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of AAD for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, bool _issave = false, T *_ig = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);

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

                    if (_ig) {
                        int offsetf = Q<T>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, bool _issave = false, T *_ig = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_ig) {
                        int offsetf = Q<T>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, T _gx, T _gy, bool _issave = false, T *_ig = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], imx, imy, _tem[idx], iqx, iqy, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, _gx, _gy, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

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

                    if (_ig) {
                        int offsetf = Q<T>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], ip, iux, iuy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, _ux[idx], _uy[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of AAD for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, T _gx, T _gy, T _gz, bool _issave = false, T *_ig = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho[idx], _ux[idx], _uy[idx], _uz[idx], imx, imy, imz, _tem[idx], iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha[idx], idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho[idx], _ux[idx], _uy[idx], _uz[idx], _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;

                    if (_ig) {
                        int offsetf = Q<T>::nc*idx;
                        _ig[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _ig[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                ANS::Equilibrium<T, P>(feq, _ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of setting initial condition of AAD for 2D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_ux, const T *_uy, const T *_item, const T *_iqx, const T *_iqy) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _item[idx], _iqx[idx], _iqy[idx], _ux[idx], _uy[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting initial condition of AAD for 3D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _item[idx], _iqx[idx], _iqy[idx], _iqz[idx], _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D2Q9
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetT(Q<T>& _q, const T *_ux, const T *_uy, Ff _bctype) {
            iBoundaryConditionSetTAlongXEdge(_q, 0, -1, _ux, _uy, _bctype);         //  On xmin
            iBoundaryConditionSetTAlongXEdge(_q, _q.lx - 1, 1, _ux, _uy, _bctype);  //  On xmax
            iBoundaryConditionSetTAlongYEdge(_q, 0, -1, _ux, _uy, _bctype);         //  On ymin
            iBoundaryConditionSetTAlongYEdge(_q, _q.ly - 1, 1, _ux, _uy, _bctype);  //  On ymax
        }

        //  Function of setting boundary condition set iT of AAD for D3Q15
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetT(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            iBoundaryConditionSetTAlongXFace(_q, 0, -1, _ux, _uy, _uz, _bctype);        //  On xmin
            iBoundaryConditionSetTAlongXFace(_q, _q.lx - 1, 1, _ux, _uy, _uz, _bctype); //  On xmax
            iBoundaryConditionSetTAlongYFace(_q, 0, -1, _ux, _uy, _uz, _bctype);        //  On ymin
            iBoundaryConditionSetTAlongYFace(_q, _q.ly - 1, 1, _ux, _uy, _uz, _bctype); //  On ymax
            iBoundaryConditionSetTAlongZFace(_q, 0, -1, _ux, _uy, _uz, _bctype);        //  On zmin
            iBoundaryConditionSetTAlongZFace(_q, _q.lz - 1, 1, _ux, _uy, _uz, _bctype); //  On zmax
        }
    
        //  Function of setting boundary condition set iQ of AAD for D2Q9
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQ(Q<T>& _q, const T *_ux, const T *_uy, Ff _bctype, T _eps = T()) {
            iBoundaryConditionSetQAlongXEdge(_q, 0, -1, _ux, _uy, _bctype, _eps);           //  On xmin
            iBoundaryConditionSetQAlongXEdge(_q, _q.lx - 1, 1, _ux, _uy, _bctype, _eps);    //  On xmax
            iBoundaryConditionSetQAlongYEdge(_q, 0, -1, _ux, _uy, _bctype, _eps);           //  On ymin
            iBoundaryConditionSetQAlongYEdge(_q, _q.ly - 1, 1, _ux, _uy, _bctype, _eps);    //  On ymax
        }

        //  Function of setting boundary condition set iQ of AAD for D3Q15
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQ(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, Ff _bctype, T _eps = T()) {
            iBoundaryConditionSetQAlongXFace(_q, 0, -1, _ux, _uy, _uz, _bctype, _eps);          //  On xmin
            iBoundaryConditionSetQAlongXFace(_q, _q.lx - 1, 1, _ux, _uy, _uz, _bctype, _eps);   //  On xmax
            iBoundaryConditionSetQAlongYFace(_q, 0, -1, _ux, _uy, _uz, _bctype, _eps);          //  On ymin
            iBoundaryConditionSetQAlongYFace(_q, _q.ly - 1, 1, _ux, _uy, _uz, _bctype, _eps);   //  On ymax
            iBoundaryConditionSetQAlongZFace(_q, 0, -1, _ux, _uy, _uz, _bctype, _eps);          //  On zmin
            iBoundaryConditionSetQAlongZFace(_q, _q.lz - 1, 1, _ux, _uy, _uz, _bctype, _eps);   //  On zmax
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D2Q9
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRho(P<T>& _p, Q<T>& _q, const T *_rho, const T *_ux, const T *_uy, const T *_tem, Ff _bctype, T _eps = T()) { 
            iBoundaryConditionSetRhoAlongXEdge(_p, _q, 0, -1, _rho, _ux, _uy, _tem, _bctype, _eps);         //  On xmin
            iBoundaryConditionSetRhoAlongXEdge(_p, _q, _p.lx - 1, 1, _rho, _ux, _uy, _tem, _bctype, _eps);  //  On xmax
            iBoundaryConditionSetRhoAlongYEdge(_p, _q, 0, -1, _rho, _ux, _uy, _tem, _bctype, _eps);         //  On ymin
            iBoundaryConditionSetRhoAlongYEdge(_p, _q, _p.ly - 1, 1, _rho, _ux, _uy, _tem, _bctype, _eps);  //  On ymax
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRho(P<T>& _p, Q<T>& _q, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, Ff _bctypeg, T _eps = T()) {
            iBoundaryConditionSetRhoAlongXFace(_p, _q, 0, -1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);           //  On xmin
            iBoundaryConditionSetRhoAlongXFace(_p, _q, _p.lx - 1, 1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);    //  On xmax
            iBoundaryConditionSetRhoAlongYFace(_p, _q, 0, -1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);           //  On ymin
            iBoundaryConditionSetRhoAlongYFace(_p, _q, _p.ly - 1, 1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);    //  On ymax
            iBoundaryConditionSetRhoAlongZFace(_p, _q, 0, -1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);           //  On zmin
            iBoundaryConditionSetRhoAlongZFace(_p, _q, _p.lz - 1, 1, _rho, _ux, _uy, _uz, _tem, _bctypeg, _eps);    //  On zmax
        }
    
        //  Function of getting sensitivity of heat exchange
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityHeatExchange(Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_imx, const T *_imy, const T *_dads, const T *_tem, const T *_item, const T *_dbds) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]) - _dbds[idx]*(1.0 - _tem[idx])*(1.0 + _item[idx]);
            }
        }

        //  Function of getting sensitivity of heat exchange
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityHeatExchange(Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz, const T *_dads, const T *_tem, const T *_item, const T *_dbds) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]) - _dbds[idx]*(1.0 - _tem[idx])*(1.0 + _item[idx]);
            }
        }

        //  Function of getting sensitivity of Brinkman model and diffusivity term
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityBrinkmanDiffusivity(
            Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_imx, const T *_imy, const T *_dads,
            const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_g, const T *_ig,
            const T *_diffusivity, const T *_dkds
        ) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                int offsetf = Q<T>::nc*idx;
                T sumg = T();
                for (int c = 0; c < Q<T>::nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
            }
        }

        //  Function of getting sensitivity of Brinkman model and diffusivity term
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityBrinkmanDiffusivity(
            Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz, const T *_dads,
            const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz, const T *_g, const T *_ig,
            const T *_diffusivity, const T *_dkds
        ) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                int offsetf = Q<T>::nc*idx;
                T sumg = T();
                for (int c = 0; c < Q<T>::nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D2Q9
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_imx, const T *_imy, const T *_dads,
            const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_g, const T *_ig,
            const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype
        ) {
            //  Brinkman term and diffusivity term
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                int offsetf = Q<T>::nc*idx;
                T sumg = T();
                for (int c = 0; c < Q<T>::nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
            }

            SensitivityTemperatureAtHeatSourceAlongXEdge(_q, 0, -1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);         //  Boundary term along xmin
            SensitivityTemperatureAtHeatSourceAlongXEdge(_q, _q.lx - 1, 1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);  //  Boundary term along xmax
            SensitivityTemperatureAtHeatSourceAlongYEdge(_q, 0, -1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);         //  Boundary term along ymin
            SensitivityTemperatureAtHeatSourceAlongYEdge(_q, _q.ly - 1, 1, _ux, _uy, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);  //  Boundary term along ymax
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            Q<T>& _q, T *_dfds, const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz, const T *_dads,
            const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz, const T *_g, const T *_ig,
            const T *_diffusivity, const T *_dkds, Fv _qnbc, Ff _bctype
        ) {
            //  Brinkman term and diffusivity term
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                int offsetf = Q<T>::nc*idx;
                T sumg = T();
                for (int c = 0; c < Q<T>::nc; ++c) {
                    sumg += _g[offsetf + c]*_ig[offsetf + c];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
            }

            SensitivityTemperatureAtHeatSourceAlongXFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along xmin
            SensitivityTemperatureAtHeatSourceAlongXFace(_q, _q.lx - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along xmax
            SensitivityTemperatureAtHeatSourceAlongYFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along ymin
            SensitivityTemperatureAtHeatSourceAlongYFace(_q, _q.ly - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along ymax
            SensitivityTemperatureAtHeatSourceAlongZFace(_q, 0, -1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype);        //  Boundary term along zmin
            SensitivityTemperatureAtHeatSourceAlongZFace(_q, _q.lz - 1, 1, _ux, _uy, _uz, _ig, _dfds, _diffusivity, _dkds, _qnbc, _bctype); //  Boundary term along zmax
        }
    }
}