//*****************************************************************************
//  Title       :   src/equation/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include "navierstokes.h"
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/advection_avx.h"
#endif

namespace PANSLBM2 {
    namespace AD {
        //  Function of updating macroscopic values of AD for 2D
        template<class T, template<class>class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T _ux, T _uy, const T *_g0, const T *_g, T _omegag, int _idx) {
            _tem = _g0[_idx];
            _qx = T();
            _qy = T();
            for (int c = 1; c < Q<T>::nc; ++c) {
                T g = _g[Q<T>::IndexF(_idx, c)];
                _tem += g;
                _qx += Q<T>::cx[c]*g;
                _qy += Q<T>::cy[c]*g;
            }
            T coef = 1.0 - 0.5*_omegag;
            _qx = coef*(_qx - _tem*_ux);
            _qy = coef*(_qy - _tem*_uy);
        }

        //  Function of updating macroscopic values of AD for 3D
        template<class T, template<class>class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T &_qz, T _ux, T _uy, T _uz, const T *_g0, const T *_g, T _omegag, int _idx) {
            _tem = _g0[_idx];
            _qx = T();
            _qy = T();
            _qz = T();
            for (int c = 1; c < Q<T>::nc; ++c) {
                T g = _g[Q<T>::IndexF(_idx, c)];
                _tem += g;
                _qx += Q<T>::cx[c]*g;
                _qy += Q<T>::cy[c]*g;
                _qz += Q<T>::cz[c]*g;
            }
            T coef = 1.0 - 0.5*_omegag;
            _qx = coef*(_qx - _tem*_ux);
            _qy = coef*(_qy - _tem*_uy);
            _qz = coef*(_qz - _tem*_uz);
        }

        //  Function of getting equilibrium of AD for 2D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _tem, T _ux, T _uy) {
            for (int c = 0; c < Q<T>::nc; ++c) {
                T ciu = Q<T>::cx[c]*_ux + Q<T>::cy[c]*_uy;
                _geq[c] = Q<T>::ei[c]*_tem*(1.0 + 3.0*ciu);
            }
        }

        //  Function of getting equilibrium of AD for 3D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _tem, T _ux, T _uy, T _uz) {
            for (int c = 0; c < Q<T>::nc; ++c) {
                T ciu = Q<T>::cx[c]*_ux + Q<T>::cy[c]*_uy + Q<T>::cz[c]*_uz;
                _geq[c] = Q<T>::ei[c]*_tem*(1.0 + 3.0*ciu);
            }
        }

        //  Function of applying external force of AD with natural convection for 2D
        template<class T, template<class>class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _tem0, T *_f, int _idx) {
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*P<T>::ei[c]*(P<T>::cx[c]*_gx + P<T>::cy[c]*_gy)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with natural convection for 3D
        template<class T, template<class>class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _gz, T _tem0, T *_f, int _idx) {
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*P<T>::ei[c]*(P<T>::cx[c]*_gx + P<T>::cy[c]*_gy + P<T>::cz[c]*_gz)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with heat exchange for 2D/3D
        template<class T, template<class>class Q>
        void ExternalForceHeatExchange(T _tem, T _beta, T *_g0, T *_g, int _idx) {
            T coef = _beta*(1.0 - _tem)/(1.0 + _beta);
            _g0[_idx] += Q<T>::ei[0]*coef;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += Q<T>::ei[c]*coef;
            }
        }

        //  Function of applying external force of AD with heat source for 2D/3D
        template<class T, template<class>class Q>
        void ExternalForceHeatSource(T _heatsource, T *_g0, T *_g, int _idx) {
            _g0[_idx] += Q<T>::ei[0]*_heatsource;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += Q<T>::ei[c]*_heatsource;
            }
        }

        //  Function of setting boundary condition set T of AD for D2Q9 along x edge
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetTAlongXEdge(Q<T>& _q, int _i, int _directionx, Fv _tembc, const T *_ux, const T *_uy, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T tem = _tembc(i + _q.offsetx, j + _q.offsety);
                        if (_directionx == -1) {
                            T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        } else if (_directionx == 1) {
                            T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D2Q9 along y edge
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetTAlongYEdge(Q<T>& _q, int _j, int _directiony, Fv _tembc, const T *_ux, const T *_uy, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T tem = _tembc(i + _q.offsetx, j + _q.offsety);
                        if (_directiony == -1) {
                            T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        } else if (_directiony == 1) {
                            T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D3Q15 along x face
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetTAlongXFace(Q<T>& _q, int _i, int _directionx, Fv _tembc, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T tem = _tembc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionx == -1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directionx == 1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(1.0 - 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D3Q15 along y face
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetTAlongYFace(Q<T>& _q, int _j, int _directiony, Fv _tembc, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T tem = _tembc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directiony == -1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directiony == 1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D3Q15 along z face
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetTAlongZFace(Q<T>& _q, int _k, int _directionz, Fv _tembc, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T tem = _tembc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionz == -1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            } else if (_directionz == 1) {
                                T tem0 = 6.0*(tem - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set q of AD for D2Q9 along x edge (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongXEdge(Q<T>& _q, int _i, int _directionx, Fv _qnbc, const T *_ux, const T *_uy, T _diffusivity, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directionx == -1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        } else if (_directionx == 1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D2Q9 along y edge (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongYEdge(Q<T>& _q, int _j, int _directiony, Fv _qnbc, const T *_ux, const T *_uy, T _diffusivity, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directiony == -1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        } else if (_directiony) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 along x face (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongXFace(Q<T>& _q, int _i, int _directionx, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionx == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directionx == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 along y face (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongYFace(Q<T>& _q, int _j, int _directiony, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directiony == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directiony == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }
            
        //  Function of setting boundary condition set q of AD for D3Q15 along z face (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongZFace(Q<T>& _q, int _k, int _directionz, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, Ff _bctype) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionz == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            } else if (_directionz == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*qn + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D2Q9 along x edge (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongXEdge(Q<T>& _q, int _i, int _directionx, Fv _qnbc, const T *_ux, const T *_uy, const T *_diffusivity, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directionx == -1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        } else if (_directionx == 1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D2Q9 along y edge (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongYEdge(Q<T>& _q, int _j, int _directiony, Fv _qnbc, const T *_ux, const T *_uy, const T *_diffusivity, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(i, j);
                        T qn = _qnbc(i + _q.offsetx, j + _q.offsety);
                        if (_directiony == -1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        } else if (_directiony == 1) {
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        }
                    }
                }
            }
        }
        
        //  Function of setting boundary condition set q of AD for D3Q15 along x face (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongXFace(Q<T>& _q, int _i, int _directionx, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, Ff _bctype) {
            int i = _i - _q.offsetx;
            if (0 <= i && i < _q.nx) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionx == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directionx == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                                _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }
            
        //  Function of setting boundary condition set q of AD for D3Q15 along y face (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongYFace(Q<T>& _q, int _j, int _directiony, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, Ff _bctype) {
            int j = _j - _q.offsety;
            if (0 <= j && j < _q.ny) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directiony == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            } else if (_directiony == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                                _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 along z face (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQAlongZFace(Q<T>& _q, int _k, int _directionz, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, Ff _bctype) {
            int k = _k - _q.offsetz;
            if (0 <= k && k < _q.nz) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, j, k);
                            T qn = _qnbc(i + _q.offsetx, j + _q.offsety, k + _q.offsetz);
                            if (_directionz == -1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            } else if (_directionz == 1) {
                                T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*qn + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                                _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                                _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                                _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            }
                        }
                    }
                }
            }
        }
    }

    namespace AD {
        //  Function of Update macro and Collide of force convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of force convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of natural convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc]; 
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of natural convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c); 
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of Brinkman and heat exchange for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and heat exchange for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {            
                    int idxf = Q<T>::IndexF(idx, c);      
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and force convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            bool _issave = false, T *_g = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;

                    if (_g) {
                        int offsetf = Q<T>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and force convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            bool _issave = false, T *_g = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;

                    if (_g) {
                        int offsetf = Q<T>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false, T *_g = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;

                    if (_g) {
                        int offsetf = Q<T>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false, T *_g = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;

                    if (_g) {
                        int offsetf = Q<T>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection with heat source for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvectionHeatSource(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, const T *_heatsource, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false, T *_g = nullptr
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                ExternalForceHeatSource<T, Q>(_heatsource[idx], _q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _uz[idx] = uz;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                    _qz[idx] = qz;

                    if (_g) {
                        int offsetf = Q<T>::nc*idx;
                        _g[offsetf] = _q.f0[idx];
                        for (int c = 1; c < Q<T>::nc; ++c) {
                            _g[offsetf + c] = _q.f[Q<T>::IndexF(idx, c)];
                        }
                    }
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of setting initial condition of AD for 2D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_tem, const T *_ux, const T *_uy) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _tem[idx], _ux[idx], _uy[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting initial condition of AD for 3D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_tem, const T *_ux, const T *_uy, const T *_uz) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _tem[idx], _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D2Q9
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetT(Q<T>& _q, Fv _tembc, const T *_ux, const T *_uy, Ff _bctype) {
            BoundaryConditionSetTAlongXEdge(_q, 0, -1, _tembc, _ux, _uy, _bctype);          //  On xmin
            BoundaryConditionSetTAlongXEdge(_q, _q.lx - 1, 1, _tembc, _ux, _uy, _bctype);   //  On xmax
            BoundaryConditionSetTAlongYEdge(_q, 0, -1, _tembc, _ux, _uy, _bctype);          //  On ymin
            BoundaryConditionSetTAlongYEdge(_q, _q.ly - 1, 1, _tembc, _ux, _uy, _bctype);   //  On ymax
        }

        //  Function of setting boundary condition set T of AD for D3Q15
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetT(Q<T>& _q, Fv _tembc, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            BoundaryConditionSetTAlongXFace(_q, 0, -1, _tembc, _ux, _uy, _uz, _bctype);         //  On xmin
            BoundaryConditionSetTAlongXFace(_q, _q.lx - 1, 1, _tembc, _ux, _uy, _uz, _bctype);  //  On xmax
            BoundaryConditionSetTAlongYFace(_q, 0, -1, _tembc, _ux, _uy, _uz, _bctype);         //  On ymin
            BoundaryConditionSetTAlongYFace(_q, _q.ly - 1, 1, _tembc, _ux, _uy, _uz, _bctype);  //  On ymax
            BoundaryConditionSetTAlongZFace(_q, 0, -1, _tembc, _ux, _uy, _uz, _bctype);         //  On zmin
            BoundaryConditionSetTAlongZFace(_q, _q.lz - 1, 1, _tembc, _ux, _uy, _uz, _bctype);  //  On zmax
        }
    
        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, T _diffusivity, Ff _bctype) {
            BoundaryConditionSetQAlongXEdge(_q, 0, -1, _qnbc, _ux, _uy, _diffusivity, _bctype);          //  On xmin
            BoundaryConditionSetQAlongXEdge(_q, _q.lx - 1, 1, _qnbc, _ux, _uy, _diffusivity, _bctype);   //  On xmax
            BoundaryConditionSetQAlongYEdge(_q, 0, -1, _qnbc, _ux, _uy, _diffusivity, _bctype);          //  On ymin
            BoundaryConditionSetQAlongYEdge(_q, _q.ly - 1, 1, _qnbc, _ux, _uy, _diffusivity, _bctype);   //  On ymax
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, Ff _bctype) {
            BoundaryConditionSetQAlongXFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On xmin
            BoundaryConditionSetQAlongXFace(_q, _q.lx - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On xmax
            BoundaryConditionSetQAlongYFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On ymin
            BoundaryConditionSetQAlongYFace(_q, _q.ly - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On ymax
            BoundaryConditionSetQAlongZFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On zmin
            BoundaryConditionSetQAlongZFace(_q, _q.lz - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On zmax
        }

        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_diffusivity, Ff _bctype) {
            BoundaryConditionSetQAlongXEdge(_q, 0, -1, _qnbc, _ux, _uy, _diffusivity, _bctype);          //  On xmin
            BoundaryConditionSetQAlongXEdge(_q, _q.lx - 1, 1, _qnbc, _ux, _uy, _diffusivity, _bctype);   //  On xmax
            BoundaryConditionSetQAlongYEdge(_q, 0, -1, _qnbc, _ux, _uy, _diffusivity, _bctype);          //  On ymin
            BoundaryConditionSetQAlongYEdge(_q, _q.ly - 1, 1, _qnbc, _ux, _uy, _diffusivity, _bctype);   //  On ymax
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, Ff _bctype) {
            BoundaryConditionSetQAlongXFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On xmin
            BoundaryConditionSetQAlongXFace(_q, _q.lx - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On xmax
            BoundaryConditionSetQAlongYFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On ymin
            BoundaryConditionSetQAlongYFace(_q, _q.ly - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On ymax
            BoundaryConditionSetQAlongZFace(_q, 0, -1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);         //  On zmin
            BoundaryConditionSetQAlongZFace(_q, _q.lz - 1, 1, _qnbc, _ux, _uy, _uz, _diffusivity, _bctype);  //  On zmax
        }
    }
}