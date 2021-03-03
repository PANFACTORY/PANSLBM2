//*****************************************************************************
//  Title       :   src/equation/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   Update macroscopic values, rho, ux, uy
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _rho[i] = T();
                _ux[i] = T();
                _uy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _rho[i] += _particle.ft[j][i];
                    _ux[i] += P<T>::cx[j]*_particle.ft[j][i];
                    _uy[i] += P<T>::cy[j]*_particle.ft[j][i];
                }
                _ux[i] /= _rho[i];
                _uy[i] /= _rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   Update macroscopic values, rho, ux, uy, uz
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _uz) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                _rho[i] = T();
                _ux[i] = T();
                _uy[i] = T();
                _uz[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _rho[i] += _particle.ft[j][i];
                    _ux[i] += P<T>::cx[j]*_particle.ft[j][i];
                    _uy[i] += P<T>::cy[j]*_particle.ft[j][i];
                    _uz[i] += P<T>::cz[j]*_particle.ft[j][i];
                }
                _ux[i] /= _rho[i];
                _uy[i] /= _rho[i];
                _uz[i] /= _rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i];
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _uz) {
            assert(P<T>::nd == 3);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i] + P<T>::cz[j]*_uz[i];
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i] + _uz[i]*_uz[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Collision term with multiple-relaxation-time model
        //*********************************************************************
        template<class T, template<class>class P>
        void CollisionMRT(T _s1, T _s2, T _s4, T _s6, T _s7, T _s8, P<T>& _p, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);

            const static T M[P<T>::nc][P<T>::nc] = {
                {   1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0,    1.0 },
                {   -4.0,   -1.0,   -1.0,   -1.0,   -1.0,   2.0,    2.0,    2.0,    2.0 },
                {   4.0,    -2.0,   -2.0,   -2.0,   -2.0,   1.0,    1.0,    1.0,    1.0 },
                {   0.0,    1.0,    0.0,    -1.0,   0.0,    1.0,    -1.0,   -1.0,   1.0 },
                {   0.0,    -2.0,   0.0,    2.0,    0.0,    1.0,    -1.0,   -1.0,   1.0 },
                {   0.0,    0.0,    1.0,    0.0,    -1.0,   1.0,    1.0,    -1.0,   -1.0    },
                {   0.0,    0.0,    -2.0,   0.0,    2.0,    1.0,    1.0,    -1.0,   -1.0    },
                {   0.0,    1.0,    -1.0,   1.0,    -1.0,   0.0,    0.0,    0.0,    0.0 },
                {   0.0,    0.0,    0.0,    0.0,    0.0,    1.0,    -1.0,   1.0,    -1.0    }
            };

            const static T iMS[P<T>::nc][P<T>::nc] = {
                {   0.0,    -_s1/9.0,   _s2/9.0,    0.0,    0.0,        0.0,    0.0,        0.0,        0.0 },
                {   0.0,    -_s1/36.0,  -_s2/18.0,  0.0,    -_s4/6.0,   0.0,    0.0,        _s7/4.0,    0.0 },
                {   0.0,    -_s1/36.0,  -_s2/18.0,  0.0,    0.0,        0.0,    -_s6/6.0,   -_s7/4.0,   0.0 },
                {   0.0,    -_s1/36.0,  -_s2/18.0,  0.0,    _s4/6.0,    0.0,    0.0,        _s7/4.0,    0.0 },
                {   0.0,    -_s1/36.0,  -_s2/18.0,  0.0,    0.0,        0.0,    _s6/6.0,    -_s7/4.0,   0.0 },
                {   0.0,    _s1/18.0,   _s2/36.0,   0.0,    _s4/12.0,   0.0,    _s6/12.0,   0.0,        _s8/4.0 },
                {   0.0,    _s1/18.0,   _s2/36.0,   0.0,    -_s4/12.0,  0.0,    _s6/12.0,   0.0,        -_s8/4.0    },
                {   0.0,    _s1/18.0,   _s2/36.0,   0.0,    -_s4/12.0,  0.0,    -_s6/12.0,  0.0,        _s8/4.0 },
                {   0.0,    _s1/18.0,   _s2/36.0,   0.0,    _s4/12.0,   0.0,    -_s6/12.0,  0.0,        -_s8/4.0    }
            };

            for (int i = 0; i < _p.np; i++) {
                T meq[P<T>::nc] = {
                    _rho[i],
                    -2.0*_rho[i] + _rho[i]*(pow(_ux[i], 2.0) + pow(_uy[i], 2.0)),
                    9.0*_rho[i]*pow(_ux[i]*_uy[i], 2.0) - 3.0*_rho[i]*(pow(_ux[i], 2.0) + pow(_uy[i], 2.0)) + _rho[i],
                    _rho[i]*_ux[i],
                    3.0*_rho[i]*pow(_ux[i], 3.0) - _rho[i]*_ux[i],
                    _rho[i]*_uy[i],
                    3.0*_rho[i]*pow(_uy[i], 3.0) - _rho[i]*_uy[i],
                    _rho[i]*(pow(_ux[i], 2.0) - pow(_uy[i], 2.0)),
                    _rho[i]*_ux[i]*_uy[i]
                };

                T mmmeq[P<T>::nc];
                for (int j = 0; j < P<T>::nc; j++) {
                    mmmeq[j] = -meq[j];
                    for (int k = 0; k < P<T>::nc; k++) {
                        mmmeq[j] += M[j][k]*_p.ft[k][i];
                    }
                }

                for (int j = 0; j < P<T>::nc; j++) {
                    _p.ftp1[j][i] = _p.ft[j][i];
                    for (int k = 0; k < P<T>::nc; k++) {
                        _p.ftp1[j][i] -= iMS[j][k]*mmmeq[k];
                    }
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Collision term with two-relaxation-time model
        //*********************************************************************
        template<class T, template<class>class P>
        void CollisionTRT(T _viscosity, T _lambda, P<T>& _particle, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            T omegap = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            T omegam = 1.0/(_lambda/(3.0*_viscosity) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                double feq[P<T>::nc];
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i];
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    feq[j] = P<T>::ei[j]*_rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                }
                double dfp[P<T>::nc] = {
                    _particle.ft[0][i] - feq[0],
                    0.5*(_particle.ft[1][i] + _particle.ft[3][i]) - 0.5*(feq[1] + feq[3]),
                    0.5*(_particle.ft[2][i] + _particle.ft[4][i]) - 0.5*(feq[2] + feq[4]),
                    0.5*(_particle.ft[3][i] + _particle.ft[1][i]) - 0.5*(feq[3] + feq[1]),
                    0.5*(_particle.ft[4][i] + _particle.ft[2][i]) - 0.5*(feq[4] + feq[2]),
                    0.5*(_particle.ft[5][i] + _particle.ft[7][i]) - 0.5*(feq[5] + feq[7]),
                    0.5*(_particle.ft[6][i] + _particle.ft[8][i]) - 0.5*(feq[6] + feq[8]),
                    0.5*(_particle.ft[7][i] + _particle.ft[5][i]) - 0.5*(feq[7] + feq[5]),
                    0.5*(_particle.ft[8][i] + _particle.ft[6][i]) - 0.5*(feq[8] + feq[6])
                };
                double dfm[P<T>::nc] = {
                    T(),
                    0.5*(_particle.ft[1][i] - _particle.ft[3][i]) - 0.5*(feq[1] - feq[3]),
                    0.5*(_particle.ft[2][i] - _particle.ft[4][i]) - 0.5*(feq[2] - feq[4]),
                    0.5*(_particle.ft[3][i] - _particle.ft[1][i]) - 0.5*(feq[3] - feq[1]),
                    0.5*(_particle.ft[4][i] - _particle.ft[2][i]) - 0.5*(feq[4] - feq[2]),
                    0.5*(_particle.ft[5][i] - _particle.ft[7][i]) - 0.5*(feq[5] - feq[7]),
                    0.5*(_particle.ft[6][i] - _particle.ft[8][i]) - 0.5*(feq[6] - feq[8]),
                    0.5*(_particle.ft[7][i] - _particle.ft[5][i]) - 0.5*(feq[7] - feq[5]),
                    0.5*(_particle.ft[8][i] - _particle.ft[6][i]) - 0.5*(feq[8] - feq[6])
                };
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ftp1[j][i] = _particle.ft[j][i] - omegap*dfp[j] - omegam*dfm[j];
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _rho, T _ux, T _uy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy; 
                T uu = _ux*_ux + _uy*_uy;
                _particle.ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _rho, T _ux, T _uy, T _uz) {
            assert(P<T>::nd == 3 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy + P<T>::cz[j]*_uz; 
                T uu = _ux*_ux + _uy*_uy + _uz*_uz;
                _particle.ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }
    }
}