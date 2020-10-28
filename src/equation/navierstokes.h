//*****************************************************************************
//  Title       :   src/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes   :   Update macroscopic values, rho, u, v, w
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void UpdateMacro(P<T>& _particle, T* _rho, Ts ..._u) {
            assert(P<T>::nd == sizeof...(Ts));
            T* u[P<T>::nd] = { _u... };
            for (int i = 0; i < _particle.np; i++) {
                _rho[i] = T();
                for (int d = 0; d < P<T>::nd; d++) {
                    u[d][i] = T();
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _rho[i] += _particle.ft[j][i];
                    for (int d = 0; d < P<T>::nd; d++) {
                        u[d][i] += P<T>::ci[j][d]*_particle.ft[j][i];
                    }
                }
                for (int d = 0; d < P<T>::nd; d++) {
                    u[d][i] /= _rho[i];
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes   :   Collision term
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void Collision(T _viscosity, P<T>& _particle, T* _rho, Ts ..._u) {
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            T* u[P<T>::nd] = { _u... };
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = T(), uu = T();
                    for (int d = 0; d < P<T>::nd; d++) {
                        ciu += P<T>::ci[j][d]*u[d][i];
                        uu += u[d][i]*u[d][i];
                    }
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes   :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void InitialCondition(int _i, P<T>& _particle, T _rho, Ts ..._u) {
            assert(P<T>::nd == sizeof...(Ts) && 0 <= _i && _i < _particle.np);
            T us[P<T>::nd] = { _u... };
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T(), uu = T();
                for (int d = 0; d < P<T>::nd; d++) {
                    ciu += P<T>::ci[j][d]*us[d];
                    uu += us[d]*us[d];
                }
                _particle.ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }
    }
}