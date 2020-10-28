//*****************************************************************************
//  Title       :   src/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection   :   Update macroscopic values, rho, u, v, w
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _temperature) {
            for (int i = 0; i < _particle.np; i++) {
                _temperature[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _temperature[i] += _particle.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Advection   :   Collision term
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void Collision(T _diffusivity, P<T>& _particle, T* _temperature, Ts ..._u) {
            T omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            T* u[P<T>::nd] = { _u... };
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = T();
                    for (int d = 0; d < P<T>::nd; d++) {
                        ciu += P<T>::ci[j][d]*u[d][i];
                    }
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_temperature[i]*(1.0 + 3.0*ciu);
                }
            }
        }

        //*********************************************************************
        //  Advection   :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void InitialCondition(int _i, P<T>& _particle, T _temperature, Ts ..._u) {
            assert(P<T>::nd == sizeof...(Ts) && 0 <= _i && _i < _particle.np);
            T us[P<T>::nd] = { _u... };
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                for (int d = 0; d < P<T>::nd; d++) {
                    ciu += P<T>::ci[j][d]*us[d];
                }
                _particle.ft[j][_i] = P<T>::ei[j]*_temperature*(1.0 + 3.0*ciu);
            }
        }
    }

    namespace NS {
        //*********************************************************************
        //  Navier-Stokes   :   External force with natural convection
        //*********************************************************************
        template<class T, template<class>class P, template<class>class Q>
        void ExternalForceNaturalConvection(T _rhog, T _temperature0, P<T>& _particlef, Q<T>& _particleg) {
            for (int i = 0; i < _particlef.np; i++) {
                T temperature = T();
                for (int j = 0; j < Q<T>::nc; j++) {
                    temperature += _particleg.ft[j][i];
                }            
                T rhogdt = _rhog*(temperature - _temperature0);
                for (int j = 0; j < P<T>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_particlef.dx*P<T>::ei[j]*Q<T>::ci[j][1]*rhogdt;
                }
            }
        }
    }
}