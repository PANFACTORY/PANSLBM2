//*****************************************************************************
//  Title       :   src/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace NS2 {
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
        //  Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
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
        //  Navier-Stokes 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void InitialCondition(int _i, P<T>& _particle, T _rho, T _ux, T _uy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy; 
                T uu = _ux*_ux + _uy*_uy;
                _particle.ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }
    }
}