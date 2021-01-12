//*****************************************************************************
//  Title       :   src/equation/adjointadvection.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/11/02
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>

namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection   :   External force with heat generation
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceHeatgeneration(P<T>& _particle, T* _beta) {
            for (int i = 0; i < _particle.np; i++) {
                T tmpT = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    tmpT += _particle.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] += _particle.dx*P<T>::ei[j]*_beta[i]*(1.0 - tmpT)/(1.0 + _particle.dx*_beta[i]);
                }
            }
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint navier-stokes 2D    :   External force with heat exchange
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceHeatexchange(T _diffusivity, P<T>& _particlef, T* _rho, T* _ux, T* _uy, T* _T, T* _qtildex, T* _qtildey) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_diffusivity*_particlef.dt/(_particlef.dx*_particlef.dx) + 0.5);
            for (int i = 0; i < _particlef.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_T[i]*omega*((P<T>::cx[j] - _ux[i])*_qtildex[i] + (P<T>::cy[j] - _uy[i])*_qtildey[i])/_rho[i];
                }
            }
        }
    }

    namespace AAD {
        //*********************************************************************
        //  Adjoint advection 2D    :   Update macroscopic values
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _Ttilde, T* _qtildex, T* _qtildey) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _Ttilde[i] = T();
                _qtildex[i] = T();
                _qtildey[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _Ttilde[i] += P<T>::ei[j]*_particle.ft[j][i];
                    _qtildex[i] += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    _qtildey[i] += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _diffusivity, P<T>& _particle, T* _ux, T* _uy, T* _Ttilde, T* _qtildex, T* _qtildey) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                T feq = _Ttilde[i] + 3.0*(_ux[i]*_qtildex[i] + _uy[i]*_qtildey[i]);
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   External force
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceHeatexchange(P<T>& _particle, T* _beta) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T Ttilde = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    Ttilde += P<T>::ei[j]*_particle.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] -= _beta[i]*(1.0 + _particle.dx*Ttilde)/(1.0 + _particle.dx*_beta[i]);
                }
            }
        }
    }
}