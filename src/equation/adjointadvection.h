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
        void ExternalForceHeatgeneration(P<T>& _particle, T* _tem, T* _beta) {
            for (int i = 0; i < _particle.np; i++) {
                _tem[i] = (_tem[i] + _particle.dx*_beta[i])/(1.0 + _particle.dx*_beta[i]);
            }
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint navier-stokes 2D    :   External force with heat exchange
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceHeatexchange(T _diffusivity, P<T>& _particlef, T* _rho, T* _ux, T* _uy, T* _tem, T* _iqx, T* _iqy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_diffusivity*_particlef.dt/(_particlef.dx*_particlef.dx) + 0.5);
            for (int i = 0; i < _particlef.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_tem[i]*omega*((P<T>::cx[j] - _ux[i])*_iqx[i] + (P<T>::cy[j] - _uy[i])*_iqy[i])/_rho[i];
                }
            }
        }
    }

    namespace AAD {
        //*********************************************************************
        //  Adjoint advection 2D    :   Update macroscopic values
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _item, T* _iqx, T* _iqy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _item[i] = T();
                _iqx[i] = T();
                _iqy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _item[i] += P<T>::ei[j]*_particle.ft[j][i];
                    _iqx[i] += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    _iqy[i] += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _diffusivity, P<T>& _particle, T* _ux, T* _uy, T* _item, T* _iqx, T* _iqy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                T feq = _item[i] + 3.0*(_ux[i]*_iqx[i] + _uy[i]*_iqy[i]);
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   External force
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceHeatexchange(P<T>& _particle, T* _item, T* _beta) {
            for (int i = 0; i < _particle.np; i++) {
                _item[i] = (_item[i] - _beta[i])/(1.0 + _particle.dx*_beta[i]);
            }
        }

        //*********************************************************************
        //  Adjoint advection 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _ux, T _uy, T _item, T _iqx, T _iqy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                _particle.ft[j][_i] = _item + 3.0*(_ux*_iqx + _uy*_iqy);
            }
        }
    }
}