//*****************************************************************************
//  Title       :   src/equation/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection 2D    :   Update macroscopic values, T and flux
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _temperature, T* _qx, T* _qy, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _temperature[i] = T();
                _qx[i] = T();
                _qy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _temperature[i] += _particle.ft[j][i];
                    _qx[i] += P<T>::cx[j]*_particle.ft[j][i];
                    _qy[i] += P<T>::cy[j]*_particle.ft[j][i];
                }
                _qx[i] -= _temperature[i]*_ux[i];
                _qy[i] -= _temperature[i]*_uy[i];
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Update macroscopic values, T and flux
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _temperature, T* _qx, T* _qy, T* _qz, T* _ux, T* _uy, T* _uz) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                _temperature[i] = T();
                _qx[i] = T();
                _qy[i] = T();
                _qz[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _temperature[i] += _particle.ft[j][i];
                    _qx[i] += P<T>::cx[j]*_particle.ft[j][i];
                    _qy[i] += P<T>::cy[j]*_particle.ft[j][i];
                    _qz[i] += P<T>::cz[j]*_particle.ft[j][i];
                }
                _qx[i] -= _temperature[i]*_ux[i];
                _qy[i] -= _temperature[i]*_uy[i];
                _qz[i] -= _temperature[i]*_uz[i];
            }
        }

        //*********************************************************************
        //  Advection 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _diffusivity, P<T>& _particle, T* _temperature, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_temperature[i]*(1.0 + 3.0*ciu);
                }
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _diffusivity, P<T>& _particle, T* _temperature, T* _ux, T* _uy, T* _uz) {
            assert(P<T>::nd == 3);
            T omega = 1.0/(3.0*_diffusivity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i] + P<T>::cz[j]*_uz[i];
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*P<T>::ei[j]*_temperature[i]*(1.0 + 3.0*ciu);
                }
            }
        }

        //*********************************************************************
        //  Advection 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _temperature, T _ux, T _uy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy;
                _particle.ft[j][_i] = P<T>::ei[j]*_temperature*(1.0 + 3.0*ciu);
            }
        }

        //*********************************************************************
        //  Advection 3D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _temperature, T _ux, T _uy, T _uz) {
            assert(P<T>::nd == 3 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy + P<T>::cz[j]*_uz;
                _particle.ft[j][_i] = P<T>::ei[j]*_temperature*(1.0 + 3.0*ciu);
            }
        }
    }

    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with natural convection
        //*********************************************************************
        template<class T, template<class>class P, template<class>class Q>
        void ExternalForceNaturalConvection(T _gx, T _gy, T _temperature0, P<T>& _particlef, Q<T>& _particleg) {
            assert(P<T>::nd == 2 && Q<T>::nd == 2 && _particlef.nx == _particleg.nx && _particlef.ny == _particleg.ny);
            for (int i = 0; i < _particlef.np; i++) {
                T temperature = T();
                for (int j = 0; j < Q<T>::nc; j++) {
                    temperature += _particleg.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_particlef.dx*P<T>::ei[j]*(Q<T>::cx[j]*_gx + Q<T>::cy[j]*_gy)*(temperature - _temperature0);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   External force with natural convection
        //*********************************************************************
        template<class T, template<class>class P, template<class>class Q>
        void ExternalForceNaturalConvection(T _gx, T _gy, T _gz, T _temperature0, P<T>& _particlef, Q<T>& _particleg) {
            assert(P<T>::nd == 3 && Q<T>::nd == 3 && _particlef.nx == _particleg.nx && _particlef.ny == _particleg.ny && _particlef.nz == _particleg.nz);
            for (int i = 0; i < _particlef.np; i++) {
                T temperature = T();
                for (int j = 0; j < Q<T>::nc; j++) {
                    temperature += _particleg.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particlef.ft[j][i] += 3.0*_particlef.dx*P<T>::ei[j]*(Q<T>::cx[j]*_gx + Q<T>::cy[j]*_gy + Q<T>::cz[j]*_gz)*(temperature - _temperature0);
                }
            }
        }
    }
}