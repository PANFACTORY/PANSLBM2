//*****************************************************************************
//  Title       :   src/adjointnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>

namespace PANSLBM2 {
    namespace NS2 {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _alpha) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T rho = T(), ux = T(), uy = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    rho += _particle.ft[j][i];
                    ux += P<T>::ci[j][0]*_particle.ft[j][i];
                    uy += P<T>::ci[j][1]*_particle.ft[j][i];
                }
                ux /= rho;
                uy /= rho;
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] += -3.0*_particle.dx*P<T>::ei[j]*_alpha[i]*(P<T>::ci[j][0]*ux + P<T>::ci[j][1]*uy);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Check convergence
        //*********************************************************************
        template<class T, template<class>class P>
        bool CheckConvergence(P<T>& _particle, T _eps, T* _uxt, T* _uyt, T* _uxtm1, T* _uytm1) {
            T unorm = T(), dunorm = T();
            for (int i = 0; i < _particle.np; i++) {
                unorm += _uxt[i]*_uxt[i] + _uyt[i]*_uyt[i];
                dunorm += (_uxt[i] - _uxtm1[i])*(_uxt[i] - _uxtm1[i]) + (_uyt[i] - _uytm1[i])*(_uyt[i] - _uytm1[i]);
            }
            return sqrt(dunorm/unorm) < _eps ? true : false;
        }
    }

    namespace ANS2 {
        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Update macroscopic values, rho*, u*, v*, w*
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _q, T* _vx, T* _vy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _q[i] = T();
                _vx[i] = T();
                _vy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::ci[j][0]*_ux[i] + P<T>::ci[j][1]*_uy[i]; 
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    _q[i] += _particle.ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                    _vx[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::ci[j][0] + 3.0*ciu*P<T>::ci[j][0] - _ux[i]);
                    _vy[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::ci[j][1] + 3.0*ciu*P<T>::ci[j][1] - _uy[i]);
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _ux, T* _uy, T* _q, T* _vx, T* _vy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T feq = _q[i] + 3.0*(_vx[i]*(P<T>::ci[j][0] - _ux[i]) + _vy[i]*(P<T>::ci[j][1] - _uy[i]));
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P, class ...Ts>
        void ExternalForce(P<T>& _particle, T* _alpha, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T mx = T(), my = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    mx += P<T>::ei[j]*P<T>::ci[j][0]*_particle.ft[j][i];
                    my += P<T>::ei[j]*P<T>::ci[j][1]*_particle.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] -= 3.0*_alpha[i]*(mx*(P<T>::ci[j][0] - _ux[i]) + my*(P<T>::ci[j][1] - _uy[i]))/_rho[i];
                }
            }
        }
    }  
}