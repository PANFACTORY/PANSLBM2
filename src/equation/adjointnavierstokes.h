//*****************************************************************************
//  Title       :   src/equation/adjointnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/28
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T *_rho, T* _ux, T* _uy, T *_alphax, T *_alphay) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _ux[i] /= 1.0 + _particle.dx*_alphax[i]/_rho[i];
                _uy[i] /= 1.0 + _particle.dx*_alphay[i]/_rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Check convergence
        //*********************************************************************
        template<class T, template<class>class P>
        bool CheckConvergence(P<T>& _particle, T _eps, T* _uxt, T* _uyt, T* _uxtm1, T* _uytm1) {
            assert(P<T>::nd == 2);
            T unorm = T(), dunorm = T();
            for (int i = 0; i < _particle.np; i++) {
                unorm += _uxt[i]*_uxt[i] + _uyt[i]*_uyt[i];
                dunorm += (_uxt[i] - _uxtm1[i])*(_uxt[i] - _uxtm1[i]) + (_uyt[i] - _uytm1[i])*(_uyt[i] - _uytm1[i]);
            }
            return sqrt(dunorm/unorm) < _eps ? true : false;
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Update macroscopic values, rho*, u*, v*
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _ip, T* _iux, T* _iuy, T* _imx, T* _imy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _ip[i] = T();
                _iux[i] = T();
                _iuy[i] = T();
                _imx[i] = T();
                _imy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i]; 
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    _ip[i] += _particle.ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                    _iux[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cx[j] + 3.0*ciu*P<T>::cx[j] - _ux[i]);
                    _iuy[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cy[j] + 3.0*ciu*P<T>::cy[j] - _uy[i]);
                    _imx[i] += _particle.ft[j][i]*P<T>::ei[j]*P<T>::cx[j];
                    _imy[i] += _particle.ft[j][i]*P<T>::ei[j]*P<T>::cy[j];
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _ux, T* _uy, T* _ip, T* _iux, T* _iuy) {
            assert(P<T>::nd == 2);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T feq = _ip[i] + 3.0*(_iux[i]*(P<T>::cx[j] - _ux[i]) + _iuy[i]*(P<T>::cy[j] - _uy[i]));
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _rho, T* _iux, T* _iuy, T* _imx, T* _imy, T* _alphax, T* _alphay) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _imx[i] /= 1.0 + _particle.dx*_alphax[i]/_rho[i];
                _imy[i] /= 1.0 + _particle.dx*_alphay[i]/_rho[i];
                _iux[i] -= _particle.dx*_alphax[i]*_imx[i]/_rho[i];
                _iuy[i] -= _particle.dx*_alphay[i]*_imy[i]/_rho[i];
            }
        }

        //*********************************************************************
        //  Navier-Stokes 2D    :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _particle, T _ux, T _uy, T _ip, T _iux, T _iuy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _particle.np);
            for (int j = 0; j < P<T>::nc; j++) {
                _particle.ft[j][_i] = _ip[_i] + 3.0*(_iux[_i]*(P<T>::cx[j] - _ux[_i]) + _iuy[_i]*(P<T>::cy[j] - _uy[_i]));
            }
        }
    }  
}