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
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _alphax, T* _alphay) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T rho = T(), ux = T(), uy = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    rho += _particle.ft[j][i];
                    ux += P<T>::cx[j]*_particle.ft[j][i];
                    uy += P<T>::cy[j]*_particle.ft[j][i];
                }
                ux /= rho;
                uy /= rho;
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] += -3.0*_particle.dx*P<T>::ei[j]*(_alphax[i]*P<T>::cx[j]*ux + _alphay[i]*P<T>::cy[j]*uy);
                }
            }
        }

        //*********************************************************************
        //  Navier-Stokes 3D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _alphax, T* _alphay, T* _alphaz) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                T rho = T(), ux = T(), uy = T(), uz = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    rho += _particle.ft[j][i];
                    ux += P<T>::cx[j]*_particle.ft[j][i];
                    uy += P<T>::cy[j]*_particle.ft[j][i];
                    uz += P<T>::cz[j]*_particle.ft[j][i];
                }
                ux /= rho;
                uy /= rho;
                uz /= rho;
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] += -3.0*_particle.dx*P<T>::ei[j]*(_alphax[i]*P<T>::cx[j]*ux + _alphay[i]*P<T>::cy[j]*uy + _alphaz[i]*P<T>::cz[j]*uz);
                }
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

        //*********************************************************************
        //  Navier-Stokes 3D    :   Check convergence
        //*********************************************************************
        template<class T, template<class>class P>
        bool CheckConvergence(P<T>& _particle, T _eps, T* _uxt, T* _uyt, T* _uzt, T* _uxtm1, T* _uytm1, T* _uztm1) {
            assert(P<T>::nd == 3);
            T unorm = T(), dunorm = T();
            for (int i = 0; i < _particle.np; i++) {
                unorm += _uxt[i]*_uxt[i] + _uyt[i]*_uyt[i] + _uzt[i]*_uzt[i];
                dunorm += (_uxt[i] - _uxtm1[i])*(_uxt[i] - _uxtm1[i]) + (_uyt[i] - _uytm1[i])*(_uyt[i] - _uytm1[i]) + (_uzt[i] - _uztm1[i])*(_uzt[i] - _uztm1[i]);
            }
            return sqrt(dunorm/unorm) < _eps ? true : false;
        }
    }

    namespace ANS {
        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Update macroscopic values, rho*, u*, v*
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _q, T* _vx, T* _vy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                _q[i] = T();
                _vx[i] = T();
                _vy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i]; 
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i];
                    _q[i] += _particle.ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                    _vx[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cx[j] + 3.0*ciu*P<T>::cx[j] - _ux[i]);
                    _vy[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cy[j] + 3.0*ciu*P<T>::cy[j] - _uy[i]);
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 3D    :   Update macroscopic values, rho*, u*, v*, w*
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _particle, T* _rho, T* _ux, T* _uy, T* _uz, T* _q, T* _vx, T* _vy, T* _vz) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                _q[i] = T();
                _vx[i] = T();
                _vy[i] = T();
                _vz[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    T ciu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i] + P<T>::cz[j]*_uz[i]; 
                    T uu = _ux[i]*_ux[i] + _uy[i]*_uy[i] + _uz[i]*_uz[i];
                    _q[i] += _particle.ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                    _vx[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cx[j] + 3.0*ciu*P<T>::cx[j] - _ux[i]);
                    _vy[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cy[j] + 3.0*ciu*P<T>::cy[j] - _uy[i]);
                    _vz[i] += _particle.ft[j][i]*P<T>::ei[j]*(P<T>::cz[j] + 3.0*ciu*P<T>::cz[j] - _uz[i]);
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
                    T feq = _q[i] + 3.0*(_vx[i]*(P<T>::cx[j] - _ux[i]) + _vy[i]*(P<T>::cy[j] - _uy[i]));
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 3D    :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _viscosity, P<T>& _particle, T* _ux, T* _uy, T* _uz, T* _q, T* _vx, T* _vy, T* _vz) {
            assert(P<T>::nd == 3);
            T omega = 1.0/(3.0*_viscosity*_particle.dt/(_particle.dx*_particle.dx) + 0.5);
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T feq = _q[i] + 3.0*(_vx[i]*(P<T>::cx[j] - _ux[i]) + _vy[i]*(P<T>::cy[j] - _uy[i]) + _vz[i]*(P<T>::cz[j] - _uz[i]));
                    _particle.ftp1[j][i] = (1.0 - omega)*_particle.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _alphax, T* _alphay, T* _rho, T* _ux, T* _uy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T mx = T(), my = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    mx += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    my += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] -= 3.0*(_alphax[i]*mx*(P<T>::cx[j] - _ux[i]) + _alphay[i]*my*(P<T>::cy[j] - _uy[i]))/_rho[i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 3D    :   External force with Brinkman model
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceBrinkman(P<T>& _particle, T* _alphax, T* _alphay, T* _alphaz, T* _rho, T* _ux, T* _uy, T* _uz) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                T mx = T(), my = T(), mz = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    mx += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    my += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                    mz += P<T>::ei[j]*P<T>::cz[j]*_particle.ft[j][i];
                }
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j][i] -= 3.0*(_alphax[i]*mx*(P<T>::cx[j] - _ux[i]) + _alphay[i]*my*(P<T>::cy[j] - _uy[i]) + _alphaz[i]*mz*(P<T>::cz[j] - _uz[i]))/_rho[i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 2D    :   Update sensitivity
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateSensitivity(P<T>& _particle, T* _ux, T* _uy, T* _sensitivity) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _particle.np; i++) {
                T mx = T(), my = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    mx += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    my += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                }
                _sensitivity[i] += 3.0*(mx*_ux[i] + my*_uy[i]);
            }
        }

        //*********************************************************************
        //  Adjoint Navier-Stokes 3D    :   Update sensitivity
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateSensitivity(P<T>& _particle, T* _ux, T* _uy, T* _uz, T* _sensitivity) {
            assert(P<T>::nd == 3);
            for (int i = 0; i < _particle.np; i++) {
                T mx = T(), my = T(), mz = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    mx += P<T>::ei[j]*P<T>::cx[j]*_particle.ft[j][i];
                    my += P<T>::ei[j]*P<T>::cy[j]*_particle.ft[j][i];
                    mz += P<T>::ei[j]*P<T>::cz[j]*_particle.ft[j][i];
                }
                _sensitivity[i] += 3.0*(mx*_ux[i] + my*_uy[i] + mz*_uz[i]);
            }
        }
    }  
}