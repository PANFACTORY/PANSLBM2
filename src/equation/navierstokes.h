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
        //  Navier-Stokes 2D    :   External force with immersed boundary method
        //*********************************************************************
//  要修正
        template<class T, template<class>class P>
        void ExternalForceImmersedBoundary(P<T>& _particle, int _nb, T* _Xbx, T* _Xby, T* _Ubx, T* _Uby) {
            //  Step0：境界点上の体積力を初期化する
            T *Gbx = new T[_nb], *Gby = new T[_nb];
            for (int k = 0; k < _nb; k++) {
                Gbx[k] = T();   Gby[k] = T();
            }

            T *Gx = new T[_particle.np], *Gy = new T[_particle.np], *ux = new T[_particle.np], *uy = new T[_particle.np];
            for (int i = 0; i < _particle.np; i++) {
                Gx[i] = T();    Gy[i] = T();
                ux[i] = T();    uy[i] = T();
            }

            for (int l = 0; l < 5; l++) {
                //  格子点 → 境界点
                for (int k = 0; k < _nb; k ++) {
                    //  Step3：格子点上の流速 → 境界点上の流速
                    T ubx = T(), uby = T();
                    for (int i = 0; i < _particle.nx; i++) {
                        for (int j = 0; j < _particle.ny; j++) {
                            int ij = _particle.ny*i + j;
                            T w = 1.0, rx = fabs(_Xbx[k] - i), ry = fabs(_Xby[k] - j);
                            w *= rx < 1.0 ? (3.0 - 2.0*rx + sqrt(1.0 + 4.0*rx - 4.0*rx*rx))/8.0 : (1.0 <= rx && rx < 2.0 ? (5.0 - 2.0*rx - sqrt(-7.0 + 12.0*rx - 4.0*rx*rx))/8.0 : T());
                            w *= ry < 1.0 ? (3.0 - 2.0*ry + sqrt(1.0 + 4.0*ry - 4.0*ry*ry))/8.0 : (1.0 <= ry && ry < 2.0 ? (5.0 - 2.0*ry - sqrt(-7.0 + 12.0*ry - 4.0*ry*ry))/8.0 : T());
                            ubx += ux[ij]*w*pow(_particle.dx, P<T>::nd);    uby += uy[ij]*w*pow(_particle.dx, P<T>::nd);
                        }
                    }

                    //  Step4：境界点上の流速 → 境界点上の体積力
                    Gbx[k] += _particle.dx*(_Ubx[k] - ubx);     Gby[k] += _particle.dx*(_Uby[k] - uby);
                }

                //  境界点 → 格子点
                for (int i = 0; i < _particle.nx; i++) {
                    for (int j = 0; j < _particle.ny; j++) {
                        int ij = _particle.ny*i + j;

                        //  Step1：境界点上の体積力 → 格子点上の体積力
                        Gx[ij] = T();
                        Gy[ij] = T();
                        for (int k = 0; k < _nb; k++) {
                            T w = 1.0, rx = fabs(_Xbx[k] - i), ry = fabs(_Xby[k] - j);
                            w *= rx < 1.0 ? (3.0 - 2.0*rx + sqrt(1.0 + 4.0*rx - 4.0*rx*rx))/8.0 : (1.0 <= rx && rx < 2.0 ? (5.0 - 2.0*rx - sqrt(-7.0 + 12.0*rx - 4.0*rx*rx))/8.0 : T());
                            w *= ry < 1.0 ? (3.0 - 2.0*ry + sqrt(1.0 + 4.0*ry - 4.0*ry*ry))/8.0 : (1.0 <= ry && ry < 2.0 ? (5.0 - 2.0*ry - sqrt(-7.0 + 12.0*ry - 4.0*ry*ry))/8.0 : T());
                            Gx[ij] += Gbx[k]*w*pow(_particle.dx, P<T>::nd);     Gy[ij] += Gby[k]*w*pow(_particle.dx, P<T>::nd);
                        }

                        //  Step2：格子点上の体積力 → 格子点上の流速
                        ux[ij] += _particle.dx*Gx[ij];  uy[ij] += _particle.dx*Gy[ij];
                    }
                }
            }
            
            //  Step5：速度分布関数を更新する
            for (int i = 0; i < _particle.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _particle.ft[j] += 3.0*_particle.dx*P<T>::ei[j]*(P<T>::cx[j]*Gx[i] + P<T>::cy[j]*Gy[i]);
                }
            }

            delete[] Gbx;   delete[] Gby;
            delete[] Gx;    delete[] Gy;
            delete[] ux;    delete[] uy;
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