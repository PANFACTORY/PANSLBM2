//*****************************************************************************
//  Title       :   src/lbmns.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/10
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "lbmbase.h"


namespace PANSLBM2 {
    template<class T>
    class LBMNS : public LBMBASE<T> {
public:
        LBMNS(int _nx, int _ny, T _viscosity);
        ~LBMNS() {};

        void Inlet(T _u, T _v);
        void UpdateMacro();
        void Collision();

        std::vector<std::vector<T> > rho, u, v;

protected:
        T omega;
    };


    template<class T>
    LBMNS<T>::LBMNS(int _nx, int _ny, T _viscosity) : LBMBASE<T>(_nx, _ny) {
        this->omega = 1.0/(3.0*_viscosity*this->dt/pow(this->dx, 2.0) + 0.5);

        this->rho = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->u = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->v = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
    }


    //----------------------
    //要修正
    //----------------------
    template<class T>
    void LBMNS<T>::Inlet(T _u, T _v) {
        for (int j = 0; j < this->ny; j++) {
            T rho0 = (this->f0t[0][j] + this->f2t[0][j] + this->f4t[0][j] + 2.0*(this->f3t[0][j] + this->f6t[0][j] + this->f7t[0][j]))/(1.0 - _u);
            this->f1t[0][j] = this->f3t[0][j] + 2.0*rho0*_u/3.0;
            this->f5t[0][j] = this->f7t[0][j] - 0.5*(this->f2t[0][j] - this->f4t[0][j]) + rho0*_u/6.0 + rho0*_v/2.0;
            this->f8t[0][j] = this->f6t[0][j] + 0.5*(this->f2t[0][j] - this->f4t[0][j]) + rho0*_u/6.0 - rho0*_v/2.0;
        }
    }


    template<class T>
    void LBMNS<T>::UpdateMacro() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->rho[i][j] = this->f0t[i][j] + this->f1t[i][j] + this->f2t[i][j] + this->f3t[i][j] + this->f4t[i][j] + this->f5t[i][j] + this->f6t[i][j] + this->f7t[i][j] + this->f8t[i][j];
                this->u[i][j] = (this->f1t[i][j] - this->f3t[i][j] + this->f5t[i][j] - this->f6t[i][j] - this->f7t[i][j] + this->f8t[i][j])/this->rho[i][j];
                this->v[i][j] = (this->f2t[i][j] - this->f4t[i][j] + this->f5t[i][j] + this->f6t[i][j] - this->f7t[i][j] - this->f8t[i][j])/this->rho[i][j];
            }
        }
    }


    template<class T>
    void LBMNS<T>::Collision() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                T u2 = pow(this->u[i][j], 2.0);
                T v2 = pow(this->v[i][j], 2.0);
                T u2v2 = u2 + v2;
                T uv = this->u[i][j]*this->v[i][j];
                T omu215 = 1.0 - 1.5*u2v2;

                this->f0tp1[i][j] = (1.0 - this->omega)*this->f0t[i][j] + this->omega*this->t0*this->rho[i][j]*omu215;
                this->f1tp1[i][j] = (1.0 - this->omega)*this->f1t[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 + 3.0*this->u[i][j] + 4.5*u2);
                this->f2tp1[i][j] = (1.0 - this->omega)*this->f2t[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 + 3.0*this->v[i][j] + 4.5*v2);
                this->f3tp1[i][j] = (1.0 - this->omega)*this->f3t[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 - 3.0*this->u[i][j] + 4.5*u2);
                this->f4tp1[i][j] = (1.0 - this->omega)*this->f4t[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 - 3.0*this->v[i][j] + 4.5*v2);
                this->f5tp1[i][j] = (1.0 - this->omega)*this->f5t[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 + 3.0*(this->u[i][j] + this->v[i][j]) + 4.5*(u2v2 + 2.0*uv));
                this->f6tp1[i][j] = (1.0 - this->omega)*this->f6t[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 - 3.0*(this->u[i][j] - this->v[i][j]) + 4.5*(u2v2 - 2.0*uv));
                this->f7tp1[i][j] = (1.0 - this->omega)*this->f7t[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 - 3.0*(this->u[i][j] + this->v[i][j]) + 4.5*(u2v2 + 2.0*uv));
                this->f8tp1[i][j] = (1.0 - this->omega)*this->f8t[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 + 3.0*(this->u[i][j] - this->v[i][j]) + 4.5*(u2v2 - 2.0*uv));
            }
        }
    }
}