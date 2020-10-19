//*****************************************************************************
//  Title       :   src/ns.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/06
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "d2q9.h"


namespace PANSLBM2 {
    template<class T>
    class NSd2q9 {
public:
        NSd2q9() = delete;
        NSd2q9(D2Q9<T>* _f, T _viscosity);
        NSd2q9(const NSd2q9<T>& _e);
        virtual ~NSd2q9();

        virtual void Inlet(T _u, T _v);
        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        virtual T GetRho(int _i, int _j) const;
        virtual T GetUx(int _i, int _j) const;
        virtual T GetUy(int _i, int _j) const;
        
        const int nx, ny;

protected:
        T omega;
        D2Q9<T>* f;
        T *rho, *ux, *uy;
    };


    template<class T>
    NSd2q9<T>::NSd2q9(D2Q9<T>* _f, T _viscosity) : nx(_f->nx), ny(_f->ny) {
        assert(0 < _f->nx && 0 < _f->ny);
        this->f = _f;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/pow(this->f->dx, 2.0) + 0.5);

        this->rho = new T[this->nx*this->ny];
        this->ux = new T[this->nx*this->ny];
        this->uy = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            this->ux[i] = T();
            this->uy[i] = T();
        }
    }


    template<class T>
    NSd2q9<T>::NSd2q9(const NSd2q9<T>& _e) : nx(_e.nx), ny(_e.ny), f(_e.f) {
        this->f = _e.f;
        this->omega = _e.omega;

        this->rho = new T[this->nx*this->ny];
        this->ux = new T[this->nx*this->ny];
        this->uy = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = _e.rho[i];
            this->ux[i] = _e.ux[i];
            this->uy[i] = _e.uy[i];
        }
    }


    template<class T>
    NSd2q9<T>::~NSd2q9() {
        delete[] this->rho;
        delete[] this->ux;
        delete[] this->uy;
    }


    template<class T>
    void NSd2q9<T>::Inlet(T _ux, T _uy) {
        for (int j = 0; j < this->ny; j++) {
            T rho0 = (this->f->f0t[j] + this->f->f2t[j] + this->f->f4t[j] + 2.0*(this->f->f3t[j] + this->f->f6t[j] + this->f->f7t[j]))/(1.0 - _ux);
            this->f->f1t[j] = this->f->f3t[j] + 2.0*rho0*_ux/3.0;
            this->f->f5t[j] = this->f->f7t[j] - 0.5*(this->f->f2t[j] - this->f->f4t[j]) + rho0*_ux/6.0 + rho0*_uy/2.0;
            this->f->f8t[j] = this->f->f6t[j] + 0.5*(this->f->f2t[j] - this->f->f4t[j]) + rho0*_ux/6.0 - rho0*_uy/2.0;
        }
    }


    template<class T>
    void NSd2q9<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = this->f->f0t[i] + this->f->f1t[i] + this->f->f2t[i] + this->f->f3t[i] + this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] + this->f->f7t[i] + this->f->f8t[i];
            this->ux[i] = (this->f->f1t[i] - this->f->f3t[i] + this->f->f5t[i] - this->f->f6t[i] - this->f->f7t[i] + this->f->f8t[i])/this->rho[i];
            this->uy[i] = (this->f->f2t[i] - this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] - this->f->f7t[i] - this->f->f8t[i])/this->rho[i];
        }
    }


    template<class T>
    void NSd2q9<T>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T u2 = pow(this->ux[i], 2.0);
            T v2 = pow(this->uy[i], 2.0);
            T u2v2 = u2 + v2;
            T uv = this->ux[i]*this->uy[i];
            T omu215 = 1.0 - 1.5*u2v2;

            this->f->f0tp1[i] = (1.0 - this->omega)*this->f->f0t[i] + this->omega*this->f->t0*this->rho[i]*omu215;
            this->f->f1tp1[i] = (1.0 - this->omega)*this->f->f1t[i] + this->omega*this->f->t1*this->rho[i]*(omu215 + 3.0*this->ux[i] + 4.5*u2);
            this->f->f2tp1[i] = (1.0 - this->omega)*this->f->f2t[i] + this->omega*this->f->t1*this->rho[i]*(omu215 + 3.0*this->uy[i] + 4.5*v2);
            this->f->f3tp1[i] = (1.0 - this->omega)*this->f->f3t[i] + this->omega*this->f->t1*this->rho[i]*(omu215 - 3.0*this->ux[i] + 4.5*u2);
            this->f->f4tp1[i] = (1.0 - this->omega)*this->f->f4t[i] + this->omega*this->f->t1*this->rho[i]*(omu215 - 3.0*this->uy[i] + 4.5*v2);
            this->f->f5tp1[i] = (1.0 - this->omega)*this->f->f5t[i] + this->omega*this->f->t2*this->rho[i]*(omu215 + 3.0*(this->ux[i] + this->uy[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f->f6tp1[i] = (1.0 - this->omega)*this->f->f6t[i] + this->omega*this->f->t2*this->rho[i]*(omu215 - 3.0*(this->ux[i] - this->uy[i]) + 4.5*(u2v2 - 2.0*uv));
            this->f->f7tp1[i] = (1.0 - this->omega)*this->f->f7t[i] + this->omega*this->f->t2*this->rho[i]*(omu215 - 3.0*(this->ux[i] + this->uy[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f->f8tp1[i] = (1.0 - this->omega)*this->f->f8t[i] + this->omega*this->f->t2*this->rho[i]*(omu215 + 3.0*(this->ux[i] - this->uy[i]) + 4.5*(u2v2 - 2.0*uv));
        }
    }


    template<class T>
    void NSd2q9<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            
        }
    }


    template<class T>
    T NSd2q9<T>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->ny*_i + _j];
    }


    template<class T>
    T NSd2q9<T>::GetUx(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->ux[this->ny*_i + _j];
    }


    template<class T>
    T NSd2q9<T>::GetUy(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->uy[this->ny*_i + _j];
    }
}