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

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void SetRho(int _i, int _j, T _ux, T _uy);
        virtual void SetUxUy(int _i, int _j, T _ux, T _uy);
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
    void NSd2q9<T>::SetRho(int _i, int _j, T _rho, T _u) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T ux0 = 1.0 - (this->f->f0t[ij] + this->f->f2t[ij] + this->f->f4t[ij] + 2.0*(this->f->f3t[ij] + this->f->f6t[ij] + this->f->f7t[ij]))/_rho;
            this->f->f1t[ij] = this->f->f3t[ij] + 2.0*_rho*ux0/3.0;
            this->f->f5t[ij] = this->f->f7t[ij] - 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) + _rho*ux0/6.0 + _rho*_u/2.0;
            this->f->f8t[ij] = this->f->f6t[ij] + 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) + _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_i == this->nx - 1) {
            T ux0 = -1.0 + (this->f->f0t[ij] + this->f->f2t[ij] + this->f->f4t[ij] + 2.0*(this->f->f1t[ij] + this->f->f5t[ij] + this->f->f8t[ij]))/_rho;
            this->f->f3t[ij] = this->f->f1t[ij] - 2.0*_rho*ux0/3.0;
            this->f->f6t[ij] = this->f->f8t[ij] - 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) - _rho*ux0/6.0 + _rho*_u/2.0;
            this->f->f7t[ij] = this->f->f6t[ij] + 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) - _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_j == 0) {
            T uy0 = 1.0 - (this->f->f0t[ij] + this->f->f1t[ij] + this->f->f3t[ij] + 2.0*(this->f->f4t[ij] + this->f->f7t[ij] + this->f->f8t[ij]))/_rho;
            this->f->f2t[ij] = this->f->f4t[ij] + 2.0*_rho*uy0/3.0;
            this->f->f5t[ij] = this->f->f7t[ij] - 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) + _rho*_u/2.0 + _rho*uy0/6.0;
            this->f->f6t[ij] = this->f->f8t[ij] + 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) - _rho*_u/2.0 + _rho*uy0/6.0;
        } else if (_j == this->ny - 1) {
            T uy0 = -1.0 - (this->f->f0t[ij] + this->f->f1t[ij] + this->f->f3t[ij] + 2.0*(this->f->f2t[ij] + this->f->f5t[ij] + this->f->f6t[ij]))/_rho;
            this->f->f4t[ij] = this->f->f2t[ij] - 2.0*_rho*uy0/3.0;
            this->f->f7t[ij] = this->f->f5t[ij] + 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) - _rho*_u/2.0 - _rho*uy0/6.0;
            this->f->f8t[ij] = this->f->f6t[ij] - 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) + _rho*_u/2.0 - _rho*uy0/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void NSd2q9<T>::SetUxUy(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (this->f->f0t[ij] + this->f->f2t[ij] + this->f->f4t[ij] + 2.0*(this->f->f3t[ij] + this->f->f6t[ij] + this->f->f7t[ij]))/(1.0 - _ux);
            this->f->f1t[ij] = this->f->f3t[ij] + 2.0*rho0*_ux/3.0;
            this->f->f5t[ij] = this->f->f7t[ij] - 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) + rho0*_ux/6.0 + rho0*_uy/2.0;
            this->f->f8t[ij] = this->f->f6t[ij] + 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) + rho0*_ux/6.0 - rho0*_uy/2.0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (this->f->f0t[ij] + this->f->f2t[ij] + this->f->f4t[ij] + 2.0*(this->f->f1t[ij] + this->f->f5t[ij] + this->f->f8t[ij]))/(1.0 + _ux);
            this->f->f3t[ij] = this->f->f1t[ij] - 2.0*rho0*_ux/3.0;
            this->f->f7t[ij] = this->f->f5t[ij] + 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) - rho0*_ux/6.0 - rho0*_uy/2.0;
            this->f->f6t[ij] = this->f->f8t[ij] - 0.5*(this->f->f2t[ij] - this->f->f4t[ij]) - rho0*_ux/6.0 + rho0*_uy/2.0;
        } else if (_j == 0) {
            T rho0 = (this->f->f0t[ij] + this->f->f1t[ij] + this->f->f3t[ij] + 2.0*(this->f->f4t[ij] + this->f->f7t[ij] + this->f->f8t[ij]))/(1.0 - _uy);
            this->f->f2t[ij] = this->f->f4t[ij] + 2.0*rho0*_uy/3.0;
            this->f->f5t[ij] = this->f->f7t[ij] - 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) + rho0*_ux/2.0 + rho0*_uy/6.0;
            this->f->f6t[ij] = this->f->f8t[ij] + 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) - rho0*_ux/2.0 + rho0*_uy/6.0;
        } else if (_j == this->ny - 1) {
            T rho0 = (this->f->f0t[ij] + this->f->f1t[ij] + this->f->f3t[ij] + 2.0*(this->f->f2t[ij] + this->f->f5t[ij] + this->f->f6t[ij]))/(1.0 + _uy);
            this->f->f4t[ij] = this->f->f2t[ij] - 2.0*rho0*_uy/3.0;
            this->f->f7t[ij] = this->f->f5t[ij] + 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) - rho0*_ux/2.0 - rho0*_uy/6.0;
            this->f->f8t[ij] = this->f->f6t[ij] - 0.5*(this->f->f1t[ij] - this->f->f3t[ij]) + rho0*_ux/2.0 - rho0*_uy/6.0;
        } else {
            //  境界に沿っていないことを警告する
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