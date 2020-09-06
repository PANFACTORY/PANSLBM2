//*****************************************************************************
//  Title       :   src/nsd2q9.h
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
        NSd2q9(int _nx, int _ny, T _viscosity);
        NSd2q9(const NSd2q9<T>& _p);
        virtual ~NSd2q9();

        template<class F>
        void SetBarrier(F _f);
        template<class F>
        void SetBoundary(F _f);
        virtual void Inlet(T _u, T _v);
        virtual void Stream();
        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        bool GetBarrier(int _i, int _j) const;
        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _i, int _j) const;
        virtual T GetV(int _i, int _j) const;
        
        const int nx, ny;

protected:
        T dx, dt, t0, t1, t2, omega;
        D2Q9<T> f;
        T *rho, *u, *v;
    };


    template<class T>
    NSd2q9<T>::NSd2q9(int _nx, int _ny, T _viscosity) : nx(_nx), ny(_ny), f(_nx, _ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;     this->dt = 1.0;
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;
        this->omega = 1.0/(3.0*_viscosity*this->dt/pow(this->dx, 2.0) + 0.5);

        this->rho = new T[this->nx*this->ny];
        this->u = new T[this->nx*this->ny];
        this->v = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            this->u[i] = T();
            this->v[i] = T();
        }
    }


    template<class T>
    NSd2q9<T>::NSd2q9(const NSd2q9<T>& _p) : nx(_p.nx), ny(_p.ny), f(_p.f) {
        this->dx = 1.0;     this->dt = 1.0;
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;
        this->omega = _p.omega;

        this->rho = new T[this->nx*this->ny];
        this->u = new T[this->nx*this->ny];
        this->v = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = _p.rho[i];
            this->u[i] = _p.u[i];
            this->v[i] = _p.v[i];
        }
    }


    template<class T>
    NSd2q9<T>::~NSd2q9() {
        delete[] this->rho;
        delete[] this->u;
        delete[] this->v;
    }


    template<class T>
    template<class F>
    void NSd2q9<T>::SetBarrier(F _f) {
        this->f.SetBarrier(_f);
    }


    template<class T>
    template<class F>
    void NSd2q9<T>::SetBoundary(F _f) {
        this->f.SetBoundary(_f);
    }


    template<class T>
    void NSd2q9<T>::Inlet(T _u, T _v) {
        for (int j = 0; j < this->ny; j++) {
            T rho0 = (this->f.f0t[j] + this->f.f2t[j] + this->f.f4t[j] + 2.0*(this->f.f3t[j] + this->f.f6t[j] + this->f.f7t[j]))/(1.0 - _u);
            this->f.f1t[j] = this->f.f3t[j] + 2.0*rho0*_u/3.0;
            this->f.f5t[j] = this->f.f7t[j] - 0.5*(this->f.f2t[j] - this->f.f4t[j]) + rho0*_u/6.0 + rho0*_v/2.0;
            this->f.f8t[j] = this->f.f6t[j] + 0.5*(this->f.f2t[j] - this->f.f4t[j]) + rho0*_u/6.0 - rho0*_v/2.0;
        }
    }


    template<class T>
    void NSd2q9<T>::Stream() {
        this->f.Stream();
    }


    template<class T>
    void NSd2q9<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = this->f.f0t[i] + this->f.f1t[i] + this->f.f2t[i] + this->f.f3t[i] + this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] + this->f.f7t[i] + this->f.f8t[i];
            this->u[i] = (this->f.f1t[i] - this->f.f3t[i] + this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i])/this->rho[i];
            this->v[i] = (this->f.f2t[i] - this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i])/this->rho[i];
        }
    }


    template<class T>
    void NSd2q9<T>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T u2 = pow(this->u[i], 2.0);
            T v2 = pow(this->v[i], 2.0);
            T u2v2 = u2 + v2;
            T uv = this->u[i]*this->v[i];
            T omu215 = 1.0 - 1.5*u2v2;

            this->f.f0tp1[i] = (1.0 - this->omega)*this->f.f0t[i] + this->omega*this->t0*this->rho[i]*omu215;
            this->f.f1tp1[i] = (1.0 - this->omega)*this->f.f1t[i] + this->omega*this->t1*this->rho[i]*(omu215 + 3.0*this->u[i] + 4.5*u2);
            this->f.f2tp1[i] = (1.0 - this->omega)*this->f.f2t[i] + this->omega*this->t1*this->rho[i]*(omu215 + 3.0*this->v[i] + 4.5*v2);
            this->f.f3tp1[i] = (1.0 - this->omega)*this->f.f3t[i] + this->omega*this->t1*this->rho[i]*(omu215 - 3.0*this->u[i] + 4.5*u2);
            this->f.f4tp1[i] = (1.0 - this->omega)*this->f.f4t[i] + this->omega*this->t1*this->rho[i]*(omu215 - 3.0*this->v[i] + 4.5*v2);
            this->f.f5tp1[i] = (1.0 - this->omega)*this->f.f5t[i] + this->omega*this->t2*this->rho[i]*(omu215 + 3.0*(this->u[i] + this->v[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f.f6tp1[i] = (1.0 - this->omega)*this->f.f6t[i] + this->omega*this->t2*this->rho[i]*(omu215 - 3.0*(this->u[i] - this->v[i]) + 4.5*(u2v2 - 2.0*uv));
            this->f.f7tp1[i] = (1.0 - this->omega)*this->f.f7t[i] + this->omega*this->t2*this->rho[i]*(omu215 - 3.0*(this->u[i] + this->v[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f.f8tp1[i] = (1.0 - this->omega)*this->f.f8t[i] + this->omega*this->t2*this->rho[i]*(omu215 + 3.0*(this->u[i] - this->v[i]) + 4.5*(u2v2 - 2.0*uv));
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
    T NSd2q9<T>::GetU(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[this->ny*_i + _j];
    }


    template<class T>
    T NSd2q9<T>::GetV(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->v[this->ny*_i + _j];
    }
}