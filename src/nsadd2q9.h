//*****************************************************************************
//  Title       :   src/nsadd2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/05
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "d2q9.h"


namespace PANSLBM2 {
    template<class T>
    class NSADd2q9 {
public:
        NSADd2q9() = delete;
        NSADd2q9(int _nx, int _ny, T _viscosity);
        NSADd2q9(const NSADd2q9<T>& _p);
        virtual ~NSADd2q9();

        template<class F, class G>
        void SetBarrier(F _f, G _g);
        template<class F, class G>
        void SetBoundary(F _f, G _g);
        virtual void Inlet(T _temperature0, T _temperature1);
        virtual void Stream();
        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        bool GetBarrier(int _i, int _j) const;
        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _i, int _j) const;
        virtual T GetV(int _i, int _j) const;
        virtual T GetTemperature(int _i, int _j) const;

        const int nx, ny;

protected:
        T dx, dt, t0, t1, t2, omega;
        D2Q9<T> f, g;
        T *rho, *u, *v, *temperature;
    };


    template<class T>
    NSADd2q9<T>::NSADd2q9(int _nx, int _ny, T _viscosity) : nx(_nx), ny(_ny), f(_nx, _ny), g(_nx, _ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;     this->dt = 1.0;
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;
        this->omega = 1.0/(3.0*_viscosity*this->dt/pow(this->dx, 2.0) + 0.5);

        this->rho = new T[this->nx*this->ny];
        this->u = new T[this->nx*this->ny];
        this->v = new T[this->nx*this->ny];
        this->temperature = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            this->u[i] = T();
            this->v[i] = T();
            this->temperature[i] = T();
        }
    }


    template<class T>
    NSADd2q9<T>::NSADd2q9(const NSADd2q9<T>& _p) : nx(_p.nx), ny(_p.ny), f(_p.f), g(_p.g) {
        this->dx = 1.0;     this->dt = 1.0;
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;
        this->omega = _p.omega;

        this->rho = new T[this->nx*this->ny];
        this->u = new T[this->nx*this->ny];
        this->v = new T[this->nx*this->ny];
        this->temperature = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = _p.rho[i];
            this->u[i] = _p.u[i];
            this->v[i] = _p.v[i];
            this->temperature[i] = _p.temperature[i];
        }
    }


    template<class T>
    NSADd2q9<T>::~NSADd2q9() {
        delete[] this->rho;
        delete[] this->u;
        delete[] this->v;
        delete[] this->temperature;
    }


    template<class T>
    template<class F, class G>
    void NSADd2q9<T>::SetBarrier(F _f, G _g) {
        this->f.SetBarrier(_f);
        this->g.SetBarrier(_g);
    }


    template<class T>
    template<class F, class G>
    void NSADd2q9<T>::SetBoundary(F _f, G _g) {
        this->f.SetBoundary(_f);
        this->g.SetBoundary(_g);
    }


    template<class T>
    void NSADd2q9<T>::Inlet(T _temperature0, T _temperature1) {
        for (int i = 0; i < this->nx; i++) {
            T temperature0 = 6.0*(_temperature0 - this->g.f0t[this->ny*i] - this->g.f1t[this->ny*i] - this->g.f3t[this->ny*i] - this->g.f4t[this->ny*i] - this->g.f7t[this->ny*i] - this->g.f8t[this->ny*i]);
            this->g.f2t[this->ny*i] = temperature0/9.0;
            this->g.f5t[this->ny*i] = temperature0/36.0;
            this->g.f6t[this->ny*i] = temperature0/36.0;

            T temperature1 = 6.0*(_temperature1 - this->g.f0t[this->ny*(i + 1) - 1] - this->g.f1t[this->ny*(i + 1) - 1] - this->g.f2t[this->ny*(i + 1) - 1] - this->g.f3t[this->ny*(i + 1) - 1] - this->g.f5t[this->ny*(i + 1) - 1] - this->g.f6t[this->ny*(i + 1) - 1]);
            this->g.f4t[this->ny*(i + 1) - 1] = temperature1/9.0;
            this->g.f7t[this->ny*(i + 1) - 1] = temperature1/36.0;
            this->g.f8t[this->ny*(i + 1) - 1] = temperature1/36.0;
        }
    }


    template<class T>
    void NSADd2q9<T>::Stream() {
        this->f.Stream();
        this->g.Stream();
    }


    template<class T>
    void NSADd2q9<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = this->f.f0t[i] + this->f.f1t[i] + this->f.f2t[i] + this->f.f3t[i] + this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] + this->f.f7t[i] + this->f.f8t[i];
            this->u[i] = (this->f.f1t[i] - this->f.f3t[i] + this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i])/this->rho[i];
            this->v[i] = (this->f.f2t[i] - this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i])/this->rho[i];
            this->temperature[i] = this->g.f0t[i] + this->g.f1t[i] + this->g.f2t[i] + this->g.f3t[i] + this->g.f4t[i] + this->g.f5t[i] + this->g.f6t[i] + this->g.f7t[i] + this->g.f8t[i];
        }
    }


    template<class T>
    void NSADd2q9<T>::Collision() {
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

            this->g.f0tp1[i] = (1.0 - this->omega)*this->g.f0t[i] + this->omega*this->t0*this->temperature[i];
            this->g.f1tp1[i] = (1.0 - this->omega)*this->g.f1t[i] + this->omega*this->t1*this->temperature[i]*(1.0 + 3.0*this->u[i]);
            this->g.f2tp1[i] = (1.0 - this->omega)*this->g.f2t[i] + this->omega*this->t1*this->temperature[i]*(1.0 + 3.0*this->v[i]);
            this->g.f3tp1[i] = (1.0 - this->omega)*this->g.f3t[i] + this->omega*this->t1*this->temperature[i]*(1.0 - 3.0*this->u[i]);
            this->g.f4tp1[i] = (1.0 - this->omega)*this->g.f4t[i] + this->omega*this->t1*this->temperature[i]*(1.0 - 3.0*this->v[i]);
            this->g.f5tp1[i] = (1.0 - this->omega)*this->g.f5t[i] + this->omega*this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] + 3.0*this->v[i]);
            this->g.f6tp1[i] = (1.0 - this->omega)*this->g.f6t[i] + this->omega*this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] + 3.0*this->v[i]);
            this->g.f7tp1[i] = (1.0 - this->omega)*this->g.f7t[i] + this->omega*this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] - 3.0*this->v[i]);
            this->g.f8tp1[i] = (1.0 - this->omega)*this->g.f8t[i] + this->omega*this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] - 3.0*this->v[i]);
        }
    }


    template<class T>
    void NSADd2q9<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T temperature0 = this->g.f0t[i] + this->g.f1t[i] + this->g.f2t[i] + this->g.f3t[i] + this->g.f4t[i] + this->g.f5t[i] + this->g.f6t[i] + this->g.f7t[i] + this->g.f8t[i];
            
            T rhog = -1.0e-4*(temperature0 - 1.5);

            T dxt0alpha = 3.0*this->dx*this->t0;
            T dxt1alpha = 3.0*this->dx*this->t1;
            T dxt2alpha = 3.0*this->dx*this->t2;
    
            this->f.f2t[i] += dxt1alpha*rhog;
            this->f.f4t[i] += dxt1alpha*-rhog;
            this->f.f5t[i] += dxt2alpha*rhog;
            this->f.f6t[i] += dxt2alpha*rhog;
            this->f.f7t[i] += dxt2alpha*-rhog;
            this->f.f8t[i] += dxt2alpha*-rhog;
        }
    }


    template<class T>
    T NSADd2q9<T>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->ny*_i + _j];
    }


    template<class T>
    T NSADd2q9<T>::GetU(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[this->ny*_i + _j];
    }


    template<class T>
    T NSADd2q9<T>::GetV(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->v[this->ny*_i + _j];
    }


    template<class T>
    T NSADd2q9<T>::GetTemperature(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->temperature[this->ny*_i + _j];
    }
}