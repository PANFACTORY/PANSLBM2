//*****************************************************************************
//  Title       :   src/nsadjointd2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/06
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "d2q9.h"


namespace PANSLBM2 {
    template<class T>
    class NSAdjointd2q9 {
public:
        NSAdjointd2q9() = delete;
        NSAdjointd2q9(int _nx, int _ny, int _nt, T _viscosity);
        NSAdjointd2q9(const NSAdjointd2q9<T>& _p);
        virtual ~NSAdjointd2q9();

        template<class F>
        void SetBarrier(F _f);
        template<class F>
        void SetBoundary(F _f);
        template<class F>
        void SetAlpha(F _f);

        virtual void Inlet(T _u, T _v);
        virtual void Stream();
        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        void SwitchDirection();

        bool GetBarrier(int _i, int _j) const;
        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _i, int _j) const;
        virtual T GetV(int _i, int _j) const;
        virtual T GetSensitivity(int _i, int _j) const;
        
        const int nx, ny, nt;

        int t = 0;

protected:
        T dx, dt, t0, t1, t2, omega;
        bool direction;                         //  true : direct analyze, false : inverse analyze
        D2Q9<T> f;
        T *rho, *u, *v, *alpha, *sensitivity;
    };


    template<class T>
    NSAdjointd2q9<T>::NSAdjointd2q9(int _nx, int _ny, int _nt, T _viscosity) : nx(_nx), ny(_ny), nt(_nt), f(_nx, _ny) {
        assert(0 < _nx && 0 < _ny && 0 < _nt);
        this->t = 0;
        this->dx = 1.0;     this->dt = 1.0;
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;
        this->omega = 1.0/(3.0*_viscosity*this->dt/pow(this->dx, 2.0) + 0.5);
        this->direction = true;

        this->rho = new T[this->nx*this->ny*this->nt];
        this->u = new T[this->nx*this->ny*this->nt];
        this->v = new T[this->nx*this->ny*this->nt];
        
        for (int i = 0; i < this->nx*this->ny*this->nt; i++) {
            this->rho[i] = T();
            this->u[i] = T();
            this->v[i] = T();
        }

        this->alpha = new T[this->nx*this->ny];
        this->sensitivity = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->alpha[i] = T();
            this->sensitivity[i] = T();
        }
    }


    template<class T>
    NSAdjointd2q9<T>::NSAdjointd2q9(const NSAdjointd2q9<T>& _p) : nx(_p.nx), ny(_p.ny), nt(_p.nt), f(_p.f) {
        this->t = 0;
        this->dx = _p.dx;   this->dt = _p.dt;
        this->t0 = _p.t0;   this->t1 = _p.t1;   this->t2 = _p.t2;
        this->omega = _p.omega;
        this->direction = true;

        this->rho = new T[this->nx*this->ny];
        this->u = new T[this->nx*this->ny];
        this->v = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny*this->nt; i++) {
            this->rho[i] = _p.rho[i];
            this->u[i] = _p.u[i];
            this->v[i] = _p.v[i];
        }

        this->alpha = new T[this->nx*this->ny];
        this->sensitivity = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->alpha[i] = _p.alpha[i];
            this->sensitivity[i] = _p.sensitivity[i];
        }
    }


    template<class T>
    NSAdjointd2q9<T>::~NSAdjointd2q9() {
        delete[] this->rho;
        delete[] this->u;
        delete[] this->v;
        delete[] this->sensitivity;
    }


    template<class T>
    template<class F>
    void NSAdjointd2q9<T>::SetBarrier(F _f) {
        this->f.SetBarrier(_f);
    }


    template<class T>
    template<class F>
    void NSAdjointd2q9<T>::SetBoundary(F _f) {
        this->f.SetBoundary(_f);
    }


    template<class T>
    template<class F>
    void NSAdjointd2q9<T>::SetAlpha(F _f) {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->alpha[this->ny*i + j] = _f(i, j);
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::Inlet(T _u, T _v) {
        if (this->direction) {
            for (int j = 0; j < this->ny; j++) {
                T rho0 = (this->f.f0t[j] + this->f.f2t[j] + this->f.f4t[j] + 2.0*(this->f.f3t[j] + this->f.f6t[j] + this->f.f7t[j]))/(1.0 - _u);
                this->f.f1t[j] = this->f.f3t[j] + 2.0*rho0*_u/3.0;
                this->f.f5t[j] = this->f.f7t[j] - 0.5*(this->f.f2t[j] - this->f.f4t[j]) + rho0*_u/6.0 + rho0*_v/2.0;
                this->f.f8t[j] = this->f.f6t[j] + 0.5*(this->f.f2t[j] - this->f.f4t[j]) + rho0*_u/6.0 - rho0*_v/2.0;
            }
        } else {
            for (int j = 0; j < this->ny; j++) {
                T term2 = (4.0*_u*this->f.f3t[j] + (_u + 3.0*_v)*this->f.f7t[j] + (_u - 3.0*_v)*this->f.f6t[j])/(3.0*(1.0 - _u));
                this->f.f1t[j] = this->f.f3t[j] + term2;
                this->f.f8t[j] = this->f.f6t[j] + term2;
                this->f.f5t[j] = this->f.f7t[j] + term2;
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::Stream() {
        this->f.Stream();
    }


    template<class T>
    void NSAdjointd2q9<T>::UpdateMacro() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                this->rho[ii] = this->f.f0t[i] + this->f.f1t[i] + this->f.f2t[i] + this->f.f3t[i] + this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] + this->f.f7t[i] + this->f.f8t[i];
                this->u[ii] = (this->f.f1t[i] - this->f.f3t[i] + this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i])/this->rho[ii];
                this->v[ii] = (this->f.f2t[i] - this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i])/this->rho[ii];
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                T sensx = this->t1*(this->f.f1t[i] - this->f.f3t[i]) + this->t2*(this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i]);
                T sensy = this->t1*(this->f.f2t[i] - this->f.f4t[i]) + this->t2*(this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i]);

                this->sensitivity[i] += 3.0*this->dx*(sensx*this->u[ii] + sensy*this->v[ii]) + this->u[ii];
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::Collision() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                T u2 = pow(this->u[ii], 2.0);
                T v2 = pow(this->v[ii], 2.0);
                T u2v2 = u2 + v2;
                T uv = this->u[ii]*this->v[ii];
                T omu215 = 1.0 - 1.5*u2v2;

                this->f.f0tp1[i] = (1.0 - this->omega)*this->f.f0t[i] + this->omega*this->t0*this->rho[ii]*omu215;
                this->f.f1tp1[i] = (1.0 - this->omega)*this->f.f1t[i] + this->omega*this->t1*this->rho[ii]*(omu215 + 3.0*this->u[ii] + 4.5*u2);
                this->f.f2tp1[i] = (1.0 - this->omega)*this->f.f2t[i] + this->omega*this->t1*this->rho[ii]*(omu215 + 3.0*this->v[ii] + 4.5*v2);
                this->f.f3tp1[i] = (1.0 - this->omega)*this->f.f3t[i] + this->omega*this->t1*this->rho[ii]*(omu215 - 3.0*this->u[ii] + 4.5*u2);
                this->f.f4tp1[i] = (1.0 - this->omega)*this->f.f4t[i] + this->omega*this->t1*this->rho[ii]*(omu215 - 3.0*this->v[ii] + 4.5*v2);
                this->f.f5tp1[i] = (1.0 - this->omega)*this->f.f5t[i] + this->omega*this->t2*this->rho[ii]*(omu215 + 3.0*(this->u[ii] + this->v[ii]) + 4.5*(u2v2 + 2.0*uv));
                this->f.f6tp1[i] = (1.0 - this->omega)*this->f.f6t[i] + this->omega*this->t2*this->rho[ii]*(omu215 - 3.0*(this->u[ii] - this->v[ii]) + 4.5*(u2v2 - 2.0*uv));
                this->f.f7tp1[i] = (1.0 - this->omega)*this->f.f7t[i] + this->omega*this->t2*this->rho[ii]*(omu215 - 3.0*(this->u[ii] + this->v[ii]) + 4.5*(u2v2 + 2.0*uv));
                this->f.f8tp1[i] = (1.0 - this->omega)*this->f.f8t[i] + this->omega*this->t2*this->rho[ii]*(omu215 + 3.0*(this->u[ii] - this->v[ii]) + 4.5*(u2v2 - 2.0*uv));
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;
                T rhoii = this->rho[ii];
                T uii = this->u[ii];
                T vii = this->v[ii];

                T u2 = pow(uii, 2.0);
                T v2 = pow(vii, 2.0);
                T u2v2 = u2 + v2;
                T uv = uii*vii;
                T omu215 = 1.0 - 1.5*u2v2;

                T f0eq = this->t0*rhoii*omu215;
                T f1eq = this->t1*rhoii*(omu215 - 3.0*uii + 4.5*u2);
                T f2eq = this->t1*rhoii*(omu215 - 3.0*vii + 4.5*v2);
                T f3eq = this->t1*rhoii*(omu215 + 3.0*uii + 4.5*u2);
                T f4eq = this->t1*rhoii*(omu215 + 3.0*vii + 4.5*v2);
                T f5eq = this->t2*rhoii*(omu215 - 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
                T f6eq = this->t2*rhoii*(omu215 + 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));
                T f7eq = this->t2*rhoii*(omu215 + 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
                T f8eq = this->t2*rhoii*(omu215 - 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));
                
                T faseq1 = (f0eq*this->f.f0t[i] + f1eq*this->f.f1t[i] + f2eq*this->f.f2t[i] + f3eq*this->f.f3t[i] + f4eq*this->f.f4t[i] + f5eq*this->f.f5t[i] + f6eq*this->f.f6t[i] + f7eq*this->f.f7t[i] + f8eq*this->f.f8t[i])/this->rho[ii];
                T w0fas0 = this->t0*this->f.f0t[i];
                T w1fas1 = this->t1*this->f.f1t[i];
                T w1fas2 = this->t1*this->f.f2t[i];
                T w1fas3 = this->t1*this->f.f3t[i];
                T w1fas4 = this->t1*this->f.f4t[i];
                T w2fas5 = this->t2*this->f.f5t[i];
                T w2fas6 = this->t2*this->f.f6t[i];
                T w2fas7 = this->t2*this->f.f7t[i];
                T w2fas8 = this->t2*this->f.f8t[i];
                T faseq2x = 3.0*(-w0fas0*uii + w1fas3*(1.0 + 2.0*uii) - w1fas4*uii + w1fas1*(-1.0 + 2.0*uii) - w1fas2*uii + w2fas7*(1.0 + 2.0*uii + 3.0*vii) + w2fas8*(-1.0 + 2.0*uii - 3.0*vii) + w2fas5*(-1.0 + 2.0*uii + 3.0*vii) + w2fas6*(1.0 + 2.0*uii - 3.0*vii));
                T faseq2y = 3.0*(-w0fas0*vii - w1fas3*vii + w1fas4*(1.0 + 2.0*vii) - w1fas1*vii + w1fas2*(-1.0 + 2.0*vii) + w2fas7*(1.0 + 3.0*uii + 2.0*vii) + w2fas8*(1.0 - 3.0*uii + 2.0*vii) + w2fas5*(-1.0 + 3.0*uii + 2.0*vii) + w2fas6*(-1.0 - 3.0*uii + 2.0*vii));

                this->f.f0tp1[i] = (1.0 - this->omega)*this->f.f0t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(-this->v[ii]));
                this->f.f1tp1[i] = (1.0 - this->omega)*this->f.f1t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(-this->v[ii]));
                this->f.f2tp1[i] = (1.0 - this->omega)*this->f.f2t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
                this->f.f3tp1[i] = (1.0 - this->omega)*this->f.f3t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(-this->v[ii]));
                this->f.f4tp1[i] = (1.0 - this->omega)*this->f.f4t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
                this->f.f5tp1[i] = (1.0 - this->omega)*this->f.f5t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
                this->f.f6tp1[i] = (1.0 - this->omega)*this->f.f6t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
                this->f.f7tp1[i] = (1.0 - this->omega)*this->f.f7t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
                this->f.f8tp1[i] = (1.0 - this->omega)*this->f.f8t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::ExternalForce() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                T tmprho = this->f.f0t[i] + this->f.f1t[i] + this->f.f2t[i] + this->f.f3t[i] + this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] + this->f.f7t[i] + this->f.f8t[i];
                T tmpu = (this->f.f1t[i] - this->f.f3t[i] + this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i])/tmprho;
                T tmpv = (this->f.f2t[i] - this->f.f4t[i] + this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i])/tmprho;

                T dxt1alpha = -3.0*this->dx*this->t1*this->alpha[i];
                T dxt2alpha = -3.0*this->dx*this->t2*this->alpha[i];

                this->f.f1t[i] += dxt1alpha*tmpu;
                this->f.f2t[i] += dxt1alpha*tmpv;
                this->f.f3t[i] += dxt1alpha*(-tmpu);
                this->f.f4t[i] += dxt1alpha*(-tmpv);
                this->f.f5t[i] += dxt2alpha*(tmpu + tmpv);
                this->f.f6t[i] += dxt2alpha*(-tmpu + tmpv);
                this->f.f7t[i] += dxt2alpha*(-tmpu - tmpv);
                this->f.f8t[i] += dxt2alpha*(tmpu - tmpv);
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                T ef1x = 3.0*this->dx*this->alpha[i]*(this->t1*(this->f.f1t[i] - this->f.f3t[i]) + this->t2*(this->f.f5t[i] - this->f.f6t[i] - this->f.f7t[i] + this->f.f8t[i]));
                T ef1y = 3.0*this->dx*this->alpha[i]*(this->t1*(this->f.f2t[i] - this->f.f4t[i]) + this->t2*(this->f.f5t[i] + this->f.f6t[i] - this->f.f7t[i] - this->f.f8t[i]));
            
                T ubyrho = this->u[ii]/this->rho[ii];
                T vbyrho = this->v[ii]/this->rho[ii];
                T onebyrho = 1.0/this->rho[ii];

                this->f.f0t[i] += ef1x*(-ubyrho) + ef1y*(-vbyrho) + this->alpha[i]*(-ubyrho);
                this->f.f1t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-vbyrho) + this->alpha[i]*(-onebyrho - ubyrho);
                this->f.f2t[i] += ef1x*(-ubyrho) + ef1y*(-onebyrho - vbyrho) + this->alpha[i]*(-ubyrho);
                this->f.f3t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(-vbyrho) + this->alpha[i]*(onebyrho - ubyrho);
                this->f.f4t[i] += ef1x*(-ubyrho) + ef1y*(onebyrho - vbyrho) + this->alpha[i]*(-ubyrho);
                this->f.f5t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->alpha[i]*(-onebyrho - ubyrho);
                this->f.f6t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->alpha[i]*(onebyrho - ubyrho);
                this->f.f7t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->alpha[i]*(onebyrho - ubyrho);
                this->f.f8t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->alpha[i]*(-onebyrho - ubyrho);
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::SwitchDirection() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                this->f.f0t[i] = T();   this->f.f0tp1[i] = T();
                this->f.f1t[i] = T();   this->f.f1tp1[i] = T();
                this->f.f2t[i] = T();   this->f.f2tp1[i] = T();
                this->f.f3t[i] = T();   this->f.f3tp1[i] = T();
                this->f.f4t[i] = T();   this->f.f4tp1[i] = T();
                this->f.f5t[i] = T();   this->f.f5tp1[i] = T();
                this->f.f6t[i] = T();   this->f.f6tp1[i] = T();
                this->f.f7t[i] = T();   this->f.f7tp1[i] = T();
                this->f.f8t[i] = T();   this->f.f8tp1[i] = T();
            }
        }
        this->direction = !this->direction;
    }


    template<class T>
    bool NSAdjointd2q9<T>::GetBarrier(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->f.GetBarrier(_i, _j);
    }


    template<class T>
    T NSAdjointd2q9<T>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->nx*this->ny*(this->nt - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetU(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[this->nx*this->ny*(this->nt - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetV(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->v[this->nx*this->ny*(this->nt - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetSensitivity(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->sensitivity[this->ny*_i + _j];
    }
}