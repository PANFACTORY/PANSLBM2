//*****************************************************************************
//  Title       :   src/nsadjointd2q9_2.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/17
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "d2q9.h"


namespace PANSLBM2 {
    template<class T>
    class NSAdjointd2q9 {
public:
        NSAdjointd2q9() = delete;
        NSAdjointd2q9(D2Q9<T>* _f, int _nt, T _viscosity);
        NSAdjointd2q9(const NSAdjointd2q9<T>& _e);
        virtual ~NSAdjointd2q9();

        template<class F>
        void SetAlpha(F _f);

        virtual void Inlet(T _u, T _v, T _rho);
        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        bool CheckConvergence(T _eps);
        void SwitchDirection();

        virtual T GetRho(int _i, int _j) const;
        virtual T GetUx(int _i, int _j) const;
        virtual T GetUy(int _i, int _j) const;
        virtual T GetQ(int _i, int _j) const;
        virtual T GetVx(int _i, int _j) const;
        virtual T GetVy(int _i, int _j) const;
        virtual T GetSensitivity(int _i, int _j) const;
        
        const int nx, ny, nt;

        int t = 0;

protected:
        D2Q9<T>* f;
        int tmax;
        T omega;
        bool direction;                         //  true : direct analyze, false : inverse analyze
        T *rho, *ux, *uy, *q, *vx, *vy, *alpha, *sensitivity;
    };


    template<class T>
    NSAdjointd2q9<T>::NSAdjointd2q9(D2Q9<T>* _f, int _nt, T _viscosity) : nx(_f->nx), ny(_f->ny), nt(_nt) {
        assert(0 < _f->nx && 0 < _f->ny && 0 < _nt);
        
        this->t = 0;
        this->f = _f;
        this->tmax = 0;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/pow(this->f->dx, 2.0) + 0.5);
        this->direction = true;
        
        this->rho = new T[this->nx*this->ny*this->nt];
        this->ux = new T[this->nx*this->ny*this->nt];
        this->uy = new T[this->nx*this->ny*this->nt];
        
        for (int i = 0; i < this->nx*this->ny*this->nt; i++) {
            this->rho[i] = T();
            this->ux[i] = T();
            this->uy[i] = T();
        }

        this->q = new T[this->nx*this->ny];
        this->vx = new T[this->nx*this->ny];
        this->vy = new T[this->nx*this->ny];
        this->alpha = new T[this->nx*this->ny];
        this->sensitivity = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->q[i] = T();
            this->vx[i] = T();
            this->vy[i] = T();
            this->alpha[i] = T();
            this->sensitivity[i] = T();
        }
    }


    template<class T>
    NSAdjointd2q9<T>::NSAdjointd2q9(const NSAdjointd2q9<T>& _e) : nx(_e.nx), ny(_e.ny), nt(_e.nt) {
        this->t = _e.t;
        this->f = _e.f;
        this->tmax = _e.tmax;
        this->omega = _e.omega;
        this->direction = _e.direction;
        
        this->rho = new T[this->nx*this->ny];
        this->ux = new T[this->nx*this->ny];
        this->uy = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny*this->nt; i++) {
            this->rho[i] = _e.rho[i];
            this->ux[i] = _e.ux[i];
            this->uy[i] = _e.uy[i];
        }

        this->q = new T[this->nx*this->ny];
        this->vx = new T[this->nx*this->ny];
        this->vy = new T[this->nx*this->ny];
        this->alpha = new T[this->nx*this->ny];
        this->sensitivity = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->q[i] = _e.q[i];
            this->vx[i] = _e.vx[i];
            this->vy[i] = _e.vy[i];
            this->alpha[i] = _e.alpha[i];
            this->sensitivity[i] = _e.sensitivity[i];
        }
    }


    template<class T>
    NSAdjointd2q9<T>::~NSAdjointd2q9() {
        delete[] this->rho;
        delete[] this->ux;
        delete[] this->uy;
        delete[] this->q;
        delete[] this->vx;
        delete[] this->vy;
        delete[] this->alpha;
        delete[] this->sensitivity;
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
    void NSAdjointd2q9<T>::Inlet(T _ux, T _uy, T _rho) {
        if (this->direction) {
            //  Boundary condition for velocity
            for (int j = (int)(0.7*this->ny); j < (int)(0.9*this->ny); j++) {
                T ui = -_ux*(j - 0.7*this->ny)*(j - 0.9*this->ny)/pow(0.1*this->ny, 2.0);
                T rho0 = (this->f->f0t[j] + this->f->f2t[j] + this->f->f4t[j] + 2.0*(this->f->f3t[j] + this->f->f6t[j] + this->f->f7t[j]))/(1.0 - ui);
                this->f->f1t[j] = this->f->f3t[j] + 2.0*rho0*ui/3.0;
                this->f->f5t[j] = this->f->f7t[j] - 0.5*(this->f->f2t[j] - this->f->f4t[j]) + rho0*ui/6.0 + rho0*_uy/2.0;
                this->f->f8t[j] = this->f->f6t[j] + 0.5*(this->f->f2t[j] - this->f->f4t[j]) + rho0*ui/6.0 - rho0*_uy/2.0;
            }
            //  Boundary condition for pressure
            for (int i = (int)(0.7*this->nx); i < (int)(0.9*this->nx); i++) {
                T uy0 = 1.0 - (this->f->f0t[this->ny*i] + this->f->f1t[this->ny*i] + this->f->f3t[this->ny*i] + 2.0*(this->f->f4t[this->ny*i] + this->f->f7t[this->ny*i] + this->f->f8t[this->ny*i]))/_rho;
                this->f->f2t[this->ny*i] = this->f->f4t[this->ny*i] + 2.0*_rho*uy0/3.0;
                this->f->f5t[this->ny*i] = this->f->f7t[this->ny*i] - 0.5*(this->f->f1t[this->ny*i] - this->f->f3t[this->ny*i]) + _rho*uy0/6.0;
                this->f->f6t[this->ny*i] = this->f->f8t[this->ny*i] + 0.5*(this->f->f1t[this->ny*i] - this->f->f3t[this->ny*i]) + _rho*uy0/6.0;
            }
        } else {
            //  Boundary condition for velocity
            for (int j = (int)(0.7*this->ny); j < (int)(0.9*this->ny); j++) {
                T ui = -_ux*(j - 0.7*this->ny)*(j - 0.9*this->ny)/pow(0.1*this->ny, 2.0);
                T term2 = (-2.0 + 4.0*ui*this->f->f3t[j] + (ui + 3.0*_uy)*this->f->f7t[j] + (ui - 3.0*_uy)*this->f->f6t[j])/(3.0*(1.0 - ui));
                this->f->f1t[j] = this->f->f3t[j] + term2;
                this->f->f8t[j] = this->f->f6t[j] + term2;
                this->f->f5t[j] = this->f->f7t[j] + term2;
            }
            //  Boundary condition for pressure
            for (int i = (int)(0.7*this->nx); i < (int)(0.9*this->nx); i++) {
                T term2 = -(4.0*this->f->f4t[this->ny*i] + this->f->f7t[this->ny*i] + this->f->f8t[this->ny*i])/3.0;
                this->f->f2t[this->ny*i] = this->f->f4t[this->ny*i] + term2;
                this->f->f5t[this->ny*i] = this->f->f7t[this->ny*i] + term2;
                this->f->f6t[this->ny*i] = this->f->f8t[this->ny*i] + term2;
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::UpdateMacro() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                this->rho[ii] = this->f->f0t[i] + this->f->f1t[i] + this->f->f2t[i] + this->f->f3t[i] + this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] + this->f->f7t[i] + this->f->f8t[i];
                this->ux[ii] = (this->f->f1t[i] - this->f->f3t[i] + this->f->f5t[i] - this->f->f6t[i] - this->f->f7t[i] + this->f->f8t[i])/this->rho[ii];
                this->uy[ii] = (this->f->f2t[i] - this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] - this->f->f7t[i] - this->f->f8t[i])/this->rho[ii];
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;
                T rhoii = this->rho[ii];
                T uii = this->ux[ii];
                T vii = this->uy[ii];

                T u2 = pow(uii, 2.0);
                T v2 = pow(vii, 2.0);
                T u2v2 = u2 + v2;
                T uv = uii*vii;
                T omu215 = 1.0 - 1.5*u2v2;

                T f0eq = this->f->t0*omu215;
                T f1eq = this->f->t1*(omu215 - 3.0*uii + 4.5*u2);
                T f2eq = this->f->t1*(omu215 - 3.0*vii + 4.5*v2);
                T f3eq = this->f->t1*(omu215 + 3.0*uii + 4.5*u2);
                T f4eq = this->f->t1*(omu215 + 3.0*vii + 4.5*v2);
                T f5eq = this->f->t2*(omu215 - 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
                T f6eq = this->f->t2*(omu215 + 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));
                T f7eq = this->f->t2*(omu215 + 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
                T f8eq = this->f->t2*(omu215 - 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));
                this->q[i] = f0eq*this->f->f0t[i] + f1eq*this->f->f1t[i] + f2eq*this->f->f2t[i] + f3eq*this->f->f3t[i] + f4eq*this->f->f4t[i] + f5eq*this->f->f5t[i] + f6eq*this->f->f6t[i] + f7eq*this->f->f7t[i] + f8eq*this->f->f8t[i];
                
                T fas0 = this->f->t0*this->f->f0t[i];
                T fas1 = this->f->t1*this->f->f1t[i];
                T fas2 = this->f->t1*this->f->f2t[i];
                T fas3 = this->f->t1*this->f->f3t[i];
                T fas4 = this->f->t1*this->f->f4t[i];
                T fas5 = this->f->t2*this->f->f5t[i];
                T fas6 = this->f->t2*this->f->f6t[i];
                T fas7 = this->f->t2*this->f->f7t[i];
                T fas8 = this->f->t2*this->f->f8t[i];
                this->vx[i] = -fas0*uii + fas1*(-1.0 + 2.0*uii) - fas2*uii + fas3*(1.0 + 2.0*uii) - fas4*uii + fas5*(-1.0 + 2.0*uii + 3.0*vii) + fas6*(1.0 + 2.0*uii - 3.0*vii) + fas7*(1.0 + 2.0*uii + 3.0*vii) + fas8*(-1.0 + 2.0*uii - 3.0*vii);
                this->vy[i] = -fas0*vii - fas1*vii + fas2*(-1.0 + 2.0*vii) - fas3*vii + fas4*(1.0 + 2.0*vii) + fas5*(-1.0 + 3.0*uii + 2.0*vii) + fas6*(-1.0 - 3.0*uii + 2.0*vii) + fas7*(1.0 + 3.0*uii + 2.0*vii) + fas8*(1.0 - 3.0*uii + 2.0*vii);

                T sensx = this->f->t1*(-this->f->f1t[i] + this->f->f3t[i]) + this->f->t2*(-this->f->f5t[i] + this->f->f6t[i] + this->f->f7t[i] - this->f->f8t[i]);
                T sensy = this->f->t1*(-this->f->f2t[i] + this->f->f4t[i]) + this->f->t2*(-this->f->f5t[i] - this->f->f6t[i] + this->f->f7t[i] + this->f->f8t[i]);
                this->sensitivity[i] += 3.0*this->f->dx*(sensx*this->ux[ii] + sensy*this->uy[ii]);
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::Collision() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;

                T u2 = pow(this->ux[ii], 2.0);
                T v2 = pow(this->uy[ii], 2.0);
                T u2v2 = u2 + v2;
                T uv = this->ux[ii]*this->uy[ii];
                T omu215 = 1.0 - 1.5*u2v2;

                this->f->f0tp1[i] = (1.0 - this->omega)*this->f->f0t[i] + this->omega*this->f->t0*this->rho[ii]*omu215;
                this->f->f1tp1[i] = (1.0 - this->omega)*this->f->f1t[i] + this->omega*this->f->t1*this->rho[ii]*(omu215 + 3.0*this->ux[ii] + 4.5*u2);
                this->f->f2tp1[i] = (1.0 - this->omega)*this->f->f2t[i] + this->omega*this->f->t1*this->rho[ii]*(omu215 + 3.0*this->uy[ii] + 4.5*v2);
                this->f->f3tp1[i] = (1.0 - this->omega)*this->f->f3t[i] + this->omega*this->f->t1*this->rho[ii]*(omu215 - 3.0*this->ux[ii] + 4.5*u2);
                this->f->f4tp1[i] = (1.0 - this->omega)*this->f->f4t[i] + this->omega*this->f->t1*this->rho[ii]*(omu215 - 3.0*this->uy[ii] + 4.5*v2);
                this->f->f5tp1[i] = (1.0 - this->omega)*this->f->f5t[i] + this->omega*this->f->t2*this->rho[ii]*(omu215 + 3.0*(this->ux[ii] + this->uy[ii]) + 4.5*(u2v2 + 2.0*uv));
                this->f->f6tp1[i] = (1.0 - this->omega)*this->f->f6t[i] + this->omega*this->f->t2*this->rho[ii]*(omu215 - 3.0*(this->ux[ii] - this->uy[ii]) + 4.5*(u2v2 - 2.0*uv));
                this->f->f7tp1[i] = (1.0 - this->omega)*this->f->f7t[i] + this->omega*this->f->t2*this->rho[ii]*(omu215 - 3.0*(this->ux[ii] + this->uy[ii]) + 4.5*(u2v2 + 2.0*uv));
                this->f->f8tp1[i] = (1.0 - this->omega)*this->f->f8t[i] + this->omega*this->f->t2*this->rho[ii]*(omu215 + 3.0*(this->ux[ii] - this->uy[ii]) + 4.5*(u2v2 - 2.0*uv));
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                int ii = this->nx*this->ny*this->t + i;
                
                this->f->f0tp1[i] = (1.0 - this->omega)*this->f->f0t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-this->ux[ii]) + this->vy[i]*(-this->uy[ii])));
                this->f->f1tp1[i] = (1.0 - this->omega)*this->f->f1t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-1.0 - this->ux[ii]) + this->vy[i]*(-this->uy[ii])));
                this->f->f2tp1[i] = (1.0 - this->omega)*this->f->f2t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-this->ux[ii]) + this->vy[i]*(-1.0 - this->uy[ii])));
                this->f->f3tp1[i] = (1.0 - this->omega)*this->f->f3t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(1.0 - this->ux[ii]) + this->vy[i]*(-this->uy[ii])));
                this->f->f4tp1[i] = (1.0 - this->omega)*this->f->f4t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-this->ux[ii]) + this->vy[i]*(1.0 - this->uy[ii])));
                this->f->f5tp1[i] = (1.0 - this->omega)*this->f->f5t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-1.0 - this->ux[ii]) + this->vy[i]*(-1.0 - this->uy[ii])));
                this->f->f6tp1[i] = (1.0 - this->omega)*this->f->f6t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(1.0 - this->ux[ii]) + this->vy[i]*(-1.0 - this->uy[ii])));
                this->f->f7tp1[i] = (1.0 - this->omega)*this->f->f7t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(1.0 - this->ux[ii]) + this->vy[i]*(1.0 - this->uy[ii])));
                this->f->f8tp1[i] = (1.0 - this->omega)*this->f->f8t[i] + this->omega*(this->q[i] + 3.0*(this->vx[i]*(-1.0 - this->ux[ii]) + this->vy[i]*(1.0 - this->uy[ii])));
            }
        }
    }


    template<class T>
    void NSAdjointd2q9<T>::ExternalForce() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                T tmprho = this->f->f0t[i] + this->f->f1t[i] + this->f->f2t[i] + this->f->f3t[i] + this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] + this->f->f7t[i] + this->f->f8t[i];
                T tmpux = (this->f->f1t[i] - this->f->f3t[i] + this->f->f5t[i] - this->f->f6t[i] - this->f->f7t[i] + this->f->f8t[i])/tmprho;
                T tmpuy = (this->f->f2t[i] - this->f->f4t[i] + this->f->f5t[i] + this->f->f6t[i] - this->f->f7t[i] - this->f->f8t[i])/tmprho;

                T dxt1alpha = -3.0*this->f->dx*this->f->t1*this->alpha[i];
                T dxt2alpha = -3.0*this->f->dx*this->f->t2*this->alpha[i];

                this->f->f1t[i] += dxt1alpha*tmpux;
                this->f->f2t[i] += dxt1alpha*tmpuy;
                this->f->f3t[i] += dxt1alpha*(-tmpux);
                this->f->f4t[i] += dxt1alpha*(-tmpuy);
                this->f->f5t[i] += dxt2alpha*(tmpux + tmpuy);
                this->f->f6t[i] += dxt2alpha*(-tmpux + tmpuy);
                this->f->f7t[i] += dxt2alpha*(-tmpux - tmpuy);
                this->f->f8t[i] += dxt2alpha*(tmpux - tmpuy);
            }
        } else {
            for (int i = 0; i < this->nx*this->ny; i++) {
                T ef1x = -3.0*this->f->dx*this->alpha[i]*(this->f->t1*(this->f->f1t[i] - this->f->f3t[i]) + this->f->t2*(this->f->f5t[i] - this->f->f6t[i] - this->f->f7t[i] + this->f->f8t[i]));
                T ef1y = -3.0*this->f->dx*this->alpha[i]*(this->f->t1*(this->f->f2t[i] - this->f->f4t[i]) + this->f->t2*(this->f->f5t[i] + this->f->f6t[i] - this->f->f7t[i] - this->f->f8t[i]));
            
                int ii = this->nx*this->ny*this->t + i;

                this->f->f0t[i] -= (ef1x*(-this->ux[ii]) + ef1y*(-this->uy[ii]))/this->rho[ii];
                this->f->f1t[i] -= (ef1x*(-1.0 - this->ux[ii]) + ef1y*(-this->uy[ii]))/this->rho[ii];
                this->f->f2t[i] -= (ef1x*(-this->ux[ii]) + ef1y*(-1.0 - this->uy[ii]))/this->rho[ii];
                this->f->f3t[i] -= (ef1x*(1.0 - this->ux[ii]) + ef1y*(-this->uy[ii]))/this->rho[ii];
                this->f->f4t[i] -= (ef1x*(-this->ux[ii]) + ef1y*(1.0 - this->uy[ii]))/this->rho[ii];
                this->f->f5t[i] -= (ef1x*(-1.0 - this->ux[ii]) + ef1y*(-1.0 - this->uy[ii]))/this->rho[ii];
                this->f->f6t[i] -= (ef1x*(1.0 - this->ux[ii]) + ef1y*(-1.0 - this->uy[ii]))/this->rho[ii];
                this->f->f7t[i] -= (ef1x*(1.0 - this->ux[ii]) + ef1y*(1.0 - this->uy[ii]))/this->rho[ii];
                this->f->f8t[i] -= (ef1x*(-1.0 - this->ux[ii]) + ef1y*(1.0 - this->uy[ii]))/this->rho[ii];
            }
        }
    }


    template<class T>
    bool NSAdjointd2q9<T>::CheckConvergence(T _eps) {
        T unorm = T();
        T dunorm = T();
        for (int i = 0; i < this->nx*this->ny; i++) {
            int ii = this->nx*this->ny*this->t + i;
            unorm += pow(this->ux[ii], 2.0) + pow(this->uy[ii], 2.0);
            int jj = this->nx*this->ny*(this->t - 1) + i;
            dunorm += pow(this->ux[ii] - this->ux[jj], 2.0) + pow(this->uy[ii] - this->uy[jj], 2.0);
        }
        this->tmax = this->t;
        return sqrt(dunorm/unorm) < _eps ? true : false;
    }


    template<class T>
    void NSAdjointd2q9<T>::SwitchDirection() {
        if (this->direction) {
            for (int i = 0; i < this->nx*this->ny; i++) {
                this->f->f0t[i] = T();   this->f->f0tp1[i] = T();
                this->f->f1t[i] = T();   this->f->f1tp1[i] = T();
                this->f->f2t[i] = T();   this->f->f2tp1[i] = T();
                this->f->f3t[i] = T();   this->f->f3tp1[i] = T();
                this->f->f4t[i] = T();   this->f->f4tp1[i] = T();
                this->f->f5t[i] = T();   this->f->f5tp1[i] = T();
                this->f->f6t[i] = T();   this->f->f6tp1[i] = T();
                this->f->f7t[i] = T();   this->f->f7tp1[i] = T();
                this->f->f8t[i] = T();   this->f->f8tp1[i] = T();
            }
        }
        this->direction = !this->direction;
    }


    template<class T>
    T NSAdjointd2q9<T>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->nx*this->ny*(this->tmax - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetUx(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->ux[this->nx*this->ny*(this->tmax - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetUy(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->uy[this->nx*this->ny*(this->tmax - 1) + this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetQ(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->q[this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetVx(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->vx[this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetVy(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->vy[this->ny*_i + _j];
    }


    template<class T>
    T NSAdjointd2q9<T>::GetSensitivity(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->sensitivity[this->ny*_i + _j];
    }
}