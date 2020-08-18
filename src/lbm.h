//*****************************************************************************
//  Title       :   src/lbm.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/12
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <cassert>


namespace PANSLBM2 {
    enum BOUNDARYTYPE {
        PERIODIC, INLET, OUTLET, BARRIER, MIRROR,
    };

    template<class T>
    class LBM {
public:
        LBM() = delete;
        LBM(int _nx, int _ny, T _viscosity);
        ~LBM();

        template<class F>
        void SetBarrier(F _f);
        template<class F>
        void SetBoundary(F _f);
        template<class F>
        void SetPermeation(F _f);

        void Inlet(T _u, T _v);
        void Stream();
        void UpdateMacro();
        void Collision();
        void ExternalForce();

        bool GetBarrier(int _i, int _j);
        T GetRho(int _i, int _j);
        T GetU(int _i, int _j);
        T GetV(int _i, int _j);
        T GetPermeation(int _i, int _j);

        const int nx, ny;

protected:
        T dx, dt, omega, t0, t1, t2;
        T *f0t, *f1t, *f2t, *f3t, *f4t, *f5t, *f6t, *f7t, *f8t, *f0tp1, *f1tp1, *f2tp1, *f3tp1, *f4tp1, *f5tp1, *f6tp1, *f7tp1, *f8tp1, *rho, *u, *v, *permeation;
        bool *barrier0, *barrier1, *barrier2, *barrier3, *barrier4, *barrier5, *barrier6, *barrier7, *barrier8;
        BOUNDARYTYPE *btxmin, *btxmax, *btymin, *btymax;
    };


    template<class T>
    LBM<T>::LBM(int _nx, int _ny, T _viscosity) : nx(_nx), ny(_ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;     this->dt = 1.0;     this->omega = 1.0/(3.0*_viscosity*this->dt/pow(this->dx, 2.0) + 0.5);
        this->t0 = 4.0/9.0; this->t1 = 1.0/9.0; this->t2 = 1.0/36.0;

        this->f0t = new T[this->nx*this->ny];   this->f0tp1 = new T[this->nx*this->ny];     this->barrier0 = new bool[this->nx*this->ny];
        this->f1t = new T[this->nx*this->ny];   this->f1tp1 = new T[this->nx*this->ny];     this->barrier1 = new bool[this->nx*this->ny];        
        this->f2t = new T[this->nx*this->ny];   this->f2tp1 = new T[this->nx*this->ny];     this->barrier2 = new bool[this->nx*this->ny];     
        this->f3t = new T[this->nx*this->ny];   this->f3tp1 = new T[this->nx*this->ny];     this->barrier3 = new bool[this->nx*this->ny];    
        this->f4t = new T[this->nx*this->ny];   this->f4tp1 = new T[this->nx*this->ny];     this->barrier4 = new bool[this->nx*this->ny];
        this->f5t = new T[this->nx*this->ny];   this->f5tp1 = new T[this->nx*this->ny];     this->barrier5 = new bool[this->nx*this->ny]; 
        this->f6t = new T[this->nx*this->ny];   this->f6tp1 = new T[this->nx*this->ny];     this->barrier6 = new bool[this->nx*this->ny];
        this->f7t = new T[this->nx*this->ny];   this->f7tp1 = new T[this->nx*this->ny];     this->barrier7 = new bool[this->nx*this->ny];
        this->f8t = new T[this->nx*this->ny];   this->f8tp1 = new T[this->nx*this->ny];     this->barrier8 = new bool[this->nx*this->ny];

        this->rho = new T[this->nx*this->ny];   this->u = new T[this->nx*this->ny];         this->v = new T[this->nx*this->ny];

        this->permeation = new T[this->nx*this->ny]();

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->f0t[i] = t0;  this->f0tp1[i] = T();   this->barrier0[i] = false;
            this->f1t[i] = t1;  this->f1tp1[i] = T();   this->barrier1[i] = false;
            this->f2t[i] = t1;  this->f2tp1[i] = T();   this->barrier2[i] = false;
            this->f3t[i] = t1;  this->f3tp1[i] = T();   this->barrier3[i] = false;
            this->f4t[i] = t1;  this->f4tp1[i] = T();   this->barrier4[i] = false;
            this->f5t[i] = t2;  this->f5tp1[i] = T();   this->barrier5[i] = false;
            this->f6t[i] = t2;  this->f6tp1[i] = T();   this->barrier6[i] = false;
            this->f7t[i] = t2;  this->f7tp1[i] = T();   this->barrier7[i] = false;
            this->f8t[i] = t2;  this->f8tp1[i] = T();   this->barrier8[i] = false;

            this->rho[i] = T(); this->u[i] = T();       this->v[i] = T();

            this->permeation[i] = T();
        }

        this->btxmin = new BOUNDARYTYPE[this->ny];  this->btxmax = new BOUNDARYTYPE[this->ny];
        this->btymin = new BOUNDARYTYPE[this->nx];  this->btymax = new BOUNDARYTYPE[this->nx];

        for (int i = 0; i < this->ny; i++) {
            this->btxmin[i] = PERIODIC; this->btxmax[i] = PERIODIC;
        }

        for (int i = 0; i < this->nx; i++) {
            this->btymin[i] = PERIODIC; this->btymax[i] = PERIODIC;
        }
    }


    template<class T>
    LBM<T>::~LBM() {
        delete[] this->f0t;     delete[] this->f0tp1;   delete[] this->barrier0;
        delete[] this->f1t;     delete[] this->f1tp1;   delete[] this->barrier1;
        delete[] this->f2t;     delete[] this->f2tp1;   delete[] this->barrier2;
        delete[] this->f3t;     delete[] this->f3tp1;   delete[] this->barrier3;
        delete[] this->f4t;     delete[] this->f4tp1;   delete[] this->barrier4;
        delete[] this->f5t;     delete[] this->f5tp1;   delete[] this->barrier5;
        delete[] this->f6t;     delete[] this->f6tp1;   delete[] this->barrier6;
        delete[] this->f7t;     delete[] this->f7tp1;   delete[] this->barrier7;
        delete[] this->f8t;     delete[] this->f8tp1;   delete[] this->barrier8;

        delete[] this->rho;     delete[] this->u;       delete[] this->v;

        delete[] this->btxmin;  delete[] this->btxmax;
        delete[] this->btymin;  delete[] this->btymax;
    }


    template<class T>
    template<class F>
    void LBM<T>::SetBarrier(F _f) {
        //----------Set barrier0----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                if (_f(i, j)) {
                    this->barrier0[this->ny*i + j] = true;
                }
            }
        }

        //----------Set barrier1~8----------
        for (int i = 1; i < this->nx - 1; i++) {
            for (int j = 1; j < this->ny - 1; j++) {
                this->barrier1[this->ny*i + j] = this->barrier0[this->ny*(i - 1) + j];
                this->barrier2[this->ny*i + j] = this->barrier0[this->ny*i + (j - 1)];
                this->barrier3[this->ny*i + j] = this->barrier0[this->ny*(i + 1) + j];
                this->barrier4[this->ny*i + j] = this->barrier0[this->ny*i + (j + 1)];
                this->barrier5[this->ny*i + j] = this->barrier0[this->ny*(i - 1) + (j - 1)];
                this->barrier6[this->ny*i + j] = this->barrier0[this->ny*(i + 1) + (j - 1)];
                this->barrier7[this->ny*i + j] = this->barrier0[this->ny*(i + 1) + (j + 1)];
                this->barrier8[this->ny*i + j] = this->barrier0[this->ny*(i - 1) + (j + 1)];
            }
        }
    }


    template<class T>
    template<class F>
    void LBM<T>::SetBoundary(F _f) {
        //----------Set x boundary----------
        for (int j = 0; j < this->ny; j++) {
            this->btxmin[j] = _f(0, j);
            this->btxmax[j] = _f(this->nx - 1, j);
        }

        //----------Set y boundary----------
        for (int i = 0; i < this->nx; i++) {
            this->btymin[i] = _f(i, 0);
            this->btymax[i] = _f(i, this->ny - 1);
        }
    }


    template<class T>
    template<class F>
    void LBM<T>::SetPermeation(F _f) {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->permeation[this->ny*i + j] = _f(i, j);
            }
        }
    }


    template<class T>
    void LBM<T>::Inlet(T _u, T _v) {
        for (int j = 0; j < this->ny; j++) {
            T rho0 = (this->f0t[j] + this->f2t[j] + this->f4t[j] + 2.0*(this->f3t[j] + this->f6t[j] + this->f7t[j]))/(1.0 - _u);
            this->f1t[j] = this->f3t[j] + 2.0*rho0*_u/3.0;
            this->f5t[j] = this->f7t[j] - 0.5*(this->f2t[j] - this->f4t[j]) + rho0*_u/6.0 + rho0*_v/2.0;
            this->f8t[j] = this->f6t[j] + 0.5*(this->f2t[j] - this->f4t[j]) + rho0*_u/6.0 - rho0*_v/2.0;
        }
    }


    template<class T>
    void LBM<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->f0t[this->ny*i + j] = this->f0tp1[this->ny*i + j];
                this->f1t[this->ny*i + j] = (i == 0) ? this->f1tp1[this->ny*(this->nx - 1) + j] : this->f1tp1[this->ny*(i - 1) + j];
                this->f2t[this->ny*i + j] = (j == 0) ? this->f2tp1[this->ny*(i + 1) - 1] : this->f2tp1[this->ny*i + (j - 1)];
                this->f3t[this->ny*i + j] = (i == this->nx - 1) ? this->f3tp1[j] : this->f3tp1[this->ny*(i + 1) + j];
                this->f4t[this->ny*i + j] = (j == this->ny - 1) ? this->f4tp1[this->ny*i] : this->f4tp1[this->ny*i + (j + 1)];
                this->f5t[this->ny*i + j] = (i == 0) ? ((j == 0) ? this->f5tp1[this->ny*this->nx - 1] : this->f5tp1[this->ny*(this->nx - 1) + (j - 1)]) : ((j == 0) ? this->f5tp1[this->ny*i - 1] : this->f5tp1[this->ny*(i - 1) + (j - 1)]);
                this->f6t[this->ny*i + j] = (i == this->nx - 1) ? ((j == 0) ? this->f6tp1[this->ny - 1] : this->f6tp1[j - 1]) : ((j == 0) ? this->f6tp1[this->ny*(i + 2) - 1] : this->f6tp1[this->ny*(i + 1) + (j - 1)]);
                this->f7t[this->ny*i + j] = (i == this->nx - 1) ? ((j == this->ny - 1) ? this->f7tp1[0] : this->f7tp1[j + 1]) : ((j == this->ny - 1) ? this->f7tp1[this->ny*(i + 1)] : this->f7tp1[this->ny*(i + 1) + (j + 1)]);
                this->f8t[this->ny*i + j] = (i == 0) ? ((j == this->ny - 1) ? this->f8tp1[this->ny*(this->nx - 1)] : this->f8tp1[this->ny*(this->nx - 1) + (j + 1)]) : ((j == this->ny - 1) ? this->f8tp1[this->ny*(i - 1)] : this->f8tp1[this->ny*(i - 1) + (j + 1)]);
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->nx*this->ny; i++) {
            if (this->barrier1[i]) {    this->f1t[i] = this->f3tp1[i];  }
            if (this->barrier2[i]) {    this->f2t[i] = this->f4tp1[i];  }
            if (this->barrier3[i]) {    this->f3t[i] = this->f1tp1[i];  }
            if (this->barrier4[i]) {    this->f4t[i] = this->f2tp1[i];  }
            if (this->barrier5[i]) {    this->f5t[i] = this->f7tp1[i];  }
            if (this->barrier6[i]) {    this->f6t[i] = this->f8tp1[i];  }
            if (this->barrier7[i]) {    this->f7t[i] = this->f5tp1[i];  }
            if (this->barrier8[i]) {    this->f8t[i] = this->f6tp1[i];  }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        for (int j = 0; j < this->ny; j++) {
            //.....xmin.....
            if (this->btxmin[j] == OUTLET) {
                this->f1t[j] = this->f1t[this->ny + j];
                this->f5t[j] = this->f5t[this->ny + j];
                this->f8t[j] = this->f8t[this->ny + j];
            } else if (this->btxmin[j] == BARRIER) {
                this->f1t[j] = this->f3tp1[j];
                this->f5t[j] = this->f7tp1[j];
                this->f8t[j] = this->f6tp1[j];
            } else if (this->btxmin[j] == MIRROR) {
                this->f1t[j] = this->f3tp1[j];
                this->f5t[j] = this->f6tp1[j];
                this->f8t[j] = this->f7tp1[j];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->f3t[this->ny*(this->nx - 1) + j] = this->f3t[this->ny*(this->nx - 2) + j];
                this->f6t[this->ny*(this->nx - 1) + j] = this->f6t[this->ny*(this->nx - 2) + j];
                this->f7t[this->ny*(this->nx - 1) + j] = this->f7t[this->ny*(this->nx - 2) + j];
            } else if (this->btxmax[j] == BARRIER) {
                this->f3t[this->ny*(this->nx - 1) + j] = this->f1tp1[this->ny*(this->nx - 1) + j];
                this->f6t[this->ny*(this->nx - 1) + j] = this->f8tp1[this->ny*(this->nx - 1) + j];
                this->f7t[this->ny*(this->nx - 1) + j] = this->f5tp1[this->ny*(this->nx - 1) + j];
            } else if (this->btxmax[j] == MIRROR) {
                this->f3t[this->ny*(this->nx - 1) + j] = this->f1tp1[this->ny*(this->nx - 1) + j];
                this->f6t[this->ny*(this->nx - 1) + j] = this->f5tp1[this->ny*(this->nx - 1) + j];
                this->f7t[this->ny*(this->nx - 1) + j] = this->f8tp1[this->ny*(this->nx - 1) + j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->f2t[this->ny*i] = this->f2t[this->ny*i + 1];
                this->f5t[this->ny*i] = this->f5t[this->ny*i + 1];
                this->f6t[this->ny*i] = this->f6t[this->ny*i + 1];
            } else if (this->btymin[i] == BARRIER) {
                this->f2t[this->ny*i] = this->f4tp1[this->ny*i];
                this->f5t[this->ny*i] = this->f7tp1[this->ny*i];
                this->f6t[this->ny*i] = this->f8tp1[this->ny*i];
            } else if (this->btymin[i] == MIRROR) {
                this->f2t[this->ny*i] = this->f4tp1[this->ny*i];
                this->f5t[this->ny*i] = this->f8tp1[this->ny*i];
                this->f6t[this->ny*i] = this->f7tp1[this->ny*i];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->f4t[this->ny*(i + 1) - 1] = this->f4t[this->ny*(i + 1) - 2];
                this->f7t[this->ny*(i + 1) - 1] = this->f7t[this->ny*(i + 1) - 2];
                this->f8t[this->ny*(i + 1) - 1] = this->f8t[this->ny*(i + 1) - 2];
            } else if (this->btymax[i] == BARRIER) {
                this->f4t[this->ny*(i + 1) - 1] = this->f2tp1[this->ny*(i + 1) - 1];
                this->f7t[this->ny*(i + 1) - 1] = this->f5tp1[this->ny*(i + 1) - 1];
                this->f8t[this->ny*(i + 1) - 1] = this->f6tp1[this->ny*(i + 1) - 1];
            } else if (this->btymax[i] == MIRROR) {
                this->f4t[this->ny*(i + 1) - 1] = this->f2tp1[this->ny*(i + 1) - 1];
                this->f7t[this->ny*(i + 1) - 1] = this->f6tp1[this->ny*(i + 1) - 1];
                this->f8t[this->ny*(i + 1) - 1] = this->f5tp1[this->ny*(i + 1) - 1];
            }
        }
    }


    template<class T>
    void LBM<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = this->f0t[i] + this->f1t[i] + this->f2t[i] + this->f3t[i] + this->f4t[i] + this->f5t[i] + this->f6t[i] + this->f7t[i] + this->f8t[i];
            this->u[i] = (this->f1t[i] - this->f3t[i] + this->f5t[i] - this->f6t[i] - this->f7t[i] + this->f8t[i])/this->rho[i];
            this->v[i] = (this->f2t[i] - this->f4t[i] + this->f5t[i] + this->f6t[i] - this->f7t[i] - this->f8t[i])/this->rho[i];
        }
    }


    template<class T>
    void LBM<T>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T u2 = pow(this->u[i], 2.0);
            T v2 = pow(this->v[i], 2.0);
            T u2v2 = u2 + v2;
            T uv = this->u[i]*this->v[i];
            T omu215 = 1.0 - 1.5*u2v2;

            this->f0tp1[i] = (1.0 - this->omega)*this->f0t[i] + this->omega*this->t0*this->rho[i]*omu215;
            this->f1tp1[i] = (1.0 - this->omega)*this->f1t[i] + this->omega*this->t1*this->rho[i]*(omu215 + 3.0*this->u[i] + 4.5*u2);
            this->f2tp1[i] = (1.0 - this->omega)*this->f2t[i] + this->omega*this->t1*this->rho[i]*(omu215 + 3.0*this->v[i] + 4.5*v2);
            this->f3tp1[i] = (1.0 - this->omega)*this->f3t[i] + this->omega*this->t1*this->rho[i]*(omu215 - 3.0*this->u[i] + 4.5*u2);
            this->f4tp1[i] = (1.0 - this->omega)*this->f4t[i] + this->omega*this->t1*this->rho[i]*(omu215 - 3.0*this->v[i] + 4.5*v2);
            this->f5tp1[i] = (1.0 - this->omega)*this->f5t[i] + this->omega*this->t2*this->rho[i]*(omu215 + 3.0*(this->u[i] + this->v[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f6tp1[i] = (1.0 - this->omega)*this->f6t[i] + this->omega*this->t2*this->rho[i]*(omu215 - 3.0*(this->u[i] - this->v[i]) + 4.5*(u2v2 - 2.0*uv));
            this->f7tp1[i] = (1.0 - this->omega)*this->f7t[i] + this->omega*this->t2*this->rho[i]*(omu215 - 3.0*(this->u[i] + this->v[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f8tp1[i] = (1.0 - this->omega)*this->f8t[i] + this->omega*this->t2*this->rho[i]*(omu215 + 3.0*(this->u[i] - this->v[i]) + 4.5*(u2v2 - 2.0*uv));
        }
    }


    template<class T>
    void LBM<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T tmprho = this->f0t[i] + this->f1t[i] + this->f2t[i] + this->f3t[i] + this->f4t[i] + this->f5t[i] + this->f6t[i] + this->f7t[i] + this->f8t[i];
            T tmpu = (this->f1t[i] - this->f3t[i] + this->f5t[i] - this->f6t[i] - this->f7t[i] + this->f8t[i])/tmprho;
            T tmpv = (this->f2t[i] - this->f4t[i] + this->f5t[i] + this->f6t[i] - this->f7t[i] - this->f8t[i])/tmprho;

            T dxt1alpha = -3.0*this->dx*this->t1*this->permeation[i];
            T dxt2alpha = -3.0*this->dx*this->t2*this->permeation[i];

            this->f1t[i] += dxt1alpha*tmpu;
            this->f2t[i] += dxt1alpha*tmpv;
            this->f3t[i] += dxt1alpha*(-tmpu);
            this->f4t[i] += dxt1alpha*(-tmpv);
            this->f5t[i] += dxt2alpha*(tmpu + tmpv);
            this->f6t[i] += dxt2alpha*(-tmpu + tmpv);
            this->f7t[i] += dxt2alpha*(-tmpu - tmpv);
            this->f8t[i] += dxt2alpha*(tmpu - tmpv);
        }
    }


    template<class T>
    bool LBM<T>::GetBarrier(int _i, int _j) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->barrier0[this->ny*_i + _j];
    }


    template<class T>
    T LBM<T>::GetRho(int _i, int _j) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->ny*_i + _j];
    }


    template<class T>
    T LBM<T>::GetU(int _i, int _j) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[this->ny*_i + _j];
    }


    template<class T>
    T LBM<T>::GetV(int _i, int _j) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->v[this->ny*_i + _j];
    }


    template<class T>
    T LBM<T>::GetPermeation(int _i, int _j) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->permeation[this->ny*_i + _j];
    }
}