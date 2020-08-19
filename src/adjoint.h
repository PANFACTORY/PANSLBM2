//*****************************************************************************
//  Title       :   src/adjoint.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/08/18
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <cassert>


#include "lbm.h"


namespace PANSLBM2 {
    template<class T>
    class AdjointLBM : public LBM<T> {
public:
        AdjointLBM() = delete;
        AdjointLBM(const LBM<T>& _lbm);
        ~AdjointLBM() {};

        void Inlet(T _u, T _v);
        void Stream();
        void UpdateMacro();
        void Collision();
        void ExternalForce();
    };


    template<class T>
    AdjointLBM<T>::AdjointLBM(const LBM<T>& _lbm) : LBM<T>::LBM(_lbm) {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->f0t[i] = T(); this->f0tp1[i] = T();
            this->f1t[i] = T(); this->f1tp1[i] = T();
            this->f2t[i] = T(); this->f2tp1[i] = T();
            this->f3t[i] = T(); this->f3tp1[i] = T();
            this->f4t[i] = T(); this->f4tp1[i] = T();
            this->f5t[i] = T(); this->f5tp1[i] = T();
            this->f6t[i] = T(); this->f6tp1[i] = T();
            this->f7t[i] = T(); this->f7tp1[i] = T();
            this->f8t[i] = T(); this->f8tp1[i] = T();
        }
    }


    template<class T>
    void AdjointLBM<T>::Inlet(T _u, T _v) {
        for (int j = 0; j < this->ny; j++) {
            T term2 = (4.0*_u*this->f1t[j] + (_u + 3.0*_v)*this->f5t[j] + (_u - 3.0*_v)*this->f8t[j])/(3.0*(1.0 - _u));
            this->f3t[j] = this->f1t[j] + term2;
            this->f6t[j] = this->f8t[j] + term2;
            this->f7t[j] = this->f5t[j] + term2;
        }
    }


    template<class T>
    void AdjointLBM<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->f0t[this->ny*i + j] = this->f0tp1[this->ny*i + j];
                this->f1t[this->ny*i + j] = (i == this->nx - 1) ? this->f1tp1[j] : this->f1tp1[this->ny*(i + 1) + j];
                this->f2t[this->ny*i + j] = (j == this->ny - 1) ? this->f2tp1[this->ny*i] : this->f2tp1[this->ny*i + (j + 1)];
                this->f3t[this->ny*i + j] = (i == 0) ? this->f3tp1[this->ny*(this->nx - 1) + j] : this->f3tp1[this->ny*(i - 1) + j];
                this->f4t[this->ny*i + j] = (j == 0) ? this->f4tp1[this->ny*i + (this->ny - 1)] : this->f4tp1[this->ny*i + (j - 1)];
                this->f5t[this->ny*i + j] = (i == this->nx - 1) ? ((j == this->ny - 1) ? this->f5tp1[0] : this->f5tp1[j + 1]) : ((j == this->ny - 1) ? this->f5tp1[this->ny*(i + 1)] : this->f5tp1[this->ny*(i + 1) + (j + 1)]);
                this->f6t[this->ny*i + j] = (i == 0) ? ((j == this->ny - 1) ? this->f6tp1[this->ny*(this->nx - 1)] : this->f6tp1[this->ny*(this->nx - 1) + (j + 1)]) : ((j == this->ny - 1) ? this->f6tp1[this->ny*(i - 1)] : this->f6tp1[this->ny*(i - 1) + (j + 1)]);
                this->f7t[this->ny*i + j] = (i == 0) ? ((j == 0) ? this->f7tp1[this->nx*this->ny - 1] : this->f7tp1[this->ny*(this->nx - 1) + (j - 1)]) : ((j == 0) ? this->f7tp1[this->ny*i - 1] : this->f7tp1[this->ny*(i - 1) + (j - 1)]);
                this->f8t[this->ny*i + j] = (i == this->nx - 1) ? ((j == 0) ? this->f8tp1[this->ny - 1] : this->f8tp1[j - 1]) : ((j == 0) ? this->f8tp1[this->ny*(i + 2) - 1] : this->f8tp1[this->ny*(i + 1) + (j - 1)]);
            }
        }

        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->nx*this->ny; i++) {
            if (this->barrier1[i]) {    this->f3t[i] = this->f1tp1[i];  }
            if (this->barrier2[i]) {    this->f4t[i] = this->f2tp1[i];  }
            if (this->barrier3[i]) {    this->f1t[i] = this->f3tp1[i];  }
            if (this->barrier4[i]) {    this->f2t[i] = this->f4tp1[i];  }
            if (this->barrier5[i]) {    this->f7t[i] = this->f5tp1[i];  }
            if (this->barrier6[i]) {    this->f8t[i] = this->f6tp1[i];  }
            if (this->barrier7[i]) {    this->f5t[i] = this->f7tp1[i];  }
            if (this->barrier8[i]) {    this->f6t[i] = this->f8tp1[i];  }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        for (int j = 0; j < this->ny; j++) {
            //.....xmin.....
            if (this->btxmin[j] == OUTLET) {
                this->f3t[j] = this->f3t[this->ny + j];
                this->f6t[j] = this->f6t[this->ny + j];
                this->f7t[j] = this->f7t[this->ny + j];
            } else if (this->btxmin[j] == BARRIER) {
                this->f3t[j] = this->f1tp1[j];
                this->f6t[j] = this->f8tp1[j];
                this->f7t[j] = this->f5tp1[j];
            } else if (this->btxmin[j] == MIRROR) {
                this->f3t[j] = this->f1tp1[j];
                this->f6t[j] = this->f5tp1[j];
                this->f7t[j] = this->f8tp1[j];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->f1t[this->ny*(this->nx - 1) + j] = this->f1t[this->ny*(this->nx - 2) + j];
                this->f5t[this->ny*(this->nx - 1) + j] = this->f5t[this->ny*(this->nx - 2) + j];
                this->f8t[this->ny*(this->nx - 1) + j] = this->f8t[this->ny*(this->nx - 2) + j];
            } else if (this->btxmax[j] == BARRIER) {
                this->f1t[this->ny*(this->nx - 1) + j] = this->f3tp1[this->ny*(this->nx - 1) + j];
                this->f5t[this->ny*(this->nx - 1) + j] = this->f7tp1[this->ny*(this->nx - 1) + j];
                this->f8t[this->ny*(this->nx - 1) + j] = this->f6tp1[this->ny*(this->nx - 1) + j];
            } else if (this->btxmax[j] == MIRROR) {
                this->f1t[this->ny*(this->nx - 1) + j] = this->f3tp1[this->ny*(this->nx - 1) + j];
                this->f5t[this->ny*(this->nx - 1) + j] = this->f6tp1[this->ny*(this->nx - 1) + j];
                this->f8t[this->ny*(this->nx - 1) + j] = this->f7tp1[this->ny*(this->nx - 1) + j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->f4t[this->ny*i] = this->f4t[this->ny*i + 1];
                this->f7t[this->ny*i] = this->f7t[this->ny*i + 1];
                this->f8t[this->ny*i] = this->f8t[this->ny*i + 1];
            } else if (this->btymin[i] == BARRIER) {
                this->f4t[this->ny*i] = this->f2tp1[this->ny*i];
                this->f7t[this->ny*i] = this->f5tp1[this->ny*i];
                this->f8t[this->ny*i] = this->f6tp1[this->ny*i];
            } else if (this->btymin[i] == MIRROR) {
                this->f4t[this->ny*i] = this->f2tp1[this->ny*i];
                this->f7t[this->ny*i] = this->f6tp1[this->ny*i];
                this->f8t[this->ny*i] = this->f5tp1[this->ny*i];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->f2t[this->ny*(i + 1) - 1] = this->f2t[this->ny*(i + 1) - 2];
                this->f5t[this->ny*(i + 1) - 1] = this->f5t[this->ny*(i + 1) - 2];
                this->f6t[this->ny*(i + 1) - 1] = this->f6t[this->ny*(i + 1) - 2];
            } else if (this->btymax[i] == BARRIER) {
                this->f2t[this->ny*(i + 1) - 1] = this->f4tp1[this->ny*(i + 1) - 1];
                this->f5t[this->ny*(i + 1) - 1] = this->f7tp1[this->ny*(i + 1) - 1];
                this->f6t[this->ny*(i + 1) - 1] = this->f8tp1[this->ny*(i + 1) - 1];
            } else if (this->btymax[i] == MIRROR) {
                this->f2t[this->ny*(i + 1) - 1] = this->f4tp1[this->ny*(i + 1) - 1];
                this->f5t[this->ny*(i + 1) - 1] = this->f8tp1[this->ny*(i + 1) - 1];
                this->f6t[this->ny*(i + 1) - 1] = this->f7tp1[this->ny*(i + 1) - 1];
            }
        }
    }


    template<class T>
    void AdjointLBM<T>::UpdateMacro() {

    }


    template<class T>
    void AdjointLBM<T>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T u2 = pow(_u[i], 2.0);
            T v2 = pow(_v[i], 2.0);
            T u2v2 = u2 + v2;
            T uv = _u[i]*_v[i];
            T omu215 = 1.0 - 1.5*u2v2;

            T f0eq = this->t0*_rho[i]*omu215;
            T f1eq = this->t1*_rho[i]*(omu215 + 3.0*_u[i] + 4.5*u2);
            T f2eq = this->t1*_rho[i]*(omu215 + 3.0*_v[i] + 4.5*v2);
            T f3eq = this->t1*_rho[i]*(omu215 - 3.0*_u[i] + 4.5*u2);
            T f4eq = this->t1*_rho[i]*(omu215 - 3.0*_v[i] + 4.5*v2);
            T f5eq = this->t2*_rho[i]*(omu215 + 3.0*(_u[i] + _v[i]) + 4.5*(u2v2 + 2.0*uv));
            T f6eq = this->t2*_rho[i]*(omu215 - 3.0*(_u[i] - _v[i]) + 4.5*(u2v2 - 2.0*uv));
            T f7eq = this->t2*_rho[i]*(omu215 - 3.0*(_u[i] + _v[i]) + 4.5*(u2v2 + 2.0*uv));
            T f8eq = this->t2*_rho[i]*(omu215 + 3.0*(_u[i] - _v[i]) + 4.5*(u2v2 - 2.0*uv));

            T faseq1 = (f0eq*this->f0t[i] + f1eq*this->f1t[i] + f2eq*this->f2t[i] + f3eq*this->f3t[i] + f4eq*this->f4t[i] + f5eq*this->f5t[i] + f6eq*this->f6t[i] + f7eq*this->f7t[i] + f8eq*this->f8t[i])/_rho[i];
            T faseq2x = 3.0*(-this->t0 + 2.0*this->t1 + 8.0*this->t2)*_u[i];
            T faseq2y = 3.0*(-this->t0 + 2.0*this->t1 + 8.0*this->t2)*_v[i];

            this->f0tp1[i] = (1.0 + this->omega)*this->f0t - this->omega*(faseq1 + faseq2x*(-_u[i]) + faseq2y*(-_v[i]));
            this->f1tp1[i] = (1.0 + this->omega)*this->f1t - this->omega*(faseq1 + faseq2x*(1.0 - _u[i]) + faseq2y*(-_v[i]));
            this->f2tp1[i] = (1.0 + this->omega)*this->f2t - this->omega*(faseq1 + faseq2x*(-_u[i]) + faseq2y*(1.0 - _v[i]));
            this->f3tp1[i] = (1.0 + this->omega)*this->f3t - this->omega*(faseq1 + faseq2x*(-1.0 - _u[i]) + faseq2y*(-_v[i]));
            this->f4tp1[i] = (1.0 + this->omega)*this->f4t - this->omega*(faseq1 + faseq2x*(-_u[i]) + faseq2y*(-1.0 - _v[i]));
            this->f5tp1[i] = (1.0 + this->omega)*this->f5t - this->omega*(faseq1 + faseq2x*(1.0 - _u[i]) + faseq2y*(1.0 - _v[i]));
            this->f6tp1[i] = (1.0 + this->omega)*this->f6t - this->omega*(faseq1 + faseq2x*(-1.0 - _u[i]) + faseq2y*(1.0 - _v[i]));
            this->f7tp1[i] = (1.0 + this->omega)*this->f7t - this->omega*(faseq1 + faseq2x*(-1.0 - _u[i]) + faseq2y*(-1.0 - _v[i]));
            this->f8tp1[i] = (1.0 + this->omega)*this->f8t - this->omega*(faseq1 + faseq2x*(1.0 - _u[i]) + faseq2y*(-1.0 - _v[i]));
        }
    }


    template<class T>
    void AdjointLBM<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T ef1x = 3.0*this->dx*this->permeation[i]*(this->t1*(this->f1t[i] - this->f3t[i]) + this->t2*(this->f5t[i] - this->f6t[i] - this->f7t[i] + this->f8t[i]))/_rho[i];
            T ef1y = 3.0*this->dx*this->permeation[i]*(this->t1*(this->f2t[i] - this->f4t[i]) + this->t2*(this->f5t[i] + this->f6t[i] - this->f7t[i] - this->f8t[i]))/_rho[i];
        
            T ubyrho = _u[i]/_rho[i];
            T vbyrho = _v[i]/_rho[i];
            T onebyrho = 1.0/_rho[i]; 

            this->f0t[i] += ef1x*(-ubyrho) + ef1y*(-vbyrho) + this->permeation[i]*(-ubyrho);
            this->f1t[i] += ef1x*(onebyrho - vbyrho) + ef1y*(-vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
            this->f2t[i] += ef1x*(-ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(-ubyrho);
            this->f3t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
            this->f4t[i] += ef1x*(-ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(-ubyrho);
            this->f5t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
            this->f6t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
            this->f7t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
            this->f8t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
        }
    }
}