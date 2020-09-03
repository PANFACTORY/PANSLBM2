//*****************************************************************************
//  Title       :   src/thermal.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/03
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "lbm.h"


namespace PANSLBM2 {
    template<class T>
    class ThermalLBM : public LBM<T> {
public:
        ThermalLBM() = delete;
        ThermalLBM(int _nx, int _ny, T _viscosity);
        ~ThermalLBM();

        void Inlet(T _u, T _v);
        void Stream();
        void UpdateMacro();
        void Collision();
        void ExternalForce();

        void CaptureMacro(const LBM<T>& _lbm);

        T GetTemperature(int _i, int _j) const;

private:
        T *g0t, *g1t, *g2t, *g3t, *g4t, *g5t, *g6t, *g7t, *g8t, *g0tp1, *g1tp1, *g2tp1, *g3tp1, *g4tp1, *g5tp1, *g6tp1, *g7tp1, *g8tp1;
    };


    template<class T>
    ThermalLBM<T>::ThermalLBM(int _nx, int _ny, T _viscosity) : LBM<T>::LBM(_nx, _ny, _viscosity) {
        this->g0t = new T[this->nx*this->ny];   this->g0tp1 = new T[this->nx*this->ny];
        this->g1t = new T[this->nx*this->ny];   this->g1tp1 = new T[this->nx*this->ny];        
        this->g2t = new T[this->nx*this->ny];   this->g2tp1 = new T[this->nx*this->ny];     
        this->g3t = new T[this->nx*this->ny];   this->g3tp1 = new T[this->nx*this->ny];    
        this->g4t = new T[this->nx*this->ny];   this->g4tp1 = new T[this->nx*this->ny];
        this->g5t = new T[this->nx*this->ny];   this->g5tp1 = new T[this->nx*this->ny]; 
        this->g6t = new T[this->nx*this->ny];   this->g6tp1 = new T[this->nx*this->ny];
        this->g7t = new T[this->nx*this->ny];   this->g7tp1 = new T[this->nx*this->ny];
        this->g8t = new T[this->nx*this->ny];   this->g8tp1 = new T[this->nx*this->ny];
        
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->g0t[i] = t0;  this->g0tp1[i] = T();
            this->g1t[i] = t1;  this->g1tp1[i] = T();
            this->g2t[i] = t1;  this->g2tp1[i] = T();
            this->g3t[i] = t1;  this->g3tp1[i] = T();
            this->g4t[i] = t1;  this->g4tp1[i] = T();
            this->g5t[i] = t2;  this->g5tp1[i] = T();
            this->g6t[i] = t2;  this->g6tp1[i] = T();
            this->g7t[i] = t2;  this->g7tp1[i] = T();
            this->g8t[i] = t2;  this->g8tp1[i] = T();
        }
    }


    template<class T>
    ThermalLBM<T>::~ThermalLBM() {
        delete[] this->g0t;     delete[] this->g0tp1;
        delete[] this->g1t;     delete[] this->g1tp1;
        delete[] this->g2t;     delete[] this->g2tp1;
        delete[] this->g3t;     delete[] this->g3tp1;
        delete[] this->g4t;     delete[] this->g4tp1;
        delete[] this->g5t;     delete[] this->g5tp1;
        delete[] this->g6t;     delete[] this->g6tp1;
        delete[] this->g7t;     delete[] this->g7tp1;
        delete[] this->g8t;     delete[] this->g8tp1;
    }


    template<class T>
    void ThermalLBM<T>::Inlet(T _u, T _v) {
        for (int j = 0; j < this->ny; j++) {
            T term2 = (4.0*_u*this->f1t[j] + (_u + 3.0*_v)*this->f5t[j] + (_u - 3.0*_v)*this->f8t[j])/(3.0*(1.0 - _u));
            this->f3t[j] = this->f1t[j] + term2;
            this->f6t[j] = this->f8t[j] + term2;
            this->f7t[j] = this->f5t[j] + term2;
        }
    }


    template<class T>
    void ThermalLBM<T>::Stream() {
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
    void ThermalLBM<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            int ii = this->nx*this->ny*this->t + i;

            T sensx = this->t1*(this->f1t[i] - this->f3t[i]) + this->t2*(this->f5t[i] - this->f6t[i] - this->f7t[i] + this->f8t[i]);
            T sensy = this->t1*(this->f2t[i] - this->f4t[i]) + this->t2*(this->f5t[i] + this->f6t[i] - this->f7t[i] - this->f8t[i]);

            this->sensitivity[i] += -3.0*this->dx*(sensx*this->u[ii] + sensy*this->v[ii]) + this->u[ii];
        }
    }


    template<class T>
    void ThermalLBM<T>::Collision() {
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
            T f1eq = this->t1*rhoii*(omu215 + 3.0*uii + 4.5*u2);
            T f2eq = this->t1*rhoii*(omu215 + 3.0*vii + 4.5*v2);
            T f3eq = this->t1*rhoii*(omu215 - 3.0*uii + 4.5*u2);
            T f4eq = this->t1*rhoii*(omu215 - 3.0*vii + 4.5*v2);
            T f5eq = this->t2*rhoii*(omu215 + 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
            T f6eq = this->t2*rhoii*(omu215 - 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));
            T f7eq = this->t2*rhoii*(omu215 - 3.0*(uii + vii) + 4.5*(u2v2 + 2.0*uv));
            T f8eq = this->t2*rhoii*(omu215 + 3.0*(uii - vii) + 4.5*(u2v2 - 2.0*uv));

            T faseq1 = (f0eq*this->f0t[i] + f1eq*this->f1t[i] + f2eq*this->f2t[i] + f3eq*this->f3t[i] + f4eq*this->f4t[i] + f5eq*this->f5t[i] + f6eq*this->f6t[i] + f7eq*this->f7t[i] + f8eq*this->f8t[i])/this->rho[ii];
            T w0fas0 = this->t0*this->f0t[i];
            T w1fas1 = this->t1*this->f1t[i];
            T w1fas2 = this->t1*this->f2t[i];
            T w1fas3 = this->t1*this->f3t[i];
            T w1fas4 = this->t1*this->f4t[i];
            T w2fas5 = this->t2*this->f5t[i];
            T w2fas6 = this->t2*this->f6t[i];
            T w2fas7 = this->t2*this->f7t[i];
            T w2fas8 = this->t2*this->f8t[i];
            T faseq2x = 3.0*(-w0fas0*uii + w1fas1*(1.0 + 2.0*uii) - w1fas2*uii + w1fas3*(-1.0 + 2.0*uii) - w1fas4*uii + w2fas5*(1.0 + 2.0*uii + 3.0*vii) + w2fas6*(-1.0 + 2.0*uii - 3.0*vii) + w2fas7*(-1.0 + 2.0*uii + 3.0*vii) + w2fas8*(1.0 + 2.0*uii - 3.0*vii));
            T faseq2y = 3.0*(-w0fas0*vii - w1fas1*vii + w1fas2*(1.0 + 2.0*vii) - w1fas3*vii + w1fas4*(-1.0 + 2.0*vii) + w2fas5*(1.0 + 3.0*uii + 2.0*vii) + w2fas6*(1.0 - 3.0*uii + 2.0*vii) + w2fas7*(-1.0 + 3.0*uii + 2.0*vii) + w2fas8*(-1.0 - 3.0*uii + 2.0*vii));

            this->f0tp1[i] = (1.0 - this->omega)*this->f0t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(-this->v[ii]));
            this->f1tp1[i] = (1.0 - this->omega)*this->f1t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(-this->v[ii]));
            this->f2tp1[i] = (1.0 - this->omega)*this->f2t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
            this->f3tp1[i] = (1.0 - this->omega)*this->f3t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(-this->v[ii]));
            this->f4tp1[i] = (1.0 - this->omega)*this->f4t[i] + this->omega*(faseq1 + faseq2x*(-this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
            this->f5tp1[i] = (1.0 - this->omega)*this->f5t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
            this->f6tp1[i] = (1.0 - this->omega)*this->f6t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(1.0 - this->v[ii]));
            this->f7tp1[i] = (1.0 - this->omega)*this->f7t[i] + this->omega*(faseq1 + faseq2x*(-1.0 - this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
            this->f8tp1[i] = (1.0 - this->omega)*this->f8t[i] + this->omega*(faseq1 + faseq2x*(1.0 - this->u[ii]) + faseq2y*(-1.0 - this->v[ii]));
        }
    }


    template<class T>
    void ThermalLBM<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            int ii = this->nx*this->ny*this->t + i;

            T ef1x = -3.0*this->dx*this->permeation[i]*(this->t1*(this->f1t[i] - this->f3t[i]) + this->t2*(this->f5t[i] - this->f6t[i] - this->f7t[i] + this->f8t[i]));
            T ef1y = -3.0*this->dx*this->permeation[i]*(this->t1*(this->f2t[i] - this->f4t[i]) + this->t2*(this->f5t[i] + this->f6t[i] - this->f7t[i] - this->f8t[i]));
        
            T ubyrho = this->u[ii]/this->rho[ii];
            T vbyrho = this->v[ii]/this->rho[ii];
            T onebyrho = 1.0/this->rho[ii];

            this->f0t[i] += ef1x*(-ubyrho) + ef1y*(-vbyrho) + this->permeation[i]*(-ubyrho);
            this->f1t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(-vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
            this->f2t[i] += ef1x*(-ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(-ubyrho);
            this->f3t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
            this->f4t[i] += ef1x*(-ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(-ubyrho);
            this->f5t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
            this->f6t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(onebyrho - vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
            this->f7t[i] += ef1x*(-onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(-onebyrho - ubyrho);
            this->f8t[i] += ef1x*(onebyrho - ubyrho) + ef1y*(-onebyrho - vbyrho) + this->permeation[i]*(onebyrho - ubyrho);
        }
    }


    template<class T>
    void ThermalLBM<T>::CaptureMacro(const LBM<T>& _lbm) {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                int ii = this->nx*this->ny*this->t + this->ny*i + j;
                this->rho[ii] = _lbm.GetRho(i, j);  this->u[ii] = _lbm.GetU(i, j);  this->v[ii] = _lbm.GetV(i, j);
            }
        }
    }


    template<class T>
    T ThermalLBM<T>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->nx*this->ny*this->t + this->ny*_i + _j];
    }


    template<class T>
    T ThermalLBM<T>::GetU(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[this->nx*this->ny*this->t + this->ny*_i + _j];
    }


    template<class T>
    T ThermalLBM<T>::GetV(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->v[this->nx*this->ny*this->t + this->ny*_i + _j];
    }


    template<class T>
    T ThermalLBM<T>::GetSensitivity(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->sensitivity[this->ny*_i + _j];
    }
}