//*****************************************************************************
//  Title       :   src/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P>
    class NS {
public:
        NS() = delete;
        NS(P<T>* _f, T _viscosity);
        NS(const P<T>& _e);
        virtual ~NS();

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void SetRho(int _i, int _j, T _ux, T _uy);
        virtual void SetUxUy(int _i, int _j, T _ux, T _uy);
        virtual void ExternalForce();

        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _d, int _i, int _j) const;
        
        const int nx, ny;

protected:
        T omega;
        P<T>* f;
        T *rho, *u[P<T>::nd];
    };


    template<class T, template<class>class P>
    NS<T, P>::NS(P<T>* _f, T _viscosity) : nx(_f->nx), ny(_f->ny) {
        assert(0 < _f->nx && 0 < _f->ny);
        this->f = _f;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/this->f->dx*this->f->dx + 0.5);

        this->rho = new T[this->nx*this->ny];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->nx*this->ny];
        }

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = T();
            }
        }
    }


    template<class T, template<class>class P>
    NS<T, P>::NS(const P<T>& _e) : nx(_e.nx), ny(_e.ny), f(_e.f) {
        this->f = _e.f;
        this->omega = _e.omega;

        this->rho = new T[this->nx*this->ny];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->nx*this->ny];
        }

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = _e.rho[i];
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = _e.u[j][i];
            }
        }
    }


    template<class T, template<class>class P>
    NS<T, P>::~NS() {
        delete[] this->rho;
        for (int i = 0; i < P<T>::nd; i++) {
            delete[] this->u[i];
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = T();
            }

            for (int j = 0; j < P<T>::nc; j++) {
                this->rho[i] += this->f->ft[j][i];
                for (int k = 0; k < P<T>::nd; k++) {
                    this->u[k][i] += P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }

            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] /= this->rho[i];
            }
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();

                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][i];
                    uu += this->u[k][i]*this->u[k][i];
                }

                T fieq = P<T>::ei[j]*this->rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*fieq;
            }
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::SetRho(int _i, int _j, T _rho, T _u) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T ux0 = 1.0 - (this->f->ft[0][ij] + this->f->ft[2][ij] + this->f->ft[4][ij] + 2.0*(this->f->ft[3][ij] + this->f->ft[6][ij] + this->f->ft[7][ij]))/_rho;
            this->f->ft[1][ij] = this->f->ft[3][ij] + 2.0*_rho*ux0/3.0;
            this->f->ft[5][ij] = this->f->ft[7][ij] - 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) + _rho*ux0/6.0 + _rho*_u/2.0;
            this->f->ft[8][ij] = this->f->ft[6][ij] + 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) + _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_i == this->nx - 1) {
            T ux0 = -1.0 + (this->f->ft[0][ij] + this->f->ft[2][ij] + this->f->ft[4][ij] + 2.0*(this->f->ft[3][ij] + this->f->ft[5][ij] + this->f->ft[8][ij]))/_rho;
            this->f->ft[3][ij] = this->f->ft[1][ij] - 2.0*_rho*ux0/3.0;
            this->f->ft[6][ij] = this->f->ft[8][ij] - 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) - _rho*ux0/6.0 + _rho*_u/2.0;
            this->f->ft[7][ij] = this->f->ft[6][ij] + 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) - _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_j == 0) {
            T uy0 = 1.0 - (this->f->ft[0][ij] + this->f->ft[1][ij] + this->f->ft[3][ij] + 2.0*(this->f->ft[4][ij] + this->f->ft[7][ij] + this->f->ft[8][ij]))/_rho;
            this->f->ft[2][ij] = this->f->ft[4][ij] + 2.0*_rho*uy0/3.0;
            this->f->ft[5][ij] = this->f->ft[7][ij] - 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) + _rho*_u/2.0 + _rho*uy0/6.0;
            this->f->ft[6][ij] = this->f->ft[8][ij] + 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) - _rho*_u/2.0 + _rho*uy0/6.0;
        } else if (_j == this->ny - 1) {
            T uy0 = -1.0 - (this->f->ft[0][ij] + this->f->ft[1][ij] + this->f->ft[3][ij] + 2.0*(this->f->ft[2][ij] + this->f->ft[5][ij] + this->f->ft[6][ij]))/_rho;
            this->f->ft[4][ij] = this->f->ft[2][ij] - 2.0*_rho*uy0/3.0;
            this->f->ft[7][ij] = this->f->ft[5][ij] + 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) - _rho*_u/2.0 - _rho*uy0/6.0;
            this->f->ft[8][ij] = this->f->ft[6][ij] - 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) + _rho*_u/2.0 - _rho*uy0/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::SetUxUy(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (this->f->ft[0][ij] + this->f->ft[2][ij] + this->f->ft[4][ij] + 2.0*(this->f->ft[3][ij] + this->f->ft[6][ij] + this->f->ft[7][ij]))/(1.0 - _ux);
            this->f->ft[1][ij] = this->f->ft[3][ij] + 2.0*rho0*_ux/3.0;
            this->f->ft[5][ij] = this->f->ft[7][ij] - 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) + rho0*_ux/6.0 + rho0*_uy/2.0;
            this->f->ft[8][ij] = this->f->ft[6][ij] + 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) + rho0*_ux/6.0 - rho0*_uy/2.0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (this->f->ft[0][ij] + this->f->ft[2][ij] + this->f->ft[4][ij] + 2.0*(this->f->ft[1][ij] + this->f->ft[5][ij] + this->f->ft[8][ij]))/(1.0 + _ux);
            this->f->ft[3][ij] = this->f->ft[1][ij] - 2.0*rho0*_ux/3.0;
            this->f->ft[7][ij] = this->f->ft[5][ij] + 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) - rho0*_ux/6.0 - rho0*_uy/2.0;
            this->f->ft[6][ij] = this->f->ft[8][ij] - 0.5*(this->f->ft[2][ij] - this->f->ft[4][ij]) - rho0*_ux/6.0 + rho0*_uy/2.0;
        } else if (_j == 0) {
            T rho0 = (this->f->ft[0][ij] + this->f->ft[1][ij] + this->f->ft[3][ij] + 2.0*(this->f->ft[4][ij] + this->f->ft[7][ij] + this->f->ft[8][ij]))/(1.0 - _uy);
            this->f->ft[2][ij] = this->f->ft[4][ij] + 2.0*rho0*_uy/3.0;
            this->f->ft[5][ij] = this->f->ft[7][ij] - 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) + rho0*_ux/2.0 + rho0*_uy/6.0;
            this->f->ft[6][ij] = this->f->ft[8][ij] + 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) - rho0*_ux/2.0 + rho0*_uy/6.0;
        } else if (_j == this->ny - 1) {
            T rho0 = (this->f->ft[0][ij] + this->f->ft[1][ij] + this->f->ft[3][ij] + 2.0*(this->f->ft[2][ij] + this->f->ft[5][ij] + this->f->ft[6][ij]))/(1.0 + _uy);
            this->f->ft[4][ij] = this->f->ft[2][ij] - 2.0*rho0*_uy/3.0;
            this->f->ft[7][ij] = this->f->ft[5][ij] + 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) - rho0*_ux/2.0 - rho0*_uy/6.0;
            this->f->ft[8][ij] = this->f->ft[6][ij] - 0.5*(this->f->ft[1][ij] - this->f->ft[3][ij]) + rho0*_ux/2.0 - rho0*_uy/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            
        }
    }


    template<class T, template<class>class P>
    T NS<T, P>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->rho[this->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NS<T, P>::GetU(int _d, int _i, int _j) const {
        assert(0 < _d < P<T>::nd && 0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->u[_d][this->ny*_i + _j];
    }
}