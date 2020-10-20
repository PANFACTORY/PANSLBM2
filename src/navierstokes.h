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
        virtual T GetUx(int _i, int _j) const;
        virtual T GetUy(int _i, int _j) const;
        
        const int nx, ny;

protected:
        T omega;
        P<T>* f;
        T *rho, *ux, *uy;
    };


    template<class T, template<class>class P>
    NS<T, P>::NS(P<T>* _f, T _viscosity) : nx(_f->nx), ny(_f->ny) {
        assert(0 < _f->nx && 0 < _f->ny);
        this->f = _f;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/this->f->dx*this->f->dx + 0.5);

        this->rho = new T[this->nx*this->ny];
        this->ux = new T[this->nx*this->ny];
        this->uy = new T[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = T();
            this->ux[i] = T();
            this->uy[i] = T();
        }
    }


    template<class T, template<class>class P>
    NS<T, P>::NS(const P<T>& _e) : nx(_e.nx), ny(_e.ny), f(_e.f) {
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


    template<class T, template<class>class P>
    NS<T, P>::~NS() {
        delete[] this->rho;
        delete[] this->ux;
        delete[] this->uy;
    }


    template<class T, template<class>class P>
    void NS<T, P>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->rho[i] = this->f->ft[0][i] + this->f->ft[1][i] + this->f->ft[2][i] + this->f->ft[3][i] + this->f->ft[4][i] + this->f->ft[5][i] + this->f->ft[6][i] + this->f->ft[7][i] + this->f->ft[8][i];
            this->ux[i] = (this->f->ft[1][i] - this->f->ft[3][i] + this->f->ft[5][i] - this->f->ft[6][i] - this->f->ft[7][i] + this->f->ft[8][i])/this->rho[i];
            this->uy[i] = (this->f->ft[2][i] - this->f->ft[4][i] + this->f->ft[5][i] + this->f->ft[6][i] - this->f->ft[7][i] - this->f->ft[8][i])/this->rho[i];
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T u2 = this->ux[i]*this->ux[i];
            T v2 = this->uy[i]*this->uy[i];
            T u2v2 = u2 + v2;
            T uv = this->ux[i]*this->uy[i];
            T omu215 = 1.0 - 1.5*u2v2;

            this->f->ftp1[0][i] = (1.0 - this->omega)*this->f->ft[0][i] + this->omega*P<T>::ei[0]*this->rho[i]*omu215;
            this->f->ftp1[1][i] = (1.0 - this->omega)*this->f->ft[1][i] + this->omega*P<T>::ei[1]*this->rho[i]*(omu215 + 3.0*this->ux[i] + 4.5*u2);
            this->f->ftp1[2][i] = (1.0 - this->omega)*this->f->ft[2][i] + this->omega*P<T>::ei[2]*this->rho[i]*(omu215 + 3.0*this->uy[i] + 4.5*v2);
            this->f->ftp1[3][i] = (1.0 - this->omega)*this->f->ft[3][i] + this->omega*P<T>::ei[3]*this->rho[i]*(omu215 - 3.0*this->ux[i] + 4.5*u2);
            this->f->ftp1[4][i] = (1.0 - this->omega)*this->f->ft[4][i] + this->omega*P<T>::ei[4]*this->rho[i]*(omu215 - 3.0*this->uy[i] + 4.5*v2);
            this->f->ftp1[5][i] = (1.0 - this->omega)*this->f->ft[5][i] + this->omega*P<T>::ei[5]*this->rho[i]*(omu215 + 3.0*(this->ux[i] + this->uy[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f->ftp1[6][i] = (1.0 - this->omega)*this->f->ft[6][i] + this->omega*P<T>::ei[6]*this->rho[i]*(omu215 - 3.0*(this->ux[i] - this->uy[i]) + 4.5*(u2v2 - 2.0*uv));
            this->f->ftp1[7][i] = (1.0 - this->omega)*this->f->ft[7][i] + this->omega*P<T>::ei[7]*this->rho[i]*(omu215 - 3.0*(this->ux[i] + this->uy[i]) + 4.5*(u2v2 + 2.0*uv));
            this->f->ftp1[8][i] = (1.0 - this->omega)*this->f->ft[8][i] + this->omega*P<T>::ei[8]*this->rho[i]*(omu215 + 3.0*(this->ux[i] - this->uy[i]) + 4.5*(u2v2 - 2.0*uv));
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
    T NS<T, P>::GetUx(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->ux[this->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NS<T, P>::GetUy(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->uy[this->ny*_i + _j];
    }
}