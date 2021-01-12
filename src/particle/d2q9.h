//*****************************************************************************
//  Title       :   src/particle/d2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/21
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    enum BOUNDARYTYPE {
        PERIODIC, OUTLET, BARRIER, MIRROR, OTHER,
    };

    template<class T>
    class D2Q9 {
public:
        D2Q9() = delete;
        D2Q9(int _nx, int _ny);
        D2Q9(const D2Q9<T>& _p);
        ~D2Q9();

        void SetBarrier(int _i, int _j, bool _isbarrier);
        void SetBoundary(int _i, int _j, BOUNDARYTYPE _boundarytype);
        void SetRho(int _i, int _j, T _rho, T _u);                      //  Set boundary condition for NavierStokes  
        void SetU(int _i, int _j, T _ux, T _uy);
        void SetTemperature(int _i, int _j, T _temperature);            //  Set boundary condition for Advection
        void SetFlux(int _i, int _j, T _ux, T _uy, T _q);
        void SetiRho(int _i, int _j);                                   //  Set boundary condition for Adjoint of NavierStokes
        void SetiU(int _i, int _j, T _ux, T _uy);  
        void SetiUPressureDrop(int _i, int _j, T _ux, T _uy, T _eps = 1);
        void SetiTemperature(int _i, int _j);                           //  Set boundary condition for Adjoint of Advection
        void SetiFlux(int _i, int _j, T _ux, T _uy);
        void SetiRhoFlux(const D2Q9<T>& _g, int _i, int _j, T _rho, T _ux, T _uy, T _temperature);

        void Stream();
        void iStream();
        
        bool GetBarrier(int _i, int _j) const;
        BOUNDARYTYPE GetBoundary(int _i, int _j) const;
        int GetIndex(int _i, int _j) const;
        
        const int nx, ny, np;                               //  nx&ny : number of points along x&y coordinate, np : number of all points
        static const int nc = 9, nd = 2, cx[nc], cy[nc];    //  nc : number of particles, nd : number of dimensions
        static const T ei[nc];

        T dx, dt;
        T *ft[nc], *ftp1[nc];
        bool *barrier[nc];
        BOUNDARYTYPE *btxmin, *btxmax, *btymin, *btymax;
    };


    template<class T>
    const int D2Q9<T>::cx[D2Q9<T>::nc] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };


    template<class T>
    const int D2Q9<T>::cy[D2Q9<T>::nc] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };


    template<class T>
    const T D2Q9<T>::ei[D2Q9<T>::nc] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };


    template<class T>
    D2Q9<T>::D2Q9(int _nx, int _ny) : nx(_nx), ny(_ny), np(_nx*_ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;
        this->dt = 1.0;

        for (int i = 0; i < D2Q9<T>::nc; i++) {
            this->ft[i] = new T[this->np];
            this->ftp1[i] = new T[this->np];
            this->barrier[i] = new bool[this->np];

            for (int j = 0; j < this->np; j++) {
                this->ft[i][j] = D2Q9<T>::ei[i];
                this->ftp1[i][j] = T();
                this->barrier[i][j] = false;
            }    
        }

        this->btxmin = new BOUNDARYTYPE[this->ny];
        this->btxmax = new BOUNDARYTYPE[this->ny];
        this->btymin = new BOUNDARYTYPE[this->nx];
        this->btymax = new BOUNDARYTYPE[this->nx];

        for (int i = 0; i < this->ny; i++) {
            this->btxmin[i] = PERIODIC;
            this->btxmax[i] = PERIODIC;
        }

        for (int i = 0; i < this->nx; i++) {
            this->btymin[i] = PERIODIC;
            this->btymax[i] = PERIODIC;
        }
    }


    template<class T>
    D2Q9<T>::D2Q9(const D2Q9<T>& _p) : nx(_p.nx), ny(_p.ny), np(_p.np) {
        this->dx = _p.dx;
        this->dt = _p.dt;
        
        for (int i = 0; i < D2Q9<T>::nc; i++) {
            this->ft[i] = new T[this->np];
            this->ftp1[i] = new T[this->np];
            this->barrier[i] = new bool[this->np];

            for (int j = 0; j < this->np; j++) {
                this->ft[i][j] = _p.ft[i][j];
                this->ftp1[i][j] = _p.ftp1[i][j];
                this->barrier[i][j] = _p.barrier[i][j];
            }    
        }

        this->btxmin = new BOUNDARYTYPE[this->ny];
        this->btxmax = new BOUNDARYTYPE[this->ny];
        this->btymin = new BOUNDARYTYPE[this->nx];
        this->btymax = new BOUNDARYTYPE[this->nx];

        for (int i = 0; i < this->ny; i++) {
            this->btxmin[i] = _p.btxmin[i];
            this->btxmax[i] = _p.btxmax[i];
        }

        for (int i = 0; i < this->nx; i++) {
            this->btymin[i] = _p.btymin[i];
            this->btymax[i] = _p.btymax[i];
        }
    }


    template<class T>
    D2Q9<T>::~D2Q9() {
        for (int i = 0; i < D2Q9<T>::nc; i++) {
            delete[] this->ft[i];
            delete[] this->ftp1[i];
            delete[] this->barrier[i];
        }
        
        delete[] this->btxmin;
        delete[] this->btxmax;
        delete[] this->btymin;
        delete[] this->btymax;
    }


    template<class T>
    void D2Q9<T>::SetBarrier(int _i, int _j, bool _isbarrier) {
        for (int i = 0; i < D2Q9<T>::nc; i++) {
            this->barrier[i][this->ny*(_i + D2Q9<T>::cx[i]) + (_j + D2Q9<T>::cy[i])] = _isbarrier;
        }
    }


    template<class T>
    void D2Q9<T>::SetBoundary(int _i, int _j, BOUNDARYTYPE _boundarytype) {
        if (_i == 0) {
            this->btxmin[_j] = _boundarytype;
        } 
        if (_i == this->nx - 1) {
            this->btxmax[_j] = _boundarytype;
        } 
        if (_j == 0) {
            this->btymin[_i] = _boundarytype;
        } 
        if (_j == this->ny - 1) {
            this->btymax[_i] = _boundarytype;
        } 
        if (_i != 0 && _i != this->nx - 1 && _j != 0 && _j != this->ny - 1) {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetRho(int _i, int _j, T _rho, T _u) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T ux0 = 1.0 - (this->ft[0][ij] + this->ft[2][ij] + this->ft[4][ij] + 2.0*(this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]))/_rho;
            this->ft[1][ij] = this->ft[3][ij] + 2.0*_rho*ux0/3.0;
            this->ft[5][ij] = this->ft[7][ij] - 0.5*(this->ft[2][ij] - this->ft[4][ij]) + _rho*ux0/6.0 + _rho*_u/2.0;
            this->ft[8][ij] = this->ft[6][ij] + 0.5*(this->ft[2][ij] - this->ft[4][ij]) + _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_i == this->nx - 1) {
            T ux0 = -1.0 + (this->ft[0][ij] + this->ft[2][ij] + this->ft[4][ij] + 2.0*(this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]))/_rho;
            this->ft[3][ij] = this->ft[1][ij] - 2.0*_rho*ux0/3.0;
            this->ft[6][ij] = this->ft[8][ij] - 0.5*(this->ft[2][ij] - this->ft[4][ij]) - _rho*ux0/6.0 + _rho*_u/2.0;
            this->ft[7][ij] = this->ft[5][ij] + 0.5*(this->ft[2][ij] - this->ft[4][ij]) - _rho*ux0/6.0 - _rho*_u/2.0;
        } else if (_j == 0) {
            T uy0 = 1.0 - (this->ft[0][ij] + this->ft[1][ij] + this->ft[3][ij] + 2.0*(this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]))/_rho;
            this->ft[2][ij] = this->ft[4][ij] + 2.0*_rho*uy0/3.0;
            this->ft[5][ij] = this->ft[7][ij] - 0.5*(this->ft[1][ij] - this->ft[3][ij]) + _rho*_u/2.0 + _rho*uy0/6.0;
            this->ft[6][ij] = this->ft[8][ij] + 0.5*(this->ft[1][ij] - this->ft[3][ij]) - _rho*_u/2.0 + _rho*uy0/6.0;
        } else if (_j == this->ny - 1) {
            T uy0 = -1.0 - (this->ft[0][ij] + this->ft[1][ij] + this->ft[3][ij] + 2.0*(this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]))/_rho;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*_rho*uy0/3.0;
            this->ft[7][ij] = this->ft[5][ij] + 0.5*(this->ft[1][ij] - this->ft[3][ij]) - _rho*_u/2.0 - _rho*uy0/6.0;
            this->ft[8][ij] = this->ft[6][ij] - 0.5*(this->ft[1][ij] - this->ft[3][ij]) + _rho*_u/2.0 - _rho*uy0/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetU(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (this->ft[0][ij] + this->ft[2][ij] + this->ft[4][ij] + 2.0*(this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]))/(1.0 - _ux);
            this->ft[1][ij] = this->ft[3][ij] + 2.0*rho0*_ux/3.0;
            this->ft[5][ij] = this->ft[7][ij] - 0.5*(this->ft[2][ij] - this->ft[4][ij]) + rho0*_ux/6.0 + rho0*_uy/2.0;
            this->ft[8][ij] = this->ft[6][ij] + 0.5*(this->ft[2][ij] - this->ft[4][ij]) + rho0*_ux/6.0 - rho0*_uy/2.0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (this->ft[0][ij] + this->ft[2][ij] + this->ft[4][ij] + 2.0*(this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]))/(1.0 + _ux);
            this->ft[3][ij] = this->ft[1][ij] - 2.0*rho0*_ux/3.0;
            this->ft[7][ij] = this->ft[5][ij] + 0.5*(this->ft[2][ij] - this->ft[4][ij]) - rho0*_ux/6.0 - rho0*_uy/2.0;
            this->ft[6][ij] = this->ft[8][ij] - 0.5*(this->ft[2][ij] - this->ft[4][ij]) - rho0*_ux/6.0 + rho0*_uy/2.0;
        } else if (_j == 0) {
            T rho0 = (this->ft[0][ij] + this->ft[1][ij] + this->ft[3][ij] + 2.0*(this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]))/(1.0 - _uy);
            this->ft[2][ij] = this->ft[4][ij] + 2.0*rho0*_uy/3.0;
            this->ft[5][ij] = this->ft[7][ij] - 0.5*(this->ft[1][ij] - this->ft[3][ij]) + rho0*_ux/2.0 + rho0*_uy/6.0;
            this->ft[6][ij] = this->ft[8][ij] + 0.5*(this->ft[1][ij] - this->ft[3][ij]) - rho0*_ux/2.0 + rho0*_uy/6.0;
        } else if (_j == this->ny - 1) {
            T rho0 = (this->ft[0][ij] + this->ft[1][ij] + this->ft[3][ij] + 2.0*(this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]))/(1.0 + _uy);
            this->ft[4][ij] = this->ft[2][ij] - 2.0*rho0*_uy/3.0;
            this->ft[7][ij] = this->ft[5][ij] + 0.5*(this->ft[1][ij] - this->ft[3][ij]) - rho0*_ux/2.0 - rho0*_uy/6.0;
            this->ft[8][ij] = this->ft[6][ij] - 0.5*(this->ft[1][ij] - this->ft[3][ij]) + rho0*_ux/2.0 - rho0*_uy/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetTemperature(int _i, int _j, T _temperature) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[7][ij]);
            this->ft[1][ij] = temperature0/9.0;
            this->ft[5][ij] = temperature0/36.0;
            this->ft[8][ij] = temperature0/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[8][ij]);
            this->ft[3][ij] = temperature0/9.0;
            this->ft[6][ij] = temperature0/36.0;
            this->ft[7][ij] = temperature0/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[7][ij] - this->ft[8][ij]);
            this->ft[2][ij] = temperature0/9.0;
            this->ft[5][ij] = temperature0/36.0;
            this->ft[6][ij] = temperature0/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[5][ij] - this->ft[6][ij]);
            this->ft[4][ij] = temperature0/9.0;
            this->ft[7][ij] = temperature0/36.0;
            this->ft[8][ij] = temperature0/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }
    
    
    template<class T>
    void D2Q9<T>::SetFlux(int _i, int _j, T _ux, T _uy, T _q) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T temperature0 = 6.0*(_q + this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/(1.0 - 3.0*_ux);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(-_q + this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/(1.0 + 3.0*_ux);
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_q + this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/(1.0 - 3.0*_uy);
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(-_q + this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/(1.0 + 3.0*_uy);
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiRho(int _i, int _j) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            this->ft[3][ij] = this->ft[1][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0;
            this->ft[6][ij] = this->ft[8][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0;
            this->ft[7][ij] = this->ft[5][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0;
        } else if (_i == this->nx - 1) {
            this->ft[1][ij] = this->ft[3][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0;
            this->ft[5][ij] = this->ft[7][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0;
            this->ft[8][ij] = this->ft[6][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0;
        } else if (_j == 0) {
            this->ft[4][ij] = this->ft[2][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0;
            this->ft[7][ij] = this->ft[5][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0;
            this->ft[8][ij] = this->ft[6][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0;
        } else if (_j == this->ny - 1) {
            this->ft[2][ij] = this->ft[4][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0;
            this->ft[5][ij] = this->ft[7][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0;
            this->ft[6][ij] = this->ft[8][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiU(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (_ux*(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]) + 3.0*_uy*(this->ft[5][ij] - this->ft[8][ij]))/(3.0*(1.0 - _ux));
            this->ft[3][ij] = this->ft[1][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (-_ux*(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]) + 3.0*_uy*(this->ft[6][ij] - this->ft[7][ij]))/(3.0*(1.0 + _ux));
            this->ft[1][ij] = this->ft[3][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == 0) {
            T rho0 = (_uy*(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]) + 3.0*_ux*(this->ft[5][ij] - this->ft[6][ij]))/(3.0*(1.0 - _uy));
            this->ft[4][ij] = this->ft[2][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = (-_uy*(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]) + 3.0*_ux*(this->ft[8][ij] - this->ft[7][ij]))/(3.0*(1.0 + _uy));
            this->ft[2][ij] = this->ft[4][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiUPressureDrop(int _i, int _j, T _ux, T _uy, T _eps) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (-2.0*_eps + _ux*(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]) + 3.0*_uy*(this->ft[5][ij] - this->ft[8][ij]))/(3.0*(1.0 - _ux));
            this->ft[3][ij] = this->ft[1][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (-2.0*_eps - _ux*(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]) + 3.0*_uy*(this->ft[6][ij] - this->ft[7][ij]))/(3.0*(1.0 + _ux));
            this->ft[1][ij] = this->ft[3][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == 0) {
            T rho0 = (-2.0*_eps + _uy*(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]) + 3.0*_ux*(this->ft[5][ij] - this->ft[6][ij]))/(3.0*(1.0 - _uy));
            this->ft[4][ij] = this->ft[2][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = (-2.0*_eps - _uy*(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]) + 3.0*_ux*(this->ft[8][ij] - this->ft[7][ij]))/(3.0*(1.0 + _uy));
            this->ft[2][ij] = this->ft[4][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiTemperature(int _i, int _j) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            this->ft[3][ij] = -(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/6.0;
            this->ft[6][ij] = -(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/6.0;
            this->ft[7][ij] = -(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/6.0;
        } else if (_i == this->nx - 1) {
            this->ft[1][ij] = -(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/6.0;
            this->ft[5][ij] = -(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/6.0;
            this->ft[8][ij] = -(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/6.0;
        } else if (_j == 0) {
            this->ft[4][ij] = -(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/6.0;
            this->ft[7][ij] = -(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/6.0;
            this->ft[8][ij] = -(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/6.0;
        } else if (_j == this->ny - 1) {
            this->ft[2][ij] = -(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/6.0;
            this->ft[5][ij] = -(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/6.0;
            this->ft[6][ij] = -(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiFlux(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (1.0 + 3.0*_ux)/(6.0*(1.0 - 3.0*_ux))*(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]) + _uy/(2.0*(1.0 - 3.0*_ux))*(this->ft[5][ij] - this->ft[8][ij]);
            this->ft[3][ij] = rho0;
            this->ft[6][ij] = rho0;
            this->ft[7][ij] = rho0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (1.0 - 3.0*_ux)/(6.0*(1.0 + 3.0*_ux))*(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]) + _uy/(2.0*(1.0 + 3.0*_ux))*(this->ft[6][ij] - this->ft[7][ij]);
            this->ft[1][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == 0) {
            T rho0 = (1.0 + 3.0*_uy)/(6.0*(1.0 - 3.0*_uy))*(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]) + _ux/(2.0*(1.0 - 3.0*_uy))*(this->ft[5][ij] - this->ft[6][ij]);
            this->ft[4][ij] = rho0;
            this->ft[7][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = (1.0 - 3.0*_uy)/(6.0*(1.0 + 3.0*_uy))*(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]) + _ux/(2.0*(1.0 + 3.0*_uy))*(this->ft[8][ij] - this->ft[7][ij]);
            this->ft[2][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[6][ij] = rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiRhoFlux(const D2Q9<T>& _g, int _i, int _j, T _rho, T _ux, T _uy, T _temperature) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = -(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0;
            T flux0 = _temperature/(3.0*(1.0 - 3.0*_ux)*_rho)*(4.0*_g.ft[1][ij] + _g.ft[5][ij] + _g.ft[8][ij]) + _uy*_temperature/(2.0*(1.0 - 3.0*_ux)*_rho)*(_g.ft[5][ij] - _g.ft[8][ij]);
            this->ft[3][ij] = this->ft[1][ij] + rho0 - flux0;
            this->ft[6][ij] = this->ft[8][ij] + rho0 - flux0;
            this->ft[7][ij] = this->ft[5][ij] + rho0 - flux0;
        } else if (_i == this->nx - 1) {
            T rho0 = -(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0;
            T flux0 = _temperature/(3.0*(1.0 + 3.0*_ux)*_rho)*(4.0*_g.ft[3][ij] + _g.ft[6][ij] + _g.ft[7][ij]) + _uy*_temperature/(2.0*(1.0 + 3.0*_ux)*_rho)*(_g.ft[6][ij] - _g.ft[7][ij]);
            this->ft[1][ij] = this->ft[3][ij] + rho0 - flux0;
            this->ft[5][ij] = this->ft[7][ij] + rho0 - flux0;
            this->ft[8][ij] = this->ft[6][ij] + rho0 - flux0;
        } else if (_j == 0) {
            T rho0 = -(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0;
            T flux0 = _temperature/(3.0*(1.0 - 3.0*_uy)*_rho)*(4.0*_g.ft[2][ij] + _g.ft[5][ij] + _g.ft[6][ij]) + _ux*_temperature/(2.0*(1.0 - 3.0*_uy)*_rho)*(_g.ft[5][ij] - _g.ft[6][ij]);
            this->ft[4][ij] = this->ft[2][ij] + rho0 - flux0;
            this->ft[7][ij] = this->ft[5][ij] + rho0 - flux0;
            this->ft[8][ij] = this->ft[6][ij] + rho0 - flux0;
        } else if (_j == this->ny - 1) {
            T rho0 = -(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0;
            T flux0 = _temperature/(3.0*(1.0 + 3.0*_uy)*_rho)*(4.0*_g.ft[4][ij] + _g.ft[7][ij] + _g.ft[8][ij]) + _ux*_temperature/(2.0*(1.0 + 3.0*_uy)*_rho)*(_g.ft[8][ij] - _g.ft[7][ij]);
            this->ft[2][ij] = this->ft[4][ij] + rho0 - flux0;
            this->ft[5][ij] = this->ft[7][ij] + rho0 - flux0;
            this->ft[6][ij] = this->ft[8][ij] + rho0 - flux0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                for (int k = 0; k < D2Q9<T>::nc; k++) {
                    int ip1 = i - D2Q9<T>::cx[k] == -1 ? this->nx - 1 : (i - D2Q9<T>::cx[k] == this->nx ? 0 : i - D2Q9<T>::cx[k]);
                    int jp1 = j - D2Q9<T>::cy[k] == -1 ? this->ny - 1 : (j - D2Q9<T>::cy[k] == this->ny ? 0 : j - D2Q9<T>::cy[k]);
                    this->ft[k][this->ny*i + j] = this->ftp1[k][this->ny*ip1 + jp1];
                }
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->np; i++) {
            if (this->barrier[1][i]) {    
                this->ft[1][i] = this->ftp1[3][i];  
            }
            if (this->barrier[2][i]) {    
                this->ft[2][i] = this->ftp1[4][i];  
            }
            if (this->barrier[3][i]) {    
                this->ft[3][i] = this->ftp1[1][i];  
            }
            if (this->barrier[4][i]) {    
                this->ft[4][i] = this->ftp1[2][i];  
            }
            if (this->barrier[5][i]) {    
                this->ft[5][i] = this->ftp1[7][i];  
            }
            if (this->barrier[6][i]) {    
                this->ft[6][i] = this->ftp1[8][i];  
            }
            if (this->barrier[7][i]) {    
                this->ft[7][i] = this->ftp1[5][i];  
            }
            if (this->barrier[8][i]) {    
                this->ft[8][i] = this->ftp1[6][i];  
            }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        for (int j = 0; j < this->ny; j++) {
            //.....xmin.....
            if (this->btxmin[j] == OUTLET) {
                this->ft[1][j] = this->ft[1][this->ny + j];
                this->ft[5][j] = this->ft[5][this->ny + j];
                this->ft[8][j] = this->ft[8][this->ny + j];
            } else if (this->btxmin[j] == BARRIER) {
                this->ft[1][j] = this->ftp1[3][j];
                this->ft[5][j] = this->ftp1[7][j];
                this->ft[8][j] = this->ftp1[6][j];
            } else if (this->btxmin[j] == MIRROR) {
                this->ft[1][j] = this->ftp1[3][j];
                this->ft[5][j] = this->ftp1[6][j];
                this->ft[8][j] = this->ftp1[7][j];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->ft[3][this->ny*(this->nx - 1) + j] = this->ft[3][this->ny*(this->nx - 2) + j];
                this->ft[6][this->ny*(this->nx - 1) + j] = this->ft[6][this->ny*(this->nx - 2) + j];
                this->ft[7][this->ny*(this->nx - 1) + j] = this->ft[7][this->ny*(this->nx - 2) + j];
            } else if (this->btxmax[j] == BARRIER) {
                this->ft[3][this->ny*(this->nx - 1) + j] = this->ftp1[1][this->ny*(this->nx - 1) + j];
                this->ft[6][this->ny*(this->nx - 1) + j] = this->ftp1[8][this->ny*(this->nx - 1) + j];
                this->ft[7][this->ny*(this->nx - 1) + j] = this->ftp1[5][this->ny*(this->nx - 1) + j];
            } else if (this->btxmax[j] == MIRROR) {
                this->ft[3][this->ny*(this->nx - 1) + j] = this->ftp1[1][this->ny*(this->nx - 1) + j];
                this->ft[6][this->ny*(this->nx - 1) + j] = this->ftp1[5][this->ny*(this->nx - 1) + j];
                this->ft[7][this->ny*(this->nx - 1) + j] = this->ftp1[8][this->ny*(this->nx - 1) + j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->ft[2][this->ny*i] = this->ft[2][this->ny*i + 1];
                this->ft[5][this->ny*i] = this->ft[5][this->ny*i + 1];
                this->ft[6][this->ny*i] = this->ft[6][this->ny*i + 1];
            } else if (this->btymin[i] == BARRIER) {
                this->ft[2][this->ny*i] = this->ftp1[4][this->ny*i];
                this->ft[5][this->ny*i] = this->ftp1[7][this->ny*i];
                this->ft[6][this->ny*i] = this->ftp1[8][this->ny*i];
            } else if (this->btymin[i] == MIRROR) {
                this->ft[2][this->ny*i] = this->ftp1[4][this->ny*i];
                this->ft[5][this->ny*i] = this->ftp1[8][this->ny*i];
                this->ft[6][this->ny*i] = this->ftp1[7][this->ny*i];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->ft[4][this->ny*(i + 1) - 1] = this->ft[4][this->ny*(i + 1) - 2];
                this->ft[7][this->ny*(i + 1) - 1] = this->ft[7][this->ny*(i + 1) - 2];
                this->ft[8][this->ny*(i + 1) - 1] = this->ft[8][this->ny*(i + 1) - 2];
            } else if (this->btymax[i] == BARRIER) {
                this->ft[4][this->ny*(i + 1) - 1] = this->ftp1[2][this->ny*(i + 1) - 1];
                this->ft[7][this->ny*(i + 1) - 1] = this->ftp1[5][this->ny*(i + 1) - 1];
                this->ft[8][this->ny*(i + 1) - 1] = this->ftp1[6][this->ny*(i + 1) - 1];
            } else if (this->btymax[i] == MIRROR) {
                this->ft[4][this->ny*(i + 1) - 1] = this->ftp1[2][this->ny*(i + 1) - 1];
                this->ft[7][this->ny*(i + 1) - 1] = this->ftp1[6][this->ny*(i + 1) - 1];
                this->ft[8][this->ny*(i + 1) - 1] = this->ftp1[5][this->ny*(i + 1) - 1];
            }
        }
    }


    template<class T>
    void D2Q9<T>::iStream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                for (int k = 0; k < D2Q9<T>::nc; k++) {
                    int ip1 = i + D2Q9<T>::cx[k] == -1 ? this->nx - 1 : (i + D2Q9<T>::cx[k] == this->nx ? 0 : i + D2Q9<T>::cx[k]);
                    int jp1 = j + D2Q9<T>::cy[k] == -1 ? this->ny - 1 : (j + D2Q9<T>::cy[k] == this->ny ? 0 : j + D2Q9<T>::cy[k]);
                    this->ft[k][this->ny*i + j] = this->ftp1[k][this->ny*ip1 + jp1];
                }
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->np; i++) {
            if (this->barrier[1][i]) {    
                this->ft[3][i] = this->ftp1[1][i];  
            }
            if (this->barrier[2][i]) {    
                this->ft[4][i] = this->ftp1[2][i];  
            }
            if (this->barrier[3][i]) {    
                this->ft[1][i] = this->ftp1[3][i];  
            }
            if (this->barrier[4][i]) {    
                this->ft[2][i] = this->ftp1[4][i];  
            }
            if (this->barrier[5][i]) {    
                this->ft[7][i] = this->ftp1[5][i];  
            }
            if (this->barrier[6][i]) {    
                this->ft[8][i] = this->ftp1[6][i];  
            }
            if (this->barrier[7][i]) {    
                this->ft[5][i] = this->ftp1[7][i];  
            }
            if (this->barrier[8][i]) {    
                this->ft[6][i] = this->ftp1[8][i];  
            }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        for (int j = 0; j < this->ny; j++) {
            //.....xmin.....
//  要修正
            if (this->btxmin[j] == OUTLET) {
                this->ft[3][j] = this->ft[3][this->ny + j];
                this->ft[6][j] = this->ft[6][this->ny + j];
                this->ft[7][j] = this->ft[7][this->ny + j];
            } else if (this->btxmin[j] == BARRIER) {
                this->ft[3][j] = this->ftp1[1][j];
                this->ft[6][j] = this->ftp1[8][j];
                this->ft[7][j] = this->ftp1[5][j];
            } else if (this->btxmin[j] == MIRROR) {
                this->ft[3][j] = this->ftp1[1][j];
                this->ft[6][j] = this->ftp1[5][j];
                this->ft[7][j] = this->ftp1[8][j];
            }

            //.....xmax.....
//  要修正
            if (this->btxmax[j] == OUTLET) {
                this->ft[1][this->ny*(this->nx - 1) + j] = this->ft[1][this->ny*(this->nx - 2) + j];
                this->ft[5][this->ny*(this->nx - 1) + j] = this->ft[5][this->ny*(this->nx - 2) + j];
                this->ft[8][this->ny*(this->nx - 1) + j] = this->ft[8][this->ny*(this->nx - 2) + j];
            } else if (this->btxmax[j] == BARRIER) {
                this->ft[1][this->ny*(this->nx - 1) + j] = this->ftp1[3][this->ny*(this->nx - 1) + j];
                this->ft[5][this->ny*(this->nx - 1) + j] = this->ftp1[7][this->ny*(this->nx - 1) + j];
                this->ft[8][this->ny*(this->nx - 1) + j] = this->ftp1[6][this->ny*(this->nx - 1) + j];
            } else if (this->btxmax[j] == MIRROR) {
                this->ft[1][this->ny*(this->nx - 1) + j] = this->ftp1[3][this->ny*(this->nx - 1) + j];
                this->ft[5][this->ny*(this->nx - 1) + j] = this->ftp1[6][this->ny*(this->nx - 1) + j];
                this->ft[8][this->ny*(this->nx - 1) + j] = this->ftp1[7][this->ny*(this->nx - 1) + j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
//  要修正
            if (this->btymin[i] == OUTLET) {
                this->ft[4][this->ny*i] = this->ft[4][this->ny*i + 1];
                this->ft[7][this->ny*i] = this->ft[7][this->ny*i + 1];
                this->ft[8][this->ny*i] = this->ft[8][this->ny*i + 1];
            } else if (this->btymin[i] == BARRIER) {
                this->ft[4][this->ny*i] = this->ftp1[2][this->ny*i];
                this->ft[7][this->ny*i] = this->ftp1[5][this->ny*i];
                this->ft[8][this->ny*i] = this->ftp1[6][this->ny*i];
            } else if (this->btymin[i] == MIRROR) {
                this->ft[4][this->ny*i] = this->ftp1[2][this->ny*i];
                this->ft[7][this->ny*i] = this->ftp1[6][this->ny*i];
                this->ft[8][this->ny*i] = this->ftp1[5][this->ny*i];
            }

            //.....ymax.....
//  要修正
            if (this->btymax[i] == OUTLET) {
                this->ft[2][this->ny*(i + 1) - 1] = this->ft[2][this->ny*(i + 1) - 2];
                this->ft[5][this->ny*(i + 1) - 1] = this->ft[5][this->ny*(i + 1) - 2];
                this->ft[6][this->ny*(i + 1) - 1] = this->ft[6][this->ny*(i + 1) - 2];
            } else if (this->btymax[i] == BARRIER) {
                this->ft[2][this->ny*(i + 1) - 1] = this->ftp1[4][this->ny*(i + 1) - 1];
                this->ft[5][this->ny*(i + 1) - 1] = this->ftp1[7][this->ny*(i + 1) - 1];
                this->ft[6][this->ny*(i + 1) - 1] = this->ftp1[8][this->ny*(i + 1) - 1];
            } else if (this->btymax[i] == MIRROR) {
                this->ft[2][this->ny*(i + 1) - 1] = this->ftp1[4][this->ny*(i + 1) - 1];
                this->ft[5][this->ny*(i + 1) - 1] = this->ftp1[8][this->ny*(i + 1) - 1];
                this->ft[6][this->ny*(i + 1) - 1] = this->ftp1[7][this->ny*(i + 1) - 1];
            }
        }
    }


    template<class T>
    bool D2Q9<T>::GetBarrier(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->barrier[0][this->ny*_i + _j];
    }


    template<class T>
    BOUNDARYTYPE D2Q9<T>::GetBoundary(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        if (_i == 0) {
            return this->btxmin[_j];
        } else if (_i == this->nx - 1) {
            return this->btxmax[_j];
        } else if (_j == 0) {
            return this->btymin[_i];
        } else if (_j == this->ny - 1) {
            return this->btymax[_i];
        } else {
            if (this->barrier[0][this->ny*_i + _j]) {
                return BARRIER;
            } else {
                return OTHER;
            }
        }
    }


    template<class T>
    int D2Q9<T>::GetIndex(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->ny*_i + _j;
    }
}