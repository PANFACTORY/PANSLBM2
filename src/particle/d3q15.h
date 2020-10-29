//*****************************************************************************
//  Title       :   src/particle/d3q15.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/29
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>


namespace PANSLBM2 {
    enum BOUNDARYTYPE {
        PERIODIC, OUTLET, BARRIER, MIRROR, OTHER,
    };

    template<class T>
    class D3Q15 {
public:
        D3Q15() = delete;
        D3Q15(int _nx, int _ny, int _nz);
        D3Q15(const D3Q15<T>& _p);
        ~D3Q15();

        void SetBarrier(int _i, int _j, int _k, bool _isbarrier);
        void SetBoundary(int _i, int _j, int _k, BOUNDARYTYPE _boundarytype);
        void SetRho(int _i, int _j, int _k, T _rho, T _u0, T _u1);      //  Set boundary condition for NavierStokes  
        void SetU(int _i, int _j, int _k, T _ux, T _uy, T _uz);
        void SetTemperature(int _i, int _j, int _k, int _temperature);  //  Set boundary condition for Advection
        void SetFlux(int _i, int _j, int _k, T _ux, T _uy, T _uz, T _q);
        void SetiRho(int _i, int _j, int _k);                           //  Set boundary condition for Adjoint of NavierStokes  
        void SetiU(int _i, int _j, int _k, T _ux, T _uy, T _uz);

        void Stream();
        void iStream();
        
        bool GetBarrier(int _i, int _j, int _k) const;
        BOUNDARYTYPE GetBoundary(int _i, int _j, int _k) const;
        int GetIndex(int _i, int _j, int _k);
        
        const int nx, ny, nz, np;                                       //  nx&ny : number of points along x&y coordinate, np : number of all points
        static const int nc = 15, nd = 3, cx[nc], cy[nc], cz[nc];       //  nc : number of particles, nd : number of dimensions
        static const T ei[nc];

        T dx, dt;
        T *ft[nc], *ftp1[nc];
        bool *barrier[nc];
        BOUNDARYTYPE *btxmin, *btxmax, *btymin, *btymax, *btzmin, *btzmax;
    };


    template<class T>
    const int D3Q15<T>::cx[D3Q15<T>::nc] = { 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, -1, -1 };


    template<class T>
    const int D3Q15<T>::cy[D3Q15<T>::nc] = { 0, 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, -1, -1, 1, -1 };


    template<class T>
    const int D3Q15<T>::cz[D3Q15<T>::nc] = { 0, 0, 0, 1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1 };


    template<class T>
    const T D3Q15<T>::ei[D3Q15<T>::nc] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0 };


    template<class T>
    D3Q15<T>::D3Q15(int _nx, int _ny, int _nz) : nx(_nx), ny(_ny), nz(_nz), np(_nx*_ny*_nz) {
        assert(0 < _nx && 0 < _ny && 0 < _nz);
        this->dx = 1.0;
        this->dt = 1.0;

        for (int i = 0; i < D3Q15<T>::nc; i++) {
            this->ft[i] = new T[this->np];
            this->ftp1[i] = new T[this->np];
            this->barrier[i] = new bool[this->np];

            for (int j = 0; j < this->np; j++) {
                this->ft[i][j] = D3Q15<T>::ei[i];
                this->ftp1[i][j] = T();
                this->barrier[i][j] = false;
            }    
        }

        this->btxmin = new BOUNDARYTYPE[this->ny*this->nz];
        this->btxmax = new BOUNDARYTYPE[this->ny*this->nz];
        this->btymin = new BOUNDARYTYPE[this->nz*this->nx];
        this->btymax = new BOUNDARYTYPE[this->nz*this->nx];
        this->btzmin = new BOUNDARYTYPE[this->nx*this->ny];
        this->btzmax = new BOUNDARYTYPE[this->nx*this->ny];

        for (int i = 0; i < this->ny*this->nz; i++) {
            this->btxmin[i] = PERIODIC;
            this->btxmax[i] = PERIODIC;
        }

        for (int i = 0; i < this->nz*this->nx; i++) {
            this->btymin[i] = PERIODIC;
            this->btymax[i] = PERIODIC;
        }

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->btzmin[i] = PERIODIC;
            this->btzmax[i] = PERIODIC;
        }
    }


    template<class T>
    D3Q15<T>::D3Q15(const D3Q15<T>& _p) : nx(_p.nx), ny(_p.ny), nz(_p.nz), np(_p.np) {
        this->dx = _p.dx;
        this->dt = _p.dt;
        
        for (int i = 0; i < D3Q15<T>::nc; i++) {
            this->ft[i] = new T[this->np];
            this->ftp1[i] = new T[this->np];
            this->barrier[i] = new bool[this->np];

            for (int j = 0; j < this->np; j++) {
                this->ft[i][j] = _p.ft[i][j];
                this->ftp1[i][j] = _p.ftp1[i][j];
                this->barrier[i][j] = _p.barrier[i][j];
            }    
        }

        this->btxmin = new BOUNDARYTYPE[this->ny*this->nz];
        this->btxmax = new BOUNDARYTYPE[this->ny*this->nz];
        this->btymin = new BOUNDARYTYPE[this->nz*this->nx];
        this->btymax = new BOUNDARYTYPE[this->nz*this->nx];
        this->btzmin = new BOUNDARYTYPE[this->nx*this->ny];
        this->btzmax = new BOUNDARYTYPE[this->nx*this->ny];

        for (int i = 0; i < this->ny*this->nz; i++) {
            this->btxmin[i] = _p.btxmin[i];
            this->btxmax[i] = _p.btxmax[i];
        }

        for (int i = 0; i < this->nz*this->nx; i++) {
            this->btymin[i] = _p.btymin[i];
            this->btymax[i] = _p.btymax[i];
        }

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->btzmin[i] = _p.btzmin[i];
            this->btzmax[i] = _p.btzmax[i];
        }
    }


    template<class T>
    D3Q15<T>::~D3Q15() {
        for (int i = 0; i < D3Q15<T>::nc; i++) {
            delete[] this->ft[i];
            delete[] this->ftp1[i];
            delete[] this->barrier[i];
        }
        
        delete[] this->btxmin;
        delete[] this->btxmax;
        delete[] this->btymin;
        delete[] this->btymax;
        delete[] this->btzmin;
        delete[] this->btzmax;
    }


    template<class T>
    void D3Q15<T>::SetBarrier(int _i, int _j, int _k, bool _isbarrier) {
        for (int i = 0; i < D3Q15<T>::nc; i++) {
            this->barrier[i][(_i + D3Q15<T>::cx[i]) + this->nx*(_j + D3Q15<T>::cy[i]) + this->nx*this->ny*(_k + D3Q15<T>::cz[i])] = _isbarrier;
        }
    }


    template<class T>
    void D3Q15<T>::SetBoundary(int _i, int _j, int _k, BOUNDARYTYPE _boundarytype) {
        if (_i == 0) {
            this->btxmin[_j + this->ny*_k] = _boundarytype;
        } 
        if (_i == this->nx - 1) {
            this->btxmax[_j + this->ny*_k] = _boundarytype;
        } 
        if (_j == 0) {
            this->btymin[_k + this->nz*_i] = _boundarytype;
        } 
        if (_j == this->ny - 1) {
            this->btymax[_k + this->nz*_i] = _boundarytype;
        }
        if (_k == 0) {
            this->btzmin[_i + this->nx*_j] = _boundarytype;
        } 
        if (_k == this->nz - 1) {
            this->btzmax[_i + this->nx*_j] = _boundarytype;
        } 
        if (_i != 0 && _i != this->nx - 1 && _j != 0 && _j != this->ny - 1 && _k != 0 && _k != this->nz - 1) {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D3Q15<T>::SetRho(int _i, int _j, T _rho, T _u) {
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
    void D3Q15<T>::SetU(int _i, int _j, T _ux, T _uy) {
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
    void D3Q15<T>::SetTemperature(int _i, int _j, int _k, int _temperature) {
        int ijk = _i + this->nx*_j + this->nx*this->ny*_k;
        if (_i == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[2][ijk] - this->ft[3][ijk] - this->ft[4][ijk] - this->ft[5][ijk] - this->ft[6][ijk] - this->ft[8][ijk] - this->ft[11][ijk] - this->ft[13][ijk] - this->ft[14][ijk]);
            this->ft[1][ijk] = temperature0/9.0;
            this->ft[7][ijk] = temperature0/72.0;
            this->ft[9][ijk] = temperature0/72.0;
            this->ft[10][ijk] = temperature0/72.0;
            this->ft[12][ijk] = temperature0/72.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[1][ijk] - this->ft[2][ijk] - this->ft[3][ijk] - this->ft[5][ijk] - this->ft[6][ijk] - this->ft[7][ijk] - this->ft[9][ijk] - this->ft[10][ijk] - this->ft[12][ijk]);
            this->ft[4][ijk] = temperature0/9.0;
            this->ft[8][ijk] = temperature0/72.0;
            this->ft[11][ijk] = temperature0/72.0;
            this->ft[13][ijk] = temperature0/72.0;
            this->ft[14][ijk] = temperature0/72.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[1][ijk] - this->ft[3][ijk] - this->ft[4][ijk] - this->ft[5][ijk] - this->ft[6][ijk] - this->ft[9][ijk] - this->ft[11][ijk] - this->ft[12][ijk] - this->ft[14][ijk]);
            this->ft[2][ijk] = temperature0/9.0;
            this->ft[7][ijk] = temperature0/72.0;
            this->ft[8][ijk] = temperature0/72.0;
            this->ft[10][ijk] = temperature0/72.0;
            this->ft[13][ijk] = temperature0/72.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[1][ijk] - this->ft[2][ijk] - this->ft[3][ijk] - this->ft[4][ijk] - this->ft[6][ijk] - this->ft[7][ijk] - this->ft[8][ijk] - this->ft[10][ijk] - this->ft[13][ijk]);
            this->ft[5][ijk] = temperature0/9.0;
            this->ft[9][ijk] = temperature0/72.0;
            this->ft[11][ijk] = temperature0/72.0;
            this->ft[12][ijk] = temperature0/72.0;
            this->ft[14][ijk] = temperature0/72.0;
        } else if (_k == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[1][ijk] - this->ft[2][ijk] - this->ft[4][ijk] - this->ft[5][ijk] - this->ft[6][ijk] - this->ft[10][ijk] - this->ft[11][ijk] - this->ft[12][ijk] - this->ft[13][ijk]);
            this->ft[3][ijk] = temperature0/9.0;
            this->ft[7][ijk] = temperature0/72.0;
            this->ft[8][ijk] = temperature0/72.0;
            this->ft[9][ijk] = temperature0/72.0;
            this->ft[14][ijk] = temperature0/72.0;
        } else if (_k == this->nz - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ijk] - this->ft[1][ijk] - this->ft[2][ijk] - this->ft[3][ijk] - this->ft[4][ijk] - this->ft[5][ijk] - this->ft[7][ijk] - this->ft[8][ijk] - this->ft[9][ijk] - this->ft[14][ijk]);
            this->ft[6][ijk] = temperature0/9.0;
            this->ft[10][ijk] = temperature0/72.0;
            this->ft[11][ijk] = temperature0/72.0;
            this->ft[12][ijk] = temperature0/72.0;
            this->ft[13][ijk] = temperature0/72.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }
    
    
    template<class T>
    void D3Q15<T>::SetFlux(int _i, int _j, int _k, T _ux, T _uy, T _uz, T _q) {
        int ijk = _i + this->nx*_j + this->nx*this->ny*_k;
        if (_i == 0) {
            T temperature0 = 6.0*(_q + this->ft[4][ijk] + this->ft[8][ijk] + this->ft[11][ijk] + this->ft[13][ijk] + this->ft[14][ijk])/(1.0 - 3.0*_ux);
            this->ft[1][ijk] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[7][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[9][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[10][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[12][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(-_q + this->ft[1][ijk] + this->ft[7][ijk] + this->ft[9][ijk] + this->ft[10][ijk] + this->ft[12][ijk])/(1.0 + 3.0*_ux);
            this->ft[4][ijk] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[8][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[11][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[13][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[14][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_q + this->ft[5][ijk] + this->ft[9][ijk] + this->ft[11][ijk] + this->ft[12][ijk] + this->ft[14][ijk])/(1.0 - 3.0*_uy);
            this->ft[2][ijk] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[7][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[8][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[10][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[13][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(-_q + this->ft[2][ijk] + this->ft[7][ijk] + this->ft[8][ijk] + this->ft[10][ijk] + this->ft[13][ijk])/(1.0 + 3.0*_uy);
            this->ft[5][ijk] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[9][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[11][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[12][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[14][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
        } else if (_k == 0) {
            T temperature0 = 6.0*(_q + this->ft[6][ijk] + this->ft[10][ijk] + this->ft[11][ijk] + this->ft[12][ijk] + this->ft[13][ijk])/(1.0 - 3.0*_uz);
            this->ft[3][ijk] = temperature0*(1.0 + 3.0*_uz)/9.0;
            this->ft[7][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[8][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[9][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
            this->ft[14][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy + 3.0*_uz)/72.0;
        } else if (_k == this->nz - 1) {
            T temperature0 = 6.0*(-_q + this->ft[3][ijk] + this->ft[7][ijk] + this->ft[8][ijk] + this->ft[9][ijk] + this->ft[14][ijk])/(1.0 + 3.0*_uz);
            this->ft[6][ijk] = temperature0*(1.0 - 3.0*_uz)/9.0;
            this->ft[10][ijk] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[11][ijk] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[12][ijk] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy - 3.0*_uz)/72.0;
            this->ft[13][ijk] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy - 3.0*_uz)/72.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D3Q15<T>::SetiRho(int _i, int _j) {
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
    void D3Q15<T>::SetiU(int _i, int _j, T _ux, T _uy) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T rho0 = (-2.0 + _ux*(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]) + 3.0*_uy*(this->ft[5][ij] - this->ft[8][ij]))/(3.0*(1.0 - _ux));
            this->ft[3][ij] = this->ft[1][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
        } else if (_i == this->nx - 1) {          
            T rho0 = (-2.0 - _ux*(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]) + 3.0*_uy*(this->ft[6][ij] - this->ft[7][ij]))/(3.0*(1.0 + _ux));
            this->ft[1][ij] = this->ft[3][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == 0) {
            T rho0 = (-2.0 + _uy*(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]) + 3.0*_ux*(this->ft[5][ij] - this->ft[6][ij]))/(3.0*(1.0 - _uy));
            this->ft[4][ij] = this->ft[2][ij] + rho0;
            this->ft[7][ij] = this->ft[5][ij] + rho0;
            this->ft[8][ij] = this->ft[6][ij] + rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = (-2.0 - _uy*(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]) + 3.0*_ux*(this->ft[8][ij] - this->ft[7][ij]))/(3.0*(1.0 + _uy));
            this->ft[2][ij] = this->ft[4][ij] + rho0;
            this->ft[5][ij] = this->ft[7][ij] + rho0;
            this->ft[6][ij] = this->ft[8][ij] + rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D3Q15<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                for (int k = 0; k < this->nz; k++) {
                    for (int l = 0; l < D3Q15<T>::nc; l++) {
                        int ip1 = i - D3Q15<T>::cx[l] == -1 ? this->nx - 1 : (i - D3Q15<T>::cx[l] == this->nx ? 0 : i - D3Q15<T>::cx[l]);
                        int jp1 = j - D3Q15<T>::cy[l] == -1 ? this->ny - 1 : (j - D3Q15<T>::cy[l] == this->ny ? 0 : j - D3Q15<T>::cy[l]);
                        int kp1 = k - D3Q15<T>::cz[l] == -1 ? this->nz - 1 : (k - D3Q15<T>::cz[l] == this->nz ? 0 : k - D3Q15<T>::cz[l]);
                        this->ft[l][i + this->nx*j + this->nx*this->ny*k] = this->ftp1[k][ip1 + this->nx*jp1 + this->nx*this->ny*kp1];
                    }
                }
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->np; i++) {
            if (this->barrier[1][i])    {   this->ft[1][i] = this->ftp1[4][i];      }
            if (this->barrier[2][i])    {   this->ft[2][i] = this->ftp1[5][i];      }
            if (this->barrier[3][i])    {   this->ft[3][i] = this->ftp1[6][i];      }
            if (this->barrier[4][i])    {   this->ft[4][i] = this->ftp1[1][i];      }
            if (this->barrier[5][i])    {   this->ft[5][i] = this->ftp1[2][i];      }
            if (this->barrier[6][i])    {   this->ft[6][i] = this->ftp1[3][i];      }
            if (this->barrier[7][i])    {   this->ft[7][i] = this->ftp1[11][i];     }
            if (this->barrier[8][i])    {   this->ft[8][i] = this->ftp1[12][i];     }
            if (this->barrier[9][i])    {   this->ft[9][i] = this->ftp1[13][i];     }
            if (this->barrier[10][i])   {   this->ft[10][i] = this->ftp1[14][i];    }
            if (this->barrier[11][i])   {   this->ft[11][i] = this->ftp1[7][i];     }
            if (this->barrier[12][i])   {   this->ft[12][i] = this->ftp1[8][i];     }
            if (this->barrier[13][i])   {   this->ft[13][i] = this->ftp1[9][i];     }
            if (this->barrier[14][i])   {   this->ft[14][i] = this->ftp1[10][i];    }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        //  xmin and xmax
        for (int j = 0; j < this->ny; j++) {
            for (int k = 0; k < this->nz; k++) {
                int jk = j + this->ny*k;

                //.....xmin.....
                int minjk = 0 + this->nx*j + this->nx*this->ny*k;
                int minp1jk = 1 + this->nx*j + this->nx*this->ny*k;
                if (this->btxmin[jk] == OUTLET) {
                    this->ft[1][minjk] = this->ft[1][minp1jk];
                    this->ft[7][minjk] = this->ft[7][minp1jk];
                    this->ft[9][minjk] = this->ft[9][minp1jk];
                    this->ft[10][minjk] = this->ft[10][minp1jk];
                    this->ft[12][minjk] = this->ft[12][minp1jk];
                } else if (this->btxmin[jk] == BARRIER) {
                    this->ft[1][minjk] = this->ft[4][minjk];
                    this->ft[7][minjk] = this->ft[11][minjk];
                    this->ft[9][minjk] = this->ft[13][minjk];
                    this->ft[10][minjk] = this->ft[14][minjk];
                    this->ft[12][minjk] = this->ft[8][minjk];
                } else if (this->btxmin[jk] == MIRROR) {
                    this->ft[1][minjk] = this->ft[4][minjk];
                    this->ft[7][minjk] = this->ft[8][minjk];
                    this->ft[9][minjk] = this->ft[14][minjk];
                    this->ft[10][minjk] = this->ft[13][minjk];
                    this->ft[12][minjk] = this->ft[11][minjk];
                }

                //.....xmax.....
                int maxjk = (this->nx - 1) + this->nx*j + this->nx*this->ny*k;
                int maxm1jk = (this->nx - 2) + this->nx*j + this->nx*this->ny*k;
                if (this->btxmax[jk] == OUTLET) {
                    this->ft[4][maxjk] = this->ft[4][maxm1jk];
                    this->ft[8][maxjk] = this->ft[8][maxm1jk];
                    this->ft[11][maxjk] = this->ft[11][maxm1jk];
                    this->ft[13][maxjk] = this->ft[13][maxm1jk];
                    this->ft[14][maxjk] = this->ft[14][maxm1jk];
                } else if (this->btxmax[jk] == BARRIER) {
                    this->ft[4][maxjk] = this->ft[1][maxjk];
                    this->ft[8][maxjk] = this->ft[12][maxjk];
                    this->ft[11][maxjk] = this->ft[7][maxjk];
                    this->ft[13][maxjk] = this->ft[9][maxjk];
                    this->ft[14][maxjk] = this->ft[10][maxjk];
                } else if (this->btxmax[jk] == MIRROR) {
                    this->ft[4][maxjk] = this->ft[1][maxjk];
                    this->ft[8][maxjk] = this->ft[7][maxjk];
                    this->ft[11][maxjk] = this->ft[12][maxjk];
                    this->ft[13][maxjk] = this->ft[10][maxjk];
                    this->ft[14][maxjk] = this->ft[9][maxjk];
                }
            }
        }

        //  ymin and ymax
        for (int k = 0; k < this->nz; k++) {
            for (int i = 0; i < this->nx; i++) {
                int ki = k + this->nz*i;

                //.....ymin.....
                int imink = i + this->nx*0 + this->nx*this->ny*k;
                int iminp1k = i + this->nx*1 + this->nx*this->ny*k;
                if (this->btymin[ki] == OUTLET) {
                    this->ft[2][imink] = this->ft[2][iminp1k];
                    this->ft[7][imink] = this->ft[7][iminp1k];
                    this->ft[8][imink] = this->ft[8][iminp1k];
                    this->ft[10][imink] = this->ft[10][iminp1k];
                    this->ft[13][imink] = this->ft[13][iminp1k];
                } else if (this->btymin[ki] == BARRIER) {
                    this->ft[2][imink] = this->ft[5][imink];
                    this->ft[7][imink] = this->ft[11][imink];
                    this->ft[8][imink] = this->ft[12][imink];
                    this->ft[10][imink] = this->ft[14][imink];
                    this->ft[13][imink] = this->ft[9][imink];
                } else if (this->btymin[ki] == MIRROR) {
                    this->ft[2][imink] = this->ft[5][imink];
                    this->ft[7][imink] = this->ft[9][imink];
                    this->ft[8][imink] = this->ft[14][imink];
                    this->ft[10][imink] = this->ft[12][imink];
                    this->ft[13][imink] = this->ft[11][imink];
                }

                //.....ymax.....
                int imaxk = i + this->nx*(this->ny - 1) + this->nx*this->ny*k;
                int imaxm1k = i + this->nx*(this->ny - 2) + this->nx*this->ny*k;
                if (this->btymax[ki] == OUTLET) {
                    this->ft[5][imaxk] = this->ft[5][imaxm1k];
                    this->ft[9][imaxk] = this->ft[9][imaxm1k];
                    this->ft[11][imaxk] = this->ft[11][imaxm1k];
                    this->ft[12][imaxk] = this->ft[12][imaxm1k];
                    this->ft[14][imaxk] = this->ft[14][imaxm1k];
                } else if (this->btymax[ki] == BARRIER) {
                    this->ft[5][imaxk] = this->ft[2][imaxk];
                    this->ft[9][imaxk] = this->ft[13][imaxk];
                    this->ft[11][imaxk] = this->ft[7][imaxk];
                    this->ft[12][imaxk] = this->ft[8][imaxk];
                    this->ft[14][imaxk] = this->ft[10][imaxk];
                } else if (this->btymax[ki] == MIRROR) {
                    this->ft[5][imaxk] = this->ft[2][imaxk];
                    this->ft[9][imaxk] = this->ft[7][imaxk];
                    this->ft[11][imaxk] = this->ft[13][imaxk];
                    this->ft[12][imaxk] = this->ft[10][imaxk];
                    this->ft[14][imaxk] = this->ft[8][imaxk];
                }
            }
        }

        //  zmin and zmax
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; i++) {
                int ij = i + this->nx*j;

                //.....zmin.....
                int ijmin = i + this->nx*j + this->nx*this->ny*0;
                int ijminp1 = i + this->nx*j + this->nx*this->ny*1;
                if (this->btzmin[ki] == OUTLET) {
                    this->ft[3][ijmin] = this->ft[3][ijminp1];
                    this->ft[7][ijmin] = this->ft[7][ijminp1];
                    this->ft[8][ijmin] = this->ft[8][ijminp1];
                    this->ft[9][ijmin] = this->ft[9][ijminp1];
                    this->ft[14][ijmin] = this->ft[14][ijminp1];
                } else if (this->btzmin[ki] == BARRIER) {
                    this->ft[3][ijmin] = this->ft[6][ijmin];
                    this->ft[7][ijmin] = this->ft[11][ijmin];
                    this->ft[8][ijmin] = this->ft[12][ijmin];
                    this->ft[9][ijmin] = this->ft[13][ijmin];
                    this->ft[14][ijmin] = this->ft[10][ijmin];
                } else if (this->btzmin[ki] == MIRROR) {
                    this->ft[3][ijmin] = this->ft[6][ijmin];
                    this->ft[7][ijmin] = this->ft[10][ijmin];
                    this->ft[8][ijmin] = this->ft[13][ijmin];
                    this->ft[9][ijmin] = this->ft[12][ijmin];
                    this->ft[14][ijmin] = this->ft[11][ijmin];
                }

                //.....zmax.....
                int ijmax = i + this->nx*j + this->nx*this->ny*(this->nz - 1);
                int ijmaxm1 = i + this->nx*j + this->nx*this->ny*(this->nz - 2);
                if (this->btzmax[ki] == OUTLET) {
                    this->ft[6][ijmax] = this->ft[6][ijmaxm1];
                    this->ft[10][ijmax] = this->ft[10][ijmaxm1];
                    this->ft[11][ijmax] = this->ft[11][ijmaxm1];
                    this->ft[12][ijmax] = this->ft[12][ijmaxm1];
                    this->ft[13][ijmax] = this->ft[13][ijmaxm1];
                } else if (this->btzmax[ki] == BARRIER) {
                    this->ft[6][ijmax] = this->ft[3][ijmax];
                    this->ft[10][ijmax] = this->ft[14][ijmax];
                    this->ft[11][ijmax] = this->ft[7][ijmax];
                    this->ft[12][ijmax] = this->ft[8][ijmax];
                    this->ft[13][ijmax] = this->ft[9][ijmax];
                } else if (this->btzmax[ki] == MIRROR) {
                    this->ft[6][ijmax] = this->ft[3][ijmax];
                    this->ft[10][ijmax] = this->ft[7][ijmax];
                    this->ft[11][ijmax] = this->ft[14][ijmax];
                    this->ft[12][ijmax] = this->ft[9][ijmax];
                    this->ft[13][ijmax] = this->ft[8][ijmax];
                }
            }
        }
    }


    template<class T>
    void D3Q15<T>::iStream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                for (int k = 0; k < this->nz; k++) {
                    for (int l = 0; l < D3Q15<T>::nc; l++) {
                        int ip1 = i + D3Q15<T>::cx[l] == -1 ? this->nx - 1 : (i + D3Q15<T>::cx[l] == this->nx ? 0 : i + D3Q15<T>::cx[l]);
                        int jp1 = j + D3Q15<T>::cy[l] == -1 ? this->ny - 1 : (j + D3Q15<T>::cy[l] == this->ny ? 0 : j + D3Q15<T>::cy[l]);
                        int kp1 = k + D3Q15<T>::cz[l] == -1 ? this->nz - 1 : (k + D3Q15<T>::cz[l] == this->nz ? 0 : k + D3Q15<T>::cz[l]);
                        this->ft[l][i + this->nx*j + this->nx*this->ny*k] = this->ftp1[k][ip1 + this->nx*jp1 + this->nx*this->ny*kp1];
                    }
                }
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->np; i++) {
            if (this->barrier[1][i])    {   this->ft[4][i] = this->ftp1[1][i];      }
            if (this->barrier[2][i])    {   this->ft[5][i] = this->ftp1[2][i];      }
            if (this->barrier[3][i])    {   this->ft[6][i] = this->ftp1[3][i];      }
            if (this->barrier[4][i])    {   this->ft[1][i] = this->ftp1[4][i];      }
            if (this->barrier[5][i])    {   this->ft[2][i] = this->ftp1[5][i];      }
            if (this->barrier[6][i])    {   this->ft[3][i] = this->ftp1[6][i];      }
            if (this->barrier[7][i])    {   this->ft[11][i] = this->ftp1[7][i];     }
            if (this->barrier[8][i])    {   this->ft[12][i] = this->ftp1[8][i];     }
            if (this->barrier[9][i])    {   this->ft[13][i] = this->ftp1[9][i];     }
            if (this->barrier[10][i])   {   this->ft[14][i] = this->ftp1[10][i];    }
            if (this->barrier[11][i])   {   this->ft[7][i] = this->ftp1[11][i];     }
            if (this->barrier[12][i])   {   this->ft[8][i] = this->ftp1[12][i];     }
            if (this->barrier[13][i])   {   this->ft[9][i] = this->ftp1[13][i];     }
            if (this->barrier[14][i])   {   this->ft[10][i] = this->ftp1[14][i];    }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        //  xmin and xmax
        for (int j = 0; j < this->ny; j++) {
            for (int k = 0; k < this->nz; k++) {
                int jk = j + this->ny*k;

                //.....xmin.....
                int minjk = 0 + this->nx*j + this->nx*this->ny*k;
                int minp1jk = 1 + this->nx*j + this->nx*this->ny*k;
                if (this->btxmin[jk] == OUTLET) {
                    this->ft[4][minjk] = this->ft[4][minp1jk];
                    this->ft[8][minjk] = this->ft[8][minp1jk];
                    this->ft[11][minjk] = this->ft[11][minp1jk];
                    this->ft[13][minjk] = this->ft[13][minp1jk];
                    this->ft[14][minjk] = this->ft[14][minp1jk];
                } else if (this->btxmin[jk] == BARRIER) {
                    this->ft[4][minjk] = this->ft[1][minjk];
                    this->ft[8][minjk] = this->ft[12][minjk];
                    this->ft[11][minjk] = this->ft[7][minjk];
                    this->ft[13][minjk] = this->ft[9][minjk];
                    this->ft[14][minjk] = this->ft[10][minjk];
                } else if (this->btxmin[jk] == MIRROR) {
                    this->ft[4][minjk] = this->ft[1][minjk];
                    this->ft[8][minjk] = this->ft[7][minjk];
                    this->ft[11][minjk] = this->ft[12][minjk];
                    this->ft[13][minjk] = this->ft[10][minjk];
                    this->ft[14][minjk] = this->ft[9][minjk];
                }

                //.....xmax.....
                int maxjk = (this->nx - 1) + this->nx*j + this->nx*this->ny*k;
                int maxm1jk = (this->nx - 2) + this->nx*j + this->nx*this->ny*k;
                if (this->btxmax[jk] == OUTLET) {
                    this->ft[1][maxjk] = this->ft[1][maxm1jk];
                    this->ft[7][maxjk] = this->ft[7][maxm1jk];
                    this->ft[9][maxjk] = this->ft[9][maxm1jk];
                    this->ft[10][maxjk] = this->ft[10][maxm1jk];
                    this->ft[12][maxjk] = this->ft[12][maxm1jk];
                } else if (this->btxmax[jk] == BARRIER) {
                    this->ft[1][maxjk] = this->ft[4][maxjk];
                    this->ft[7][maxjk] = this->ft[11][maxjk];
                    this->ft[9][maxjk] = this->ft[13][maxjk];
                    this->ft[10][maxjk] = this->ft[14][maxjk];
                    this->ft[12][maxjk] = this->ft[8][maxjk];
                } else if (this->btxmax[jk] == MIRROR) {
                    this->ft[1][maxjk] = this->ft[4][maxjk];
                    this->ft[7][maxjk] = this->ft[8][maxjk];
                    this->ft[9][maxjk] = this->ft[14][maxjk];
                    this->ft[10][maxjk] = this->ft[13][maxjk];
                    this->ft[12][maxjk] = this->ft[11][maxjk];
                }
            }
        }

        //  ymin and ymax
        for (int k = 0; k < this->nz; k++) {
            for (int i = 0; i < this->nx; i++) {
                int ki = k + this->nz*i;

                //.....ymin.....
                int imink = i + this->nx*0 + this->nx*this->ny*k;
                int iminp1k = i + this->nx*1 + this->nx*this->ny*k;
                if (this->btymin[ki] == OUTLET) {
                    this->ft[5][imink] = this->ft[5][iminp1k];
                    this->ft[9][imink] = this->ft[9][iminp1k];
                    this->ft[11][imink] = this->ft[11][iminp1k];
                    this->ft[12][imink] = this->ft[12][iminp1k];
                    this->ft[14][imink] = this->ft[14][iminp1k];
                } else if (this->btymin[ki] == BARRIER) {
                    this->ft[5][imink] = this->ft[2][imink];
                    this->ft[9][imink] = this->ft[13][imink];
                    this->ft[11][imink] = this->ft[7][imink];
                    this->ft[12][imink] = this->ft[8][imink];
                    this->ft[14][imink] = this->ft[10][imink];
                } else if (this->btymin[ki] == MIRROR) {
                    this->ft[5][imink] = this->ft[2][imink];
                    this->ft[9][imink] = this->ft[7][imink];
                    this->ft[11][imink] = this->ft[13][imink];
                    this->ft[12][imink] = this->ft[10][imink];
                    this->ft[14][imink] = this->ft[8][imink];
                }

                //.....ymax.....
                int imaxk = i + this->nx*(this->ny - 1) + this->nx*this->ny*k;
                int imaxm1k = i + this->nx*(this->ny - 2) + this->nx*this->ny*k;
                if (this->btymax[ki] == OUTLET) {
                    this->ft[2][imaxk] = this->ft[2][imaxm1k];
                    this->ft[7][imaxk] = this->ft[7][imaxm1k];
                    this->ft[8][imaxk] = this->ft[8][imaxm1k];
                    this->ft[10][imaxk] = this->ft[10][imaxm1k];
                    this->ft[13][imaxk] = this->ft[13][imaxm1k];
                } else if (this->btymax[ki] == BARRIER) {
                    this->ft[2][imaxk] = this->ft[5][imaxk];
                    this->ft[7][imaxk] = this->ft[11][imaxk];
                    this->ft[8][imaxk] = this->ft[12][imaxk];
                    this->ft[10][imaxk] = this->ft[14][imaxk];
                    this->ft[13][imaxk] = this->ft[9][imaxk];
                } else if (this->btymax[ki] == MIRROR) {
                    this->ft[2][imaxk] = this->ft[5][imaxk];
                    this->ft[7][imaxk] = this->ft[9][imaxk];
                    this->ft[8][imaxk] = this->ft[14][imaxk];
                    this->ft[10][imaxk] = this->ft[12][imaxk];
                    this->ft[13][imaxk] = this->ft[11][imaxk];
                }
            }
        }

        //  zmin and zmax
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; i++) {
                int ij = i + this->nx*j;

                //.....zmin.....
                int ijmin = i + this->nx*j + this->nx*this->ny*0;
                int ijminp1 = i + this->nx*j + this->nx*this->ny*1;
                if (this->btzmin[ki] == OUTLET) {
                    this->ft[6][ijmin] = this->ft[6][ijminp1];
                    this->ft[10][ijmin] = this->ft[10][ijminp1];
                    this->ft[11][ijmin] = this->ft[11][ijminp1];
                    this->ft[12][ijmin] = this->ft[12][ijminp1];
                    this->ft[13][ijmin] = this->ft[13][ijminp1];
                } else if (this->btzmin[ki] == BARRIER) {
                    this->ft[6][ijmin] = this->ft[3][ijmin];
                    this->ft[10][ijmin] = this->ft[14][ijmin];
                    this->ft[11][ijmin] = this->ft[7][ijmin];
                    this->ft[12][ijmin] = this->ft[8][ijmin];
                    this->ft[13][ijmin] = this->ft[9][ijmin];
                } else if (this->btzmin[ki] == MIRROR) {
                    this->ft[6][ijmin] = this->ft[3][ijmin];
                    this->ft[10][ijmin] = this->ft[7][ijmin];
                    this->ft[11][ijmin] = this->ft[14][ijmin];
                    this->ft[12][ijmin] = this->ft[9][ijmin];
                    this->ft[13][ijmin] = this->ft[8][ijmin];
                }

                //.....zmax.....
                int ijmax = i + this->nx*j + this->nx*this->ny*(this->nz - 1);
                int ijmaxm1 = i + this->nx*j + this->nx*this->ny*(this->nz - 2);
                if (this->btzmax[ki] == OUTLET) {
                    this->ft[3][ijmax] = this->ft[3][ijmaxm1];
                    this->ft[7][ijmax] = this->ft[7][ijmaxm1];
                    this->ft[8][ijmax] = this->ft[8][ijmaxm1];
                    this->ft[9][ijmax] = this->ft[9][ijmaxm1];
                    this->ft[14][ijmax] = this->ft[14][ijmaxm1];
                } else if (this->btzmax[ki] == BARRIER) {
                    this->ft[3][ijmax] = this->ft[6][ijmax];
                    this->ft[7][ijmax] = this->ft[11][ijmax];
                    this->ft[8][ijmax] = this->ft[12][ijmax];
                    this->ft[9][ijmax] = this->ft[13][ijmax];
                    this->ft[14][ijmax] = this->ft[10][ijmax];
                } else if (this->btzmax[ki] == MIRROR) {
                    this->ft[3][ijmax] = this->ft[6][ijmax];
                    this->ft[7][ijmax] = this->ft[10][ijmax];
                    this->ft[8][ijmax] = this->ft[13][ijmax];
                    this->ft[9][ijmax] = this->ft[12][ijmax];
                    this->ft[14][ijmax] = this->ft[11][ijmax];
                }
            }
        }
    }


    template<class T>
    bool D3Q15<T>::GetBarrier(int _i, int _j, int _k) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny && 0 <= _k && _k < this->nz);
        return this->barrier[0][_i + this->nx*_j + this->nx*this->ny*_k];
    }


    template<class T>
    BOUNDARYTYPE D3Q15<T>::GetBoundary(int _i, int _j, int _k) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny && 0 <= _k && _k < this->nz);
        if (_i == 0) {
            return this->btxmin[_j + this->ny*_k];
        } else if (_i == this->nx - 1) {
            return this->btxmax[_j + this->ny*_k];
        } else if (_j == 0) {
            return this->btymin[_k + this->nz*_i];
        } else if (_j == this->ny - 1) {
            return this->btymax[_k + this->nz*_i];
        } else if (_k == 0) {
            return this->btzmin[_i + this->nx*_j];
        } else if (_k == this->nz - 1) {
            return this->btzmax[_i + this->nx*_j];
        } else {
            if (this->barrier[0][_i + this->nx*_j + this->nx*this->ny*_k]) {
                return BARRIER;
            } else {
                return OTHER;
            }
        }
    }


    template<class T>
    int D3Q15<T>::GetIndex(int _i, int _j, int _k) {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny && 0 <= _k & _k < this->nz);
        return _i + this->nx*_j + this->nx*this->ny*_k;
    }
}