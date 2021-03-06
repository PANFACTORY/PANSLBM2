//*****************************************************************************
//  Title       :   src/particle/d2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/21
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
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
        void SetRho(int _i, int _j, T _rho, T _u);                          //  Set boundary condition for NavierStokes  
        void SetU(int _i, int _j, T _ux, T _uy);
        void SetTemperature(int _i, int _j, T _ux, T _uy, T _temperature);  //  Set boundary condition for Advection
        void SetFlux(int _i, int _j, T _ux, T _uy, T _q);
        void SetStress(int _i, int _j, T _tx, T _ty);                       //  Set boundary condition for Elastic
        void SetiRho(int _i, int _j);                                       //  Set boundary condition for Adjoint of NavierStokes
        void SetiU(int _i, int _j, T _ux, T _uy);  
        void SetiUPressureDrop(int _i, int _j, T _ux, T _uy, T _eps = 1);
        void SetiTemperature(int _i, int _j, T _ux, T _uy);                 //  Set boundary condition for Adjoint of Advection
        void SetiFlux(int _i, int _j, T _ux, T _uy);
        void SetiRhoFlux(const D2Q9<T>& _g, int _i, int _j, T _rho, T _ux, T _uy, T _temperature); 
        void SetiStress(int _i, int _j, T _rho, T _tx, T _ty);              //  Set boundary condition for Adjoint of Elastic

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
            this->barrier[i][this->GetIndex(_i + D2Q9<T>::cx[i], _j + D2Q9<T>::cy[i])] = _isbarrier;
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
        int ij = this->GetIndex(_i, _j);
        if (_i == 0 && _j == 0) {
            int ijr = this->GetIndex(_i + 1, _j);
            T uxr = this->ft[1][ijr] - this->ft[3][ijr] + this->ft[5][ijr] - this->ft[6][ijr] - this->ft[7][ijr] + this->ft[8][ijr]; 
            T uyr = this->ft[2][ijr] - this->ft[4][ijr] + this->ft[5][ijr] + this->ft[6][ijr] - this->ft[7][ijr] - this->ft[8][ijr]; 
            int iju = this->GetIndex(_i, _j + 1);
            T uxu = this->ft[1][iju] - this->ft[3][iju] + this->ft[5][iju] - this->ft[6][iju] - this->ft[7][iju] + this->ft[8][iju]; 
            T uyu = this->ft[2][iju] - this->ft[4][iju] + this->ft[5][iju] + this->ft[6][iju] - this->ft[7][iju] - this->ft[8][iju]; 
            int ijru = this->GetIndex(_i + 1, _j + 1);
            T uxru = this->ft[1][ijru] - this->ft[3][ijru] + this->ft[5][ijru] - this->ft[6][ijru] - this->ft[7][ijru] + this->ft[8][ijru]; 
            T uyru = this->ft[2][ijru] - this->ft[4][ijru] + this->ft[5][ijru] + this->ft[6][ijru] - this->ft[7][ijru] - this->ft[8][ijru]; 
            T ux0 = (uxr + uxu + uxru)/3.0;
            T uy0 = (uyr + uyu + uyru)/3.0;
            this->ft[1][ij] = this->ft[3][ij] + 2.0*_rho*ux0/3.0;
            this->ft[2][ij] = this->ft[4][ij] + 2.0*_rho*uy0/3.0;
            this->ft[5][ij] = this->ft[7][ij] + _rho*ux0/6.0 + _rho*uy0/6.0;
            this->ft[6][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 - _rho*ux0/12.0 + _rho*uy0/12.0;
            this->ft[8][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 + _rho*ux0/12.0 - _rho*uy0/12.0;
        } else if (_i == 0 && _j == this->ny - 1) {
            int ijr = this->GetIndex(_i + 1, _j);
            T uxr = this->ft[1][ijr] - this->ft[3][ijr] + this->ft[5][ijr] - this->ft[6][ijr] - this->ft[7][ijr] + this->ft[8][ijr]; 
            T uyr = this->ft[2][ijr] - this->ft[4][ijr] + this->ft[5][ijr] + this->ft[6][ijr] - this->ft[7][ijr] - this->ft[8][ijr]; 
            int ijd = this->GetIndex(_i, _j - 1);
            T uxd = this->ft[1][ijd] - this->ft[3][ijd] + this->ft[5][ijd] - this->ft[6][ijd] - this->ft[7][ijd] + this->ft[8][ijd]; 
            T uyd = this->ft[2][ijd] - this->ft[4][ijd] + this->ft[5][ijd] + this->ft[6][ijd] - this->ft[7][ijd] - this->ft[8][ijd]; 
            int ijrd = this->GetIndex(_i + 1, _j - 1);
            T uxrd = this->ft[1][ijrd] - this->ft[3][ijrd] + this->ft[5][ijrd] - this->ft[6][ijrd] - this->ft[7][ijrd] + this->ft[8][ijrd]; 
            T uyrd = this->ft[2][ijrd] - this->ft[4][ijrd] + this->ft[5][ijrd] + this->ft[6][ijrd] - this->ft[7][ijrd] - this->ft[8][ijrd]; 
            T ux0 = (uxr + uxd + uxrd)/3.0;
            T uy0 = (uyr + uyd + uyrd)/3.0;
            this->ft[1][ij] = this->ft[3][ij] + 2.0*_rho*ux0/3.0;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*_rho*uy0/3.0;
            this->ft[8][ij] = this->ft[6][ij] + _rho*ux0/6.0 - _rho*uy0/6.0;
            this->ft[5][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 + _rho*ux0/12.0 + _rho*uy0/12.0;
            this->ft[7][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 - _rho*ux0/12.0 - _rho*uy0/12.0;
        } else if (_i == this->nx - 1 && _j == 0) {
            int ijl = this->GetIndex(_i - 1, _j);
            T uxl = this->ft[1][ijl] - this->ft[3][ijl] + this->ft[5][ijl] - this->ft[6][ijl] - this->ft[7][ijl] + this->ft[8][ijl]; 
            T uyl = this->ft[2][ijl] - this->ft[4][ijl] + this->ft[5][ijl] + this->ft[6][ijl] - this->ft[7][ijl] - this->ft[8][ijl]; 
            int iju = this->GetIndex(_i, _j + 1);
            T uxu = this->ft[1][iju] - this->ft[3][iju] + this->ft[5][iju] - this->ft[6][iju] - this->ft[7][iju] + this->ft[8][iju]; 
            T uyu = this->ft[2][iju] - this->ft[4][iju] + this->ft[5][iju] + this->ft[6][iju] - this->ft[7][iju] - this->ft[8][iju]; 
            int ijlu = this->GetIndex(_i - 1, _j + 1);
            T uxlu = this->ft[1][ijlu] - this->ft[3][ijlu] + this->ft[5][ijlu] - this->ft[6][ijlu] - this->ft[7][ijlu] + this->ft[8][ijlu]; 
            T uylu = this->ft[2][ijlu] - this->ft[4][ijlu] + this->ft[5][ijlu] + this->ft[6][ijlu] - this->ft[7][ijlu] - this->ft[8][ijlu]; 
            T ux0 = (uxl + uxu + uxlu)/3.0;
            T uy0 = (uyl + uyu + uylu)/3.0;
            this->ft[3][ij] = this->ft[1][ij] - 2.0*_rho*ux0/3.0;
            this->ft[2][ij] = this->ft[4][ij] + 2.0*_rho*uy0/3.0;
            this->ft[6][ij] = this->ft[7][ij] - _rho*ux0/6.0 + _rho*uy0/6.0;
            this->ft[5][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 + _rho*ux0/12.0 + _rho*uy0/12.0;
            this->ft[7][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 - _rho*ux0/12.0 - _rho*uy0/12.0;
        } else if (_i == this->nx - 1 && _j == this->ny - 1) {
            int ijl = this->GetIndex(_i - 1, _j);
            T uxl = this->ft[1][ijl] - this->ft[3][ijl] + this->ft[5][ijl] - this->ft[6][ijl] - this->ft[7][ijl] + this->ft[8][ijl]; 
            T uyl = this->ft[2][ijl] - this->ft[4][ijl] + this->ft[5][ijl] + this->ft[6][ijl] - this->ft[7][ijl] - this->ft[8][ijl]; 
            int ijd = this->GetIndex(_i, _j - 1);
            T uxd = this->ft[1][ijd] - this->ft[3][ijd] + this->ft[5][ijd] - this->ft[6][ijd] - this->ft[7][ijd] + this->ft[8][ijd]; 
            T uyd = this->ft[2][ijd] - this->ft[4][ijd] + this->ft[5][ijd] + this->ft[6][ijd] - this->ft[7][ijd] - this->ft[8][ijd]; 
            int ijld = this->GetIndex(_i - 1, _j - 1);
            T uxld = this->ft[1][ijld] - this->ft[3][ijld] + this->ft[5][ijld] - this->ft[6][ijld] - this->ft[7][ijld] + this->ft[8][ijld]; 
            T uyld = this->ft[2][ijld] - this->ft[4][ijld] + this->ft[5][ijld] + this->ft[6][ijld] - this->ft[7][ijld] - this->ft[8][ijld]; 
            T ux0 = (uxl + uxd + uxld)/3.0;
            T uy0 = (uyl + uyd + uyld)/3.0;
            this->ft[3][ij] = this->ft[1][ij] - 2.0*_rho*ux0/3.0;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*_rho*uy0/3.0;
            this->ft[7][ij] = this->ft[5][ij] - _rho*ux0/6.0 - _rho*uy0/6.0;
            this->ft[6][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 - _rho*ux0/12.0 + _rho*uy0/12.0;
            this->ft[8][ij] = (_rho - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 + _rho*ux0/12.0 - _rho*uy0/12.0;
        } else if (_i == 0) {
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
            T uy0 = -1.0 + (this->ft[0][ij] + this->ft[1][ij] + this->ft[3][ij] + 2.0*(this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]))/_rho;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*_rho*uy0/3.0;
            this->ft[7][ij] = this->ft[5][ij] + 0.5*(this->ft[1][ij] - this->ft[3][ij]) - _rho*_u/2.0 - _rho*uy0/6.0;
            this->ft[8][ij] = this->ft[6][ij] - 0.5*(this->ft[1][ij] - this->ft[3][ij]) + _rho*_u/2.0 - _rho*uy0/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetU(int _i, int _j, T _ux, T _uy) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0 && _j == 0) {
            int ijr = this->GetIndex(_i + 1, _j);
            T rhor = this->ft[0][ijr] + this->ft[1][ijr] + this->ft[2][ijr] + this->ft[3][ijr] + this->ft[4][ijr] + this->ft[5][ijr] + this->ft[6][ijr] + this->ft[7][ijr] + this->ft[8][ijr]; 
            int iju = this->GetIndex(_i, _j + 1);
            T rhou = this->ft[0][iju] + this->ft[1][iju] + this->ft[2][iju] + this->ft[3][iju] + this->ft[4][iju] + this->ft[5][iju] + this->ft[6][iju] + this->ft[7][iju] + this->ft[8][iju]; 
            int ijru = this->GetIndex(_i + 1, _j + 1);
            T rhoru = this->ft[0][ijru] + this->ft[1][ijru] + this->ft[2][ijru] + this->ft[3][ijru] + this->ft[4][ijru] + this->ft[5][ijru] + this->ft[6][ijru] + this->ft[7][ijru] + this->ft[8][ijru]; 
            T rho0 = (rhor + rhou + rhoru)/3.0;
            this->ft[1][ij] = this->ft[3][ij] + 2.0*rho0*_ux/3.0;
            this->ft[2][ij] = this->ft[4][ij] + 2.0*rho0*_uy/3.0;
            this->ft[5][ij] = this->ft[7][ij] + rho0*_ux/6.0 + rho0*_uy/6.0;
            this->ft[6][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 - rho0*_ux/12.0 + rho0*_uy/12.0;
            this->ft[8][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 + rho0*_ux/12.0 - rho0*_uy/12.0;
        } else if (_i == 0 && _j == this->ny - 1) {
            int ijr = this->GetIndex(_i + 1, _j);
            T rhor = this->ft[0][ijr] + this->ft[1][ijr] + this->ft[2][ijr] + this->ft[3][ijr] + this->ft[4][ijr] + this->ft[5][ijr] + this->ft[6][ijr] + this->ft[7][ijr] + this->ft[8][ijr]; 
            int ijd = this->GetIndex(_i, _j - 1);
            T rhod = this->ft[0][ijd] + this->ft[1][ijd] + this->ft[2][ijd] + this->ft[3][ijd] + this->ft[4][ijd] + this->ft[5][ijd] + this->ft[6][ijd] + this->ft[7][ijd] + this->ft[8][ijd]; 
            int ijrd = this->GetIndex(_i + 1, _j - 1);
            T rhord = this->ft[0][ijrd] + this->ft[1][ijrd] + this->ft[2][ijrd] + this->ft[3][ijrd] + this->ft[4][ijrd] + this->ft[5][ijrd] + this->ft[6][ijrd] + this->ft[7][ijrd] + this->ft[8][ijrd]; 
            T rho0 = (rhor + rhod + rhord)/3.0;
            this->ft[1][ij] = this->ft[3][ij] + 2.0*rho0*_ux/3.0;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*rho0*_uy/3.0;
            this->ft[8][ij] = this->ft[6][ij] + rho0*_ux/6.0 - rho0*_uy/6.0;
            this->ft[5][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 + rho0*_ux/12.0 + rho0*_uy/12.0;
            this->ft[7][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 - rho0*_ux/12.0 - rho0*_uy/12.0;
        } else if (_i == this->nx - 1 && _j == 0) {
            int ijl = this->GetIndex(_i - 1, _j);
            T rhol = this->ft[0][ijl] + this->ft[1][ijl] + this->ft[2][ijl] + this->ft[3][ijl] + this->ft[4][ijl] + this->ft[5][ijl] + this->ft[6][ijl] + this->ft[7][ijl] + this->ft[8][ijl]; 
            int iju = this->GetIndex(_i, _j + 1);
            T rhou = this->ft[0][iju] + this->ft[1][iju] + this->ft[2][iju] + this->ft[3][iju] + this->ft[4][iju] + this->ft[5][iju] + this->ft[6][iju] + this->ft[7][iju] + this->ft[8][iju]; 
            int ijlu = this->GetIndex(_i - 1, _j + 1);
            T rholu = this->ft[0][ijlu] + this->ft[1][ijlu] + this->ft[2][ijlu] + this->ft[3][ijlu] + this->ft[4][ijlu] + this->ft[5][ijlu] + this->ft[6][ijlu] + this->ft[7][ijlu] + this->ft[8][ijlu]; 
            T rho0 = (rhol + rhou + rholu)/3.0;
            this->ft[3][ij] = this->ft[1][ij] - 2.0*rho0*_ux/3.0;
            this->ft[2][ij] = this->ft[4][ij] + 2.0*rho0*_uy/3.0;
            this->ft[6][ij] = this->ft[7][ij] - rho0*_ux/6.0 + rho0*_uy/6.0;
            this->ft[5][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 + rho0*_ux/12.0 + rho0*_uy/12.0;
            this->ft[7][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[8][ij])/2.0 - rho0*_ux/12.0 - rho0*_uy/12.0;
        } else if (_i == this->nx - 1 && _j == this->ny - 1) {
            int ijl = this->GetIndex(_i - 1, _j);
            T rhol = this->ft[0][ijl] + this->ft[1][ijl] + this->ft[2][ijl] + this->ft[3][ijl] + this->ft[4][ijl] + this->ft[5][ijl] + this->ft[6][ijl] + this->ft[7][ijl] + this->ft[8][ijl]; 
            int ijd = this->GetIndex(_i, _j - 1);
            T rhod = this->ft[0][ijd] + this->ft[1][ijd] + this->ft[2][ijd] + this->ft[3][ijd] + this->ft[4][ijd] + this->ft[5][ijd] + this->ft[6][ijd] + this->ft[7][ijd] + this->ft[8][ijd]; 
            int ijld = this->GetIndex(_i - 1, _j - 1);
            T rhold = this->ft[0][ijld] + this->ft[1][ijld] + this->ft[2][ijld] + this->ft[3][ijld] + this->ft[4][ijld] + this->ft[5][ijld] + this->ft[6][ijld] + this->ft[7][ijld] + this->ft[8][ijld]; 
            T rho0 = (rhol + rhod + rhold)/3.0;
            this->ft[3][ij] = this->ft[1][ij] - 2.0*rho0*_ux/3.0;
            this->ft[4][ij] = this->ft[2][ij] - 2.0*rho0*_uy/3.0;
            this->ft[7][ij] = this->ft[5][ij] - rho0*_ux/6.0 - rho0*_uy/6.0;
            this->ft[6][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 - rho0*_ux/12.0 + rho0*_uy/12.0;
            this->ft[8][ij] = (rho0 - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[7][ij])/2.0 + rho0*_ux/12.0 - rho0*_uy/12.0;
        } else if (_i == 0) {
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
    void D2Q9<T>::SetTemperature(int _i, int _j, T _ux, T _uy, T _temperature) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0 && _j == 0) {
            T temperature0 = 36.0*(_temperature - (this->ft[0][ij] + this->ft[3][ij] + this->ft[4][ij] + this->ft[7][ij]))/(11.0 + 15.0*_ux + 15.0*_uy);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == 0 && _j == this->ny - 1) {
            T temperature0 = 36.0*(_temperature - (this->ft[0][ij] + this->ft[2][ij] + this->ft[3][ij] + this->ft[6][ij]))/(11.0 + 15.0*_ux - 15.0*_uy);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1 && _j == 0) {
            T temperature0 = 36.0*(_temperature - (this->ft[0][ij] + this->ft[1][ij] + this->ft[4][ij] + this->ft[8][ij]))/(11.0 - 15.0*_ux + 15.0*_uy);
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1 && _j == this->ny - 1) {
            T temperature0 = 36.0*(_temperature - (this->ft[0][ij] + this->ft[1][ij] + this->ft[2][ij] + this->ft[5][ij]))/(11.0 - 15.0*_ux - 15.0*_uy);
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[6][ij] - this->ft[7][ij])/(1.0 + 3.0*_ux);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[4][ij] - this->ft[5][ij] - this->ft[8][ij])/(1.0 - 3.0*_ux);
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[3][ij] - this->ft[4][ij] - this->ft[7][ij] - this->ft[8][ij])/(1.0 + 3.0*_uy);
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(_temperature - this->ft[0][ij] - this->ft[1][ij] - this->ft[2][ij] - this->ft[3][ij] - this->ft[5][ij] - this->ft[6][ij])/(1.0 - 3.0*_uy);
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }
    
    
    template<class T>
    void D2Q9<T>::SetFlux(int _i, int _j, T _ux, T _uy, T _q) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0 && _j == 0) {
            T temperature0 = 18.0*(sqrt(2.0)*_q + this->ft[3][ij] + this->ft[4][ij] + 2.0*this->ft[7][ij])/(5.0 - 9.0*_ux - 9.0*_uy);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == 0 && _j == this->ny - 1) {
            T temperature0 = 18.0*(sqrt(2.0)*_q + this->ft[2][ij] + this->ft[3][ij] + 2.0*this->ft[6][ij])/(5.0 - 9.0*_ux + 9.0*_uy);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1 && _j == 0) {
            T temperature0 = 18.0*(sqrt(2.0)*_q + this->ft[1][ij] + this->ft[4][ij] + 2.0*this->ft[8][ij])/(5.0 + 9.0*_ux - 9.0*_uy);
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1 && _j == this->ny - 1) {
            T temperature0 = 18.0*(sqrt(2.0)*_q + this->ft[1][ij] + this->ft[2][ij] + 2.0*this->ft[5][ij])/(5.0 + 9.0*_ux + 9.0*_uy);
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == 0) {
            T temperature0 = 6.0*(_q + this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/(1.0 - 3.0*_ux);
            this->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(_q + this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/(1.0 + 3.0*_ux);
            this->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_q + this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/(1.0 - 3.0*_uy);
            this->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(_q + this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/(1.0 + 3.0*_uy);
            this->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetStress(int _i, int _j, T _tx, T _ty) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            this->ft[1][ij] = this->ft[3][ij] - 4.0*(this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + 2.0*_tx/3.0;
            this->ft[5][ij] = this->ft[6][ij] - (this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + (_tx + 3.0*_ty)/6.0;
            this->ft[8][ij] = this->ft[7][ij] - (this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + (_tx - 3.0*_ty)/6.0;
        } else if (_i == this->nx - 1) {
            this->ft[3][ij] = this->ft[1][ij] - 4.0*(this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - 2.0*_tx/3.0;
            this->ft[6][ij] = this->ft[5][ij] - (this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - (_tx - 3.0*_ty)/6.0;
            this->ft[7][ij] = this->ft[8][ij] - (this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - (_tx + 3.0*_ty)/6.0;
        } else if (_j == 0) {
            this->ft[2][ij] = this->ft[4][ij] - 4.0*(this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + 2.0*_ty/3.0;
            this->ft[5][ij] = this->ft[8][ij] - (this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + (_ty + 3.0*_tx)/6.0;
            this->ft[6][ij] = this->ft[7][ij] - (this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + (_ty - 3.0*_tx)/6.0;
        } else if (_j == this->ny - 1) {
            this->ft[4][ij] = this->ft[2][ij] - 4.0*(this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - 2.0*_ty/3.0;
            this->ft[7][ij] = this->ft[6][ij] - (this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - (_ty + 3.0*_tx)/6.0;
            this->ft[8][ij] = this->ft[5][ij] - (this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - (_ty - 3.0*_tx)/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiRho(int _i, int _j) {
        int ij = this->GetIndex(_i, _j);
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
        int ij = this->GetIndex(_i, _j);
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
        int ij = this->GetIndex(_i, _j);
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
    void D2Q9<T>::SetiTemperature(int _i, int _j, T _ux, T _uy) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            T rho0 = -(4.0*(1.0 + 3.0*_ux)*this->ft[1][ij] + (1.0 + 3.0*_ux + 3.0*_uy)*this->ft[5][ij] + (1.0 + 3.0*_ux - 3.0*_uy)*this->ft[8][ij])/(6.0*(1.0 + 3.0*_ux));
            this->ft[3][ij] = rho0;
            this->ft[6][ij] = rho0;
            this->ft[7][ij] = rho0;
        } else if (_i == this->nx - 1) {
            T rho0 = -(4.0*(1.0 - 3.0*_ux)*this->ft[3][ij] + (1.0 - 3.0*_ux + 3.0*_uy)*this->ft[6][ij] + (1.0 - 3.0*_ux - 3.0*_uy)*this->ft[7][ij])/(6.0*(1.0 - 3.0*_ux));
            this->ft[1][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == 0) {
            T rho0 = -(4.0*(1.0 + 3.0*_uy)*this->ft[2][ij] + (1.0 + 3.0*_ux + 3.0*_uy)*this->ft[5][ij] + (1.0 - 3.0*_ux + 3.0*_uy)*this->ft[6][ij])/(6.0*(1.0 + 3.0*_uy));
            this->ft[4][ij] = rho0;
            this->ft[7][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = -(4.0*(1.0 - 3.0*_uy)*this->ft[4][ij] + (1.0 - 3.0*_ux - 3.0*_uy)*this->ft[7][ij] + (1.0 + 3.0*_ux - 3.0*_uy)*this->ft[8][ij])/(6.0*(1.0 - 3.0*_uy));
            this->ft[2][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[6][ij] = rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiFlux(int _i, int _j, T _ux, T _uy) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            T rho0 = (1.0 + 3.0*_ux)/(6.0*(1.0 - 3.0*_ux))*(4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij]) + _uy/(2.0*(1.0 - 3.0*_ux))*(this->ft[5][ij] - this->ft[8][ij]);
            this->ft[3][ij] = rho0;
            this->ft[6][ij] = rho0;
            this->ft[7][ij] = rho0;
        } else if (_i == this->nx - 1) {          
            T rho0 = -(1.0 - 3.0*_ux)/(6.0*(1.0 + 3.0*_ux))*(4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij]) - _uy/(2.0*(1.0 + 3.0*_ux))*(this->ft[6][ij] - this->ft[7][ij]);
            this->ft[1][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == 0) {
            T rho0 = (1.0 + 3.0*_uy)/(6.0*(1.0 - 3.0*_uy))*(4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij]) + _ux/(2.0*(1.0 - 3.0*_uy))*(this->ft[5][ij] - this->ft[6][ij]);
            this->ft[4][ij] = rho0;
            this->ft[7][ij] = rho0;
            this->ft[8][ij] = rho0;
        } else if (_j == this->ny - 1) {
            T rho0 = -(1.0 - 3.0*_uy)/(6.0*(1.0 + 3.0*_uy))*(4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij]) - _ux/(2.0*(1.0 + 3.0*_uy))*(this->ft[8][ij] - this->ft[7][ij]);
            this->ft[2][ij] = rho0;
            this->ft[5][ij] = rho0;
            this->ft[6][ij] = rho0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiRhoFlux(const D2Q9<T>& _g, int _i, int _j, T _rho, T _ux, T _uy, T _temperature) {
        int ij = this->GetIndex(_i, _j);
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
    void D2Q9<T>::SetiStress(int _i, int _j, T _rho, T _tx, T _ty) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            this->ft[3][ij] = this->ft[1][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*_tx/_rho;
            this->ft[6][ij] = this->ft[5][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*(_tx - _ty)/_rho;
            this->ft[7][ij] = this->ft[8][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*(_tx + _ty)/_rho;
        } else if (_i == this->nx - 1) {
            this->ft[1][ij] = this->ft[3][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*_tx/_rho;
            this->ft[5][ij] = this->ft[6][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*(_tx + _ty)/_rho;
            this->ft[8][ij] = this->ft[7][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*(_tx - _ty)/_rho;
        } else if (_j == 0) {
            this->ft[4][ij] = this->ft[2][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*_ty/_rho;
            this->ft[7][ij] = this->ft[6][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*(_ty + _tx)/_rho;
            this->ft[8][ij] = this->ft[5][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*(_ty - _tx)/_rho;
        } else if (_j == this->ny - 1) {
            this->ft[2][ij] = this->ft[4][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*_ty/_rho;
            this->ft[5][ij] = this->ft[8][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*(_ty + _tx)/_rho;
            this->ft[6][ij] = this->ft[7][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*(_ty - _tx)/_rho;
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
                    this->ft[k][this->GetIndex(i, j)] = this->ftp1[k][this->GetIndex(ip1, jp1)];
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
                this->ft[1][this->GetIndex(0, j)] = this->ft[1][this->GetIndex(1, j)];
                this->ft[5][this->GetIndex(0, j)] = this->ft[5][this->GetIndex(1, j)];
                this->ft[8][this->GetIndex(0, j)] = this->ft[8][this->GetIndex(1, j)];
            } else if (this->btxmin[j] == BARRIER) {
                this->ft[1][this->GetIndex(0, j)] = this->ftp1[3][this->GetIndex(0, j)];
                this->ft[5][this->GetIndex(0, j)] = this->ftp1[7][this->GetIndex(0, j)];
                this->ft[8][this->GetIndex(0, j)] = this->ftp1[6][this->GetIndex(0, j)];
            } else if (this->btxmin[j] == MIRROR) {
                this->ft[1][this->GetIndex(0, j)] = this->ftp1[3][this->GetIndex(0, j)];
                this->ft[5][this->GetIndex(0, j)] = this->ftp1[6][this->GetIndex(0, j)];
                this->ft[8][this->GetIndex(0, j)] = this->ftp1[7][this->GetIndex(0, j)];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->ft[3][this->GetIndex(this->nx - 1, j)] = this->ft[3][this->GetIndex(this->nx - 2, j)];
                this->ft[6][this->GetIndex(this->nx - 1, j)] = this->ft[6][this->GetIndex(this->nx - 2, j)];
                this->ft[7][this->GetIndex(this->nx - 1, j)] = this->ft[7][this->GetIndex(this->nx - 2, j)];
            } else if (this->btxmax[j] == BARRIER) {
                this->ft[3][this->GetIndex(this->nx - 1, j)] = this->ftp1[1][this->GetIndex(this->nx - 1, j)];
                this->ft[6][this->GetIndex(this->nx - 1, j)] = this->ftp1[8][this->GetIndex(this->nx - 1, j)];
                this->ft[7][this->GetIndex(this->nx - 1, j)] = this->ftp1[5][this->GetIndex(this->nx - 1, j)];
            } else if (this->btxmax[j] == MIRROR) {
                this->ft[3][this->GetIndex(this->nx - 1, j)] = this->ftp1[1][this->GetIndex(this->nx - 1, j)];
                this->ft[6][this->GetIndex(this->nx - 1, j)] = this->ftp1[5][this->GetIndex(this->nx - 1, j)];
                this->ft[7][this->GetIndex(this->nx - 1, j)] = this->ftp1[8][this->GetIndex(this->nx - 1, j)];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->ft[2][this->GetIndex(i, 0)] = this->ft[2][this->GetIndex(i, 1)];
                this->ft[5][this->GetIndex(i, 0)] = this->ft[5][this->GetIndex(i, 1)];
                this->ft[6][this->GetIndex(i, 0)] = this->ft[6][this->GetIndex(i, 1)];
            } else if (this->btymin[i] == BARRIER) {
                this->ft[2][this->GetIndex(i, 0)] = this->ftp1[4][this->GetIndex(i, 0)];
                this->ft[5][this->GetIndex(i, 0)] = this->ftp1[7][this->GetIndex(i, 0)];
                this->ft[6][this->GetIndex(i, 0)] = this->ftp1[8][this->GetIndex(i, 0)];
            } else if (this->btymin[i] == MIRROR) {
                this->ft[2][this->GetIndex(i, 0)] = this->ftp1[4][this->GetIndex(i, 0)];
                this->ft[5][this->GetIndex(i, 0)] = this->ftp1[8][this->GetIndex(i, 0)];
                this->ft[6][this->GetIndex(i, 0)] = this->ftp1[7][this->GetIndex(i, 0)];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->ft[4][this->GetIndex(i, this->ny - 1)] = this->ft[4][this->GetIndex(i, this->ny - 2)];
                this->ft[7][this->GetIndex(i, this->ny - 1)] = this->ft[7][this->GetIndex(i, this->ny - 2)];
                this->ft[8][this->GetIndex(i, this->ny - 1)] = this->ft[8][this->GetIndex(i, this->ny - 2)];
            } else if (this->btymax[i] == BARRIER) {
                this->ft[4][this->GetIndex(i, this->ny - 1)] = this->ftp1[2][this->GetIndex(i, this->ny - 1)];
                this->ft[7][this->GetIndex(i, this->ny - 1)] = this->ftp1[5][this->GetIndex(i, this->ny - 1)];
                this->ft[8][this->GetIndex(i, this->ny - 1)] = this->ftp1[6][this->GetIndex(i, this->ny - 1)];
            } else if (this->btymax[i] == MIRROR) {
                this->ft[4][this->GetIndex(i, this->ny - 1)] = this->ftp1[2][this->GetIndex(i, this->ny - 1)];
                this->ft[7][this->GetIndex(i, this->ny - 1)] = this->ftp1[6][this->GetIndex(i, this->ny - 1)];
                this->ft[8][this->GetIndex(i, this->ny - 1)] = this->ftp1[5][this->GetIndex(i, this->ny - 1)];
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
                    this->ft[k][this->GetIndex(i, j)] = this->ftp1[k][this->GetIndex(ip1, jp1)];
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
            if (this->btxmin[j] == OUTLET) {
                this->ft[3][this->GetIndex(0, j)] = this->ft[3][this->GetIndex(1, j)];
                this->ft[6][this->GetIndex(0, j)] = this->ft[6][this->GetIndex(1, j)];
                this->ft[7][this->GetIndex(0, j)] = this->ft[7][this->GetIndex(1, j)];
            } else if (this->btxmin[j] == BARRIER) {
                this->ft[3][this->GetIndex(0, j)] = this->ftp1[1][this->GetIndex(0, j)];
                this->ft[6][this->GetIndex(0, j)] = this->ftp1[8][this->GetIndex(0, j)];
                this->ft[7][this->GetIndex(0, j)] = this->ftp1[5][this->GetIndex(0, j)];
            } else if (this->btxmin[j] == MIRROR) {
                this->ft[3][this->GetIndex(0, j)] = this->ftp1[1][this->GetIndex(0, j)];
                this->ft[6][this->GetIndex(0, j)] = this->ftp1[5][this->GetIndex(0, j)];
                this->ft[7][this->GetIndex(0, j)] = this->ftp1[8][this->GetIndex(0, j)];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->ft[1][this->GetIndex(this->nx - 1, j)] = this->ft[1][this->GetIndex(this->nx - 2, j)];
                this->ft[5][this->GetIndex(this->nx - 1, j)] = this->ft[5][this->GetIndex(this->nx - 2, j)];
                this->ft[8][this->GetIndex(this->nx - 1, j)] = this->ft[8][this->GetIndex(this->nx - 2, j)];
            } else if (this->btxmax[j] == BARRIER) {
                this->ft[1][this->GetIndex(this->nx - 1, j)] = this->ftp1[3][this->GetIndex(this->nx - 1, j)];
                this->ft[5][this->GetIndex(this->nx - 1, j)] = this->ftp1[7][this->GetIndex(this->nx - 1, j)];
                this->ft[8][this->GetIndex(this->nx - 1, j)] = this->ftp1[6][this->GetIndex(this->nx - 1, j)];
            } else if (this->btxmax[j] == MIRROR) {
                this->ft[1][this->GetIndex(this->nx - 1, j)] = this->ftp1[3][this->GetIndex(this->nx - 1, j)];
                this->ft[5][this->GetIndex(this->nx - 1, j)] = this->ftp1[6][this->GetIndex(this->nx - 1, j)];
                this->ft[8][this->GetIndex(this->nx - 1, j)] = this->ftp1[7][this->GetIndex(this->nx - 1, j)];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->ft[4][this->GetIndex(i, 0)] = this->ft[4][this->GetIndex(i, 1)];
                this->ft[7][this->GetIndex(i, 0)] = this->ft[7][this->GetIndex(i, 1)];
                this->ft[8][this->GetIndex(i, 0)] = this->ft[8][this->GetIndex(i, 1)];
            } else if (this->btymin[i] == BARRIER) {
                this->ft[4][this->GetIndex(i, 0)] = this->ftp1[2][this->GetIndex(i, 0)];
                this->ft[7][this->GetIndex(i, 0)] = this->ftp1[5][this->GetIndex(i, 0)];
                this->ft[8][this->GetIndex(i, 0)] = this->ftp1[6][this->GetIndex(i, 0)];
            } else if (this->btymin[i] == MIRROR) {
                this->ft[4][this->GetIndex(i, 0)] = this->ftp1[2][this->GetIndex(i, 0)];
                this->ft[7][this->GetIndex(i, 0)] = this->ftp1[6][this->GetIndex(i, 0)];
                this->ft[8][this->GetIndex(i, 0)] = this->ftp1[5][this->GetIndex(i, 0)];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->ft[2][this->GetIndex(i, this->ny - 1)] = this->ft[2][this->GetIndex(i, this->ny - 2)];
                this->ft[5][this->GetIndex(i, this->ny - 1)] = this->ft[5][this->GetIndex(i, this->ny - 2)];
                this->ft[6][this->GetIndex(i, this->ny - 1)] = this->ft[6][this->GetIndex(i, this->ny - 2)];
            } else if (this->btymax[i] == BARRIER) {
                this->ft[2][this->GetIndex(i, this->ny - 1)] = this->ftp1[4][this->GetIndex(i, this->ny - 1)];
                this->ft[5][this->GetIndex(i, this->ny - 1)] = this->ftp1[7][this->GetIndex(i, this->ny - 1)];
                this->ft[6][this->GetIndex(i, this->ny - 1)] = this->ftp1[8][this->GetIndex(i, this->ny - 1)];
            } else if (this->btymax[i] == MIRROR) {
                this->ft[2][this->GetIndex(i, this->ny - 1)] = this->ftp1[4][this->GetIndex(i, this->ny - 1)];
                this->ft[5][this->GetIndex(i, this->ny - 1)] = this->ftp1[8][this->GetIndex(i, this->ny - 1)];
                this->ft[6][this->GetIndex(i, this->ny - 1)] = this->ftp1[7][this->GetIndex(i, this->ny - 1)];
            }
        }
    }


    template<class T>
    bool D2Q9<T>::GetBarrier(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->barrier[0][this->GetIndex(_i, _j)];
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
            if (this->barrier[0][this->GetIndex(_i, _j)]) {
                return BARRIER;
            } else {
                return OTHER;
            }
        }
    }


    template<class T>
    int D2Q9<T>::GetIndex(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return _i + this->nx*_j;
    }
}