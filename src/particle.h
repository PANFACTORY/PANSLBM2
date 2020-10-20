//*****************************************************************************
//  Title       :   src/particle.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/05
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
        virtual ~D2Q9();

        void SetBarrier(int _i, int _j, bool _isbarrier);
        void SetBoundary(int _i, int _j, BOUNDARYTYPE _boundarytype);

        void Stream();
        
        bool GetBarrier(int _i, int _j) const;
        
        const int nx, ny;
        static const int nc = 9;
        static const T ei[nc];

        T dx, dt;
        T *ft[nc], *ftp1[nc];
        bool *barrier[nc];
        BOUNDARYTYPE *btxmin, *btxmax, *btymin, *btymax;
    };


    template<class T>
    const T D2Q9<T>::ei[D2Q9<T>::nc] = {
        4.0/9.0, 
        1.0/9.0, 
        1.0/9.0, 
        1.0/9.0, 
        1.0/9.0, 
        1.0/36.0, 
        1.0/36.0, 
        1.0/36.0, 
        1.0/36.0
    };


    template<class T>
    D2Q9<T>::D2Q9(int _nx, int _ny) : nx(_nx), ny(_ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;
        this->dt = 1.0;

        for (int i = 0; i < D2Q9<T>::nc; i++) {
            this->ft[i] = new T[this->nx*this->ny];
            this->ftp1[i] = new T[this->nx*this->ny];
            this->barrier[i] = new bool[this->nx*this->ny];

            for (int j = 0; j < this->nx*this->ny; j++) {
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
    D2Q9<T>::D2Q9(const D2Q9<T>& _p) : nx(_p.nx), ny(_p.ny) {
        this->dx = _p.dx;
        this->dt = _p.dt;
        
        for (int i = 0; i < D2Q9<T>::nc; i++) {
            this->ft[i] = new T[this->nx*this->ny];
            this->ftp1[i] = new T[this->nx*this->ny];
            this->barrier[i] = new bool[this->nx*this->ny];

            for (int j = 0; j < this->nx*this->ny; j++) {
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
        this->barrier[0][this->ny*_i + _j] = _isbarrier;
        this->barrier[1][this->ny*(_i + 1) + _j] = _isbarrier;
        this->barrier[2][this->ny*_i + (_j + 1)] = _isbarrier;
        this->barrier[3][this->ny*(_i - 1) + _j] = _isbarrier;
        this->barrier[4][this->ny*_i + (_j - 1)] = _isbarrier;
        this->barrier[5][this->ny*(_i + 1) + (_j + 1)] = _isbarrier;
        this->barrier[6][this->ny*(_i - 1) + (_j + 1)] = _isbarrier;
        this->barrier[7][this->ny*(_i - 1) + (_j - 1)] = _isbarrier;
        this->barrier[8][this->ny*(_i + 1) + (_j - 1)] = _isbarrier;
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
    void D2Q9<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->ft[0][this->ny*i + j] = this->ftp1[0][this->ny*i + j];
                this->ft[1][this->ny*i + j] = (i == 0) ? this->ftp1[1][this->ny*(this->nx - 1) + j] : this->ftp1[1][this->ny*(i - 1) + j];
                this->ft[2][this->ny*i + j] = (j == 0) ? this->ftp1[2][this->ny*(i + 1) - 1] : this->ftp1[2][this->ny*i + (j - 1)];
                this->ft[3][this->ny*i + j] = (i == this->nx - 1) ? this->ftp1[3][j] : this->ftp1[3][this->ny*(i + 1) + j];
                this->ft[4][this->ny*i + j] = (j == this->ny - 1) ? this->ftp1[4][this->ny*i] : this->ftp1[4][this->ny*i + (j + 1)];
                this->ft[5][this->ny*i + j] = (i == 0) ? ((j == 0) ? this->ftp1[5][this->ny*this->nx - 1] : this->ftp1[5][this->ny*(this->nx - 1) + (j - 1)]) : ((j == 0) ? this->ftp1[5][this->ny*i - 1] : this->ftp1[5][this->ny*(i - 1) + (j - 1)]);
                this->ft[6][this->ny*i + j] = (i == this->nx - 1) ? ((j == 0) ? this->ftp1[6][this->ny - 1] : this->ftp1[6][j - 1]) : ((j == 0) ? this->ftp1[6][this->ny*(i + 2) - 1] : this->ftp1[6][this->ny*(i + 1) + (j - 1)]);
                this->ft[7][this->ny*i + j] = (i == this->nx - 1) ? ((j == this->ny - 1) ? this->ftp1[7][0] : this->ftp1[7][j + 1]) : ((j == this->ny - 1) ? this->ftp1[7][this->ny*(i + 1)] : this->ftp1[7][this->ny*(i + 1) + (j + 1)]);
                this->ft[8][this->ny*i + j] = (i == 0) ? ((j == this->ny - 1) ? this->ftp1[8][this->ny*(this->nx - 1)] : this->ftp1[8][this->ny*(this->nx - 1) + (j + 1)]) : ((j == this->ny - 1) ? this->ftp1[8][this->ny*(i - 1)] : this->ftp1[8][this->ny*(i - 1) + (j + 1)]);
            }
        }
        
        //----------Bouns-Back (inner boundary)----------
        for (int i = 0; i < this->nx*this->ny; i++) {
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
    bool D2Q9<T>::GetBarrier(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->barrier[0][this->ny*_i + _j];
    }
}