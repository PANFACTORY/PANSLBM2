//*****************************************************************************
//  Title       :   src/d2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/05
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
        virtual ~D2Q9();

        void SetBarrier(int _i, int _j, bool _isbarrier);
        void SetBoundary(int _i, int _j, BOUNDARYTYPE _boundarytype);

        void Stream();
        
        bool GetBarrier(int _i, int _j) const;
        
        const int nx, ny;
        T dx, dt, t0, t1, t2;
        T *f0t, *f1t, *f2t, *f3t, *f4t, *f5t, *f6t, *f7t, *f8t, *f0tp1, *f1tp1, *f2tp1, *f3tp1, *f4tp1, *f5tp1, *f6tp1, *f7tp1, *f8tp1;
        bool *barrier0, *barrier1, *barrier2, *barrier3, *barrier4, *barrier5, *barrier6, *barrier7, *barrier8;
        BOUNDARYTYPE *btxmin, *btxmax, *btymin, *btymax;
    };


    template<class T>
    D2Q9<T>::D2Q9(int _nx, int _ny) : nx(_nx), ny(_ny) {
        assert(0 < _nx && 0 < _ny);
        this->dx = 1.0;     this->dt = 1.0;
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
    D2Q9<T>::D2Q9(const D2Q9<T>& _p) : nx(_p.nx), ny(_p.ny) {
        this->dx = _p.dx;   this->dt = _p.dt;
        this->t0 = _p.t0;   this->t1 = _p.t1;   this->t2 = _p.t2;

        this->f0t = new T[this->nx*this->ny];   this->f0tp1 = new T[this->nx*this->ny];     this->barrier0 = new bool[this->nx*this->ny];
        this->f1t = new T[this->nx*this->ny];   this->f1tp1 = new T[this->nx*this->ny];     this->barrier1 = new bool[this->nx*this->ny];        
        this->f2t = new T[this->nx*this->ny];   this->f2tp1 = new T[this->nx*this->ny];     this->barrier2 = new bool[this->nx*this->ny];     
        this->f3t = new T[this->nx*this->ny];   this->f3tp1 = new T[this->nx*this->ny];     this->barrier3 = new bool[this->nx*this->ny];    
        this->f4t = new T[this->nx*this->ny];   this->f4tp1 = new T[this->nx*this->ny];     this->barrier4 = new bool[this->nx*this->ny];
        this->f5t = new T[this->nx*this->ny];   this->f5tp1 = new T[this->nx*this->ny];     this->barrier5 = new bool[this->nx*this->ny]; 
        this->f6t = new T[this->nx*this->ny];   this->f6tp1 = new T[this->nx*this->ny];     this->barrier6 = new bool[this->nx*this->ny];
        this->f7t = new T[this->nx*this->ny];   this->f7tp1 = new T[this->nx*this->ny];     this->barrier7 = new bool[this->nx*this->ny];
        this->f8t = new T[this->nx*this->ny];   this->f8tp1 = new T[this->nx*this->ny];     this->barrier8 = new bool[this->nx*this->ny];

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->f0t[i] = _p.f0t[i];   this->f0tp1[i] = _p.f0tp1[i];   this->barrier0[i] = _p.barrier0[i];
            this->f1t[i] = _p.f1t[i];   this->f1tp1[i] = _p.f1tp1[i];   this->barrier1[i] = _p.barrier1[i];
            this->f2t[i] = _p.f2t[i];   this->f2tp1[i] = _p.f2tp1[i];   this->barrier2[i] = _p.barrier2[i];
            this->f3t[i] = _p.f3t[i];   this->f3tp1[i] = _p.f3tp1[i];   this->barrier3[i] = _p.barrier3[i];
            this->f4t[i] = _p.f4t[i];   this->f4tp1[i] = _p.f4tp1[i];   this->barrier4[i] = _p.barrier4[i];
            this->f5t[i] = _p.f5t[i];   this->f5tp1[i] = _p.f5tp1[i];   this->barrier5[i] = _p.barrier5[i];
            this->f6t[i] = _p.f6t[i];   this->f6tp1[i] = _p.f6tp1[i];   this->barrier6[i] = _p.barrier6[i];
            this->f7t[i] = _p.f7t[i];   this->f7tp1[i] = _p.f7tp1[i];   this->barrier7[i] = _p.barrier7[i];
            this->f8t[i] = _p.f8t[i];   this->f8tp1[i] = _p.f8tp1[i];   this->barrier8[i] = _p.barrier8[i];
        }

        this->btxmin = new BOUNDARYTYPE[this->ny];  this->btxmax = new BOUNDARYTYPE[this->ny];
        this->btymin = new BOUNDARYTYPE[this->nx];  this->btymax = new BOUNDARYTYPE[this->nx];

        for (int i = 0; i < this->ny; i++) {
            this->btxmin[i] = _p.btxmin[i]; this->btxmax[i] = _p.btxmax[i];
        }

        for (int i = 0; i < this->nx; i++) {
            this->btymin[i] = _p.btymin[i]; this->btymax[i] = _p.btymax[i];
        }
    }


    template<class T>
    D2Q9<T>::~D2Q9() {
        delete[] this->f0t;     delete[] this->f0tp1;   delete[] this->barrier0;
        delete[] this->f1t;     delete[] this->f1tp1;   delete[] this->barrier1;
        delete[] this->f2t;     delete[] this->f2tp1;   delete[] this->barrier2;
        delete[] this->f3t;     delete[] this->f3tp1;   delete[] this->barrier3;
        delete[] this->f4t;     delete[] this->f4tp1;   delete[] this->barrier4;
        delete[] this->f5t;     delete[] this->f5tp1;   delete[] this->barrier5;
        delete[] this->f6t;     delete[] this->f6tp1;   delete[] this->barrier6;
        delete[] this->f7t;     delete[] this->f7tp1;   delete[] this->barrier7;
        delete[] this->f8t;     delete[] this->f8tp1;   delete[] this->barrier8;

        delete[] this->btxmin;  delete[] this->btxmax;
        delete[] this->btymin;  delete[] this->btymax;
    }


    template<class T>
    void D2Q9<T>::SetBarrier(int _i, int _j, bool _isbarrier) {
        this->barrier0[this->ny*_i + _j] = _isbarrier;
        this->barrier1[this->ny*(_i + 1) + _j] = _isbarrier;
        this->barrier2[this->ny*_i + (_j + 1)] = _isbarrier;
        this->barrier3[this->ny*(_i - 1) + _j] = _isbarrier;
        this->barrier4[this->ny*_i + (_j - 1)] = _isbarrier;
        this->barrier5[this->ny*(_i + 1) + (_j + 1)] = _isbarrier;
        this->barrier6[this->ny*(_i - 1) + (_j + 1)] = _isbarrier;
        this->barrier7[this->ny*(_i - 1) + (_j - 1)] = _isbarrier;
        this->barrier8[this->ny*(_i + 1) + (_j - 1)] = _isbarrier;
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
    bool D2Q9<T>::GetBarrier(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->barrier0[this->ny*_i + _j];
    }
}