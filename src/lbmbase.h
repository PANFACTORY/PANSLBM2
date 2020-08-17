//*****************************************************************************
//  Title       :   src/lbmbase.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/12
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>


namespace PANSLBM2 {
    enum BOUNDARYTYPE {
        PERIODIC, INLET, OUTLET, BARRIER, MIRROR,
    };

    template<class T>
    class LBMBASE {
public:
        LBMBASE(int _nx, int _ny);
        ~LBMBASE();

        template<class F>
        void SetBarrier(F _f);
        template<class F>
        void SetBoundary(F _f);

        void Stream();
        virtual void UpdateMacro() = 0;
        virtual void Collision() = 0;

        const int nx, ny;

        bool* barrier0;
        BOUNDARYTYPE* btxmin, btxmax, btymin, btymax;

protected:
        T dx, dt, t0, t1, t2;
        T* f0t, f1t, f2t, f3t, f4t, f5t, f6t, f7t, f8t, f0tp1, f1tp1, f2tp1, f3tp1, f4tp1, f5tp1, f6tp1, f7tp1, f8tp1;
        bool* barrier1, barrier2, barrier3, barrier4, barrier5, barrier6, barrier7, barrier8;
    };


    template<class T>
    LBMBASE<T>::LBMBASE(int _nx, int _ny) : nx(_nx), ny(_ny) {
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
    LBMBASE<T>::~LBMBASE() {
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
    template<class F>
    void LBMBASE<T>::SetBarrier(F _f) {
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
    void LBMBASE<T>::SetBoundary(F _f) {
        //----------Set x boundary----------
        for (int j = 0; j < this->ny; j++) {
            btxmin[j] = _f(0, j);
            btxmax[j] = _f(this->nx - 1, j);
        }

        //----------Set y boundary----------
        for (int i = 0; i < this->nx; i++) {
            btymin[i] = _f(i, 0);
            btymax[i] = _f(i, this->ny - 1);
        }
    }


    template<class T>
    void LBMBASE<T>::Stream() {
        //----------Stream and periodic boundary----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->f0t[i][j] = this->f0tp1[i][j];
                this->f1t[i][j] = (i == 0) ? this->f1tp1[this->nx - 1][j] : this->f1tp1[i - 1][j];
                this->f2t[i][j] = (j == 0) ? this->f2tp1[i][this->ny - 1] : this->f2tp1[i][j - 1];
                this->f3t[i][j] = (i == this->nx - 1) ? this->f3tp1[0][j] : this->f3tp1[i + 1][j];
                this->f4t[i][j] = (j == this->ny - 1) ? this->f4tp1[i][0] : this->f4tp1[i][j + 1];
                this->f5t[i][j] = (i == 0) ? ((j == 0) ? this->f5tp1[this->nx - 1][this->ny - 1] : this->f5tp1[this->nx - 1][j - 1]) : ((j == 0) ? this->f5tp1[i - 1][this->ny - 1] : this->f5tp1[i - 1][j - 1]);
                this->f6t[i][j] = (i == this->nx - 1) ? ((j == 0) ? this->f6tp1[0][this->ny - 1] : this->f6tp1[0][j - 1]) : ((j == 0) ? this->f6tp1[i + 1][this->ny - 1] : this->f6tp1[i + 1][j - 1]);
                this->f7t[i][j] = (i == this->nx - 1) ? ((j == this->ny - 1) ? this->f7tp1[0][0] : this->f7tp1[0][j + 1]) : ((j == this->ny - 1) ? this->f7tp1[i + 1][0] : this->f7tp1[i + 1][j + 1]);
                this->f8t[i][j] = (i == 0) ? ((j == this->ny - 1) ? this->f8tp1[this->nx - 1][0] : this->f8tp1[this->nx - 1][j + 1]) : ((j == this->ny - 1) ? this->f8tp1[i - 1][0] : this->f8tp1[i - 1][j + 1]);
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
                this->f1t[0][j] = this->f1t[1][j];
                this->f5t[0][j] = this->f5t[1][j];
                this->f8t[0][j] = this->f8t[1][j];
            } else if (this->btxmin[j] == BARRIER) {
                this->f1t[0][j] = this->f3tp1[0][j];
                this->f5t[0][j] = this->f7tp1[0][j];
                this->f8t[0][j] = this->f6tp1[0][j];
            } else if (this->btxmin[j] == MIRROR) {
                this->f1t[0][j] = this->f3tp1[0][j];
                this->f5t[0][j] = this->f6tp1[0][j];
                this->f8t[0][j] = this->f7tp1[0][j];
            }

            //.....xmax.....
            if (this->btxmax[j] == OUTLET) {
                this->f3t[this->nx - 1][j] = this->f3t[this->nx - 2][j];
                this->f6t[this->nx - 1][j] = this->f6t[this->nx - 2][j];
                this->f7t[this->nx - 1][j] = this->f7t[this->nx - 2][j];
            } else if (this->btxmax[j] == BARRIER) {
                this->f3t[this->nx - 1][j] = this->f1tp1[this->nx - 1][j];
                this->f6t[this->nx - 1][j] = this->f8tp1[this->nx - 1][j];
                this->f7t[this->nx - 1][j] = this->f5tp1[this->nx - 1][j];
            } else if (this->btxmax[j] == MIRROR) {
                this->f3t[this->nx - 1][j] = this->f1tp1[this->nx - 1][j];
                this->f6t[this->nx - 1][j] = this->f5tp1[this->nx - 1][j];
                this->f7t[this->nx - 1][j] = this->f8tp1[this->nx - 1][j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == OUTLET) {
                this->f2t[i][0] = this->f2t[i][1];
                this->f5t[i][0] = this->f5t[i][1];
                this->f6t[i][0] = this->f6t[i][1];
            } else if (this->btymin[i] == BARRIER) {
                this->f2t[i][0] = this->f4tp1[i][0];
                this->f5t[i][0] = this->f7tp1[i][0];
                this->f6t[i][0] = this->f8tp1[i][0];
            } else if (this->btymin[i] == MIRROR) {
                this->f2t[i][0] = this->f4tp1[i][0];
                this->f5t[i][0] = this->f8tp1[i][0];
                this->f6t[i][0] = this->f7tp1[i][0];
            }

            //.....ymax.....
            if (this->btymax[i] == OUTLET) {
                this->f4t[i][this->ny - 1] = this->f4t[i][this->ny - 2];
                this->f7t[i][this->ny - 1] = this->f7t[i][this->ny - 2];
                this->f8t[i][this->ny - 1] = this->f8t[i][this->ny - 2];
            } else if (this->btymax[i] == BARRIER) {
                this->f4t[i][this->ny - 1] = this->f2tp1[i][this->ny - 1];
                this->f7t[i][this->ny - 1] = this->f5tp1[i][this->ny - 1];
                this->f8t[i][this->ny - 1] = this->f6tp1[i][this->ny - 1];
            } else if (this->btymax[i] == MIRROR) {
                this->f4t[i][this->ny - 1] = this->f2tp1[i][this->ny - 1];
                this->f7t[i][this->ny - 1] = this->f6tp1[i][this->ny - 1];
                this->f8t[i][this->ny - 1] = this->f5tp1[i][this->ny - 1];
            }
        }
    }
}