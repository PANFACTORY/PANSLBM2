//*****************************************************************************
//  Title       :   src/lbmbase.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/12
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <vector>
#include <array>


namespace PANSLBM2 {
    template<class T>
    class LBMBASE {
public:
        LBMBASE(int _nx, int _ny);
        ~LBMBASE() {};

        template<class F>
        void SetBarrier(F _f);

        void SetBoundary(int _btxmin, int _btxmax, int _btymin, int _btymax);

        void Stream();
        virtual void UpdateMacro() = 0;
        virtual void Collision() = 0;

        const int nx, ny;

        std::vector<std::vector<bool> > barrier0;

        std::vector<int> btxmin;    //  0 : periodic, 1 : inlet, 2 : outlet, 3 : barrier, 4 : mirror
        std::vector<int> btxmax;
        std::vector<int> btymin;
        std::vector<int> btymax;


protected:
        T dx, dt, t0, t1, t2;

        std::vector<std::vector<T> > f0t;
        std::vector<std::vector<T> > f1t;
        std::vector<std::vector<T> > f2t;
        std::vector<std::vector<T> > f3t;
        std::vector<std::vector<T> > f4t;
        std::vector<std::vector<T> > f5t;
        std::vector<std::vector<T> > f6t;
        std::vector<std::vector<T> > f7t;
        std::vector<std::vector<T> > f8t;

        std::vector<std::vector<T> > f0tp1;
        std::vector<std::vector<T> > f1tp1;
        std::vector<std::vector<T> > f2tp1;
        std::vector<std::vector<T> > f3tp1;
        std::vector<std::vector<T> > f4tp1;
        std::vector<std::vector<T> > f5tp1;
        std::vector<std::vector<T> > f6tp1;
        std::vector<std::vector<T> > f7tp1;
        std::vector<std::vector<T> > f8tp1;

        std::vector<std::vector<bool> > barrier1;
        std::vector<std::vector<bool> > barrier2;
        std::vector<std::vector<bool> > barrier3;
        std::vector<std::vector<bool> > barrier4;
        std::vector<std::vector<bool> > barrier5;
        std::vector<std::vector<bool> > barrier6;
        std::vector<std::vector<bool> > barrier7;
        std::vector<std::vector<bool> > barrier8;
    };


    template<class T>
    LBMBASE<T>::LBMBASE(int _nx, int _ny) : nx(_nx), ny(_ny) {
        this->dx = 1.0;
        this->dt = 1.0;
        this->t0 = 4.0/9.0;
        this->t1 = 1.0/9.0;
        this->t2 = 1.0/36.0;

        this->f0t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t0));
        this->f1t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f2t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f3t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f4t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f5t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f6t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f7t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f8t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));

        this->f0tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f1tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f2tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f3tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f4tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f5tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f6tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f7tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f8tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
    
        this->barrier0 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier1 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier2 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier3 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier4 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier5 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier6 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier7 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));
        this->barrier8 = std::vector<std::vector<bool> >(this->nx, std::vector<bool>(this->ny, false));

        this->btxmin = std::vector<int>(this->ny, 0);
        this->btxmax = std::vector<int>(this->ny, 0);
        this->btymin = std::vector<int>(this->nx, 0);
        this->btymax = std::vector<int>(this->nx, 0);
    }


    template<class T>
    template<class F>
    void LBMBASE<T>::SetBarrier(F _f) {
        //----------Set barrier0----------
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                if (_f(i, j)) {
                    this->barrier0[i][j] = true;
                }
            }
        }

        //----------Set barrier1~8----------
        for (int i = 1; i < this->nx - 1; i++) {
            for (int j = 1; j < this->ny - 1; j++) {
                this->barrier1[i][j] = this->barrier0[i - 1][j];
                this->barrier2[i][j] = this->barrier0[i][j - 1];
                this->barrier3[i][j] = this->barrier0[i + 1][j];
                this->barrier4[i][j] = this->barrier0[i][j + 1];
                this->barrier5[i][j] = this->barrier0[i - 1][j - 1];
                this->barrier6[i][j] = this->barrier0[i + 1][j - 1];
                this->barrier7[i][j] = this->barrier0[i + 1][j + 1];
                this->barrier8[i][j] = this->barrier0[i - 1][j + 1];
            }
        }
    }


    template<class T>
    void LBMBASE<T>::SetBoundary(int _btxmin, int _btxmax, int _btymin, int _btymax) {
        for (int j = 0; j < this->ny; j++) {
            btxmin[j] = _btxmin;
            btxmax[j] = _btxmax;
        }

        for (int i = 0; i < this->nx; i++) {
            btymin[i] = _btymin;
            btymax[i] = _btymax;
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
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                if (this->barrier1[i][j]) {
                    this->f1t[i][j] = this->f3tp1[i][j];
                }
                if (this->barrier2[i][j]) {
                    this->f2t[i][j] = this->f4tp1[i][j];
                }
                if (this->barrier3[i][j]) {
                    this->f3t[i][j] = this->f1tp1[i][j];
                }
                if (this->barrier4[i][j]) {
                    this->f4t[i][j] = this->f2tp1[i][j];
                }
                if (this->barrier5[i][j]) {
                    this->f5t[i][j] = this->f7tp1[i][j];
                }
                if (this->barrier6[i][j]) {
                    this->f6t[i][j] = this->f8tp1[i][j];
                }
                if (this->barrier7[i][j]) {
                    this->f7t[i][j] = this->f5tp1[i][j];
                }
                if (this->barrier8[i][j]) {
                    this->f8t[i][j] = this->f6tp1[i][j];
                }
            }
        }

        //----------boundary (Bouns-Back, Outlet and Mirror)----------
        for (int j = 0; j < this->ny; j++) {
            //.....xmin.....
            if (this->btxmin[j] == 2) {
                this->f1t[0][j] = this->f1t[1][j];
                this->f5t[0][j] = this->f5t[1][j];
                this->f8t[0][j] = this->f8t[1][j];
            } else if (this->btxmin[j] == 3) {
                this->f1t[0][j] = this->f3tp1[0][j];
                this->f5t[0][j] = this->f7tp1[0][j];
                this->f8t[0][j] = this->f6tp1[0][j];
            } else if (this->btxmin[j] == 4) {
                this->f1t[0][j] = this->f3tp1[0][j];
                this->f5t[0][j] = this->f6tp1[0][j];
                this->f8t[0][j] = this->f7tp1[0][j];
            }

            //.....xmax.....
            if (this->btxmax[j] == 2) {
                this->f3t[this->nx - 1][j] = this->f3t[this->nx - 2][j];
                this->f6t[this->nx - 1][j] = this->f6t[this->nx - 2][j];
                this->f7t[this->nx - 1][j] = this->f7t[this->nx - 2][j];
            } else if (this->btxmax[j] == 3) {
                this->f3t[this->nx - 1][j] = this->f1tp1[this->nx - 1][j];
                this->f6t[this->nx - 1][j] = this->f8tp1[this->nx - 1][j];
                this->f7t[this->nx - 1][j] = this->f5tp1[this->nx - 1][j];
            } else if (this->btxmax[j] == 4) {
                this->f3t[this->nx - 1][j] = this->f1tp1[this->nx - 1][j];
                this->f6t[this->nx - 1][j] = this->f5tp1[this->nx - 1][j];
                this->f7t[this->nx - 1][j] = this->f8tp1[this->nx - 1][j];
            }
        }

        for (int i = 0; i < this->nx; i++) {
            //.....ymin.....
            if (this->btymin[i] == 2) {
                this->f2t[i][0] = this->f2t[i][1];
                this->f5t[i][0] = this->f5t[i][1];
                this->f6t[i][0] = this->f6t[i][1];
            } else if (this->btymin[i] == 3) {
                this->f2t[i][0] = this->f4tp1[i][0];
                this->f5t[i][0] = this->f7tp1[i][0];
                this->f6t[i][0] = this->f8tp1[i][0];
            } else if (this->btymin[i] == 4) {
                this->f2t[i][0] = this->f4tp1[i][0];
                this->f5t[i][0] = this->f8tp1[i][0];
                this->f6t[i][0] = this->f7tp1[i][0];
            }

            //.....ymax.....
            if (this->btymax[i] == 2) {
                this->f4t[i][this->ny - 1] = this->f4t[i][this->ny - 2];
                this->f7t[i][this->ny - 1] = this->f7t[i][this->ny - 2];
                this->f8t[i][this->ny - 1] = this->f8t[i][this->ny - 2];
            } else if (this->btymax[i] == 3) {
                this->f4t[i][this->ny - 1] = this->f2tp1[i][this->ny - 1];
                this->f7t[i][this->ny - 1] = this->f5tp1[i][this->ny - 1];
                this->f8t[i][this->ny - 1] = this->f6tp1[i][this->ny - 1];
            } else if (this->btymax[i] == 4) {
                this->f4t[i][this->ny - 1] = this->f2tp1[i][this->ny - 1];
                this->f7t[i][this->ny - 1] = this->f6tp1[i][this->ny - 1];
                this->f8t[i][this->ny - 1] = this->f5tp1[i][this->ny - 1];
            }
        }
    }
}