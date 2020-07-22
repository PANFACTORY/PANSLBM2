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
        template<class F>
        void SetBoundary(F _f);

        void Stream();
        virtual void UpdateMacro() = 0;
        virtual void Collision() = 0;

        const int nx, ny;

        std::vector<std::vector<bool> > barrier0;
        std::vector<int> btxmin, btxmax, btymin, btymax;

        static const int PERIODIC = 0;  //  0   :   periodic
        static const int INLET = 1;     //  1   :   inlet
        static const int OUTLET = 2;    //  2   :   outlet
        static const int BARRIER = 3;   //  3   :   barrier
        static const int MIRROR = 4;    //  4   :   mirror

protected:
        T dx, dt, t0, t1, t2;

        std::vector<std::vector<T> > f0t, f1t, f2t, f3t, f4t, f5t, f6t, f7t, f8t;
        std::vector<std::vector<T> > f0tp1, f1tp1, f2tp1, f3tp1, f4tp1, f5tp1, f6tp1, f7tp1, f8tp1;
        std::vector<std::vector<bool> > barrier1, barrier2, barrier3, barrier4, barrier5, barrier6, barrier7, barrier8;
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