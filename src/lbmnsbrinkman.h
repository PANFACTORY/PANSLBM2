//*****************************************************************************
//  Title       :   src/lbmnsbrinkman.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/24
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "lbmns.h"


namespace PANSLBM2 {
    template<class T>
    class LBMNSBRINKMAN : public LBMNS<T> {
public:
        LBMNSBRINKMAN(int _nx, int _ny, T _viscosity);
        ~LBMNSBRINKMAN() {};

        template<class F>
        void SetPermeation(F _f);

        void ExternalForce();

//protected:
        std::vector<std::vector<T> > permeation;
    };


    template<class T>
    LBMNSBRINKMAN<T>::LBMNSBRINKMAN(int _nx, int _ny, T _viscosity) : LBMNS<T>(_nx, _ny, _viscosity) {
        this->permeation = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
    }


    template<class T>
    template<class F>
    void LBMNSBRINKMAN<T>::SetPermeation(F _f) {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->permeation[i][j] = _f(i, j);
            }
        }
    }


    template<class T>
    void LBMNSBRINKMAN<T>::ExternalForce() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                T tmprho = this->f0t[i][j] + this->f1t[i][j] + this->f2t[i][j] + this->f3t[i][j] + this->f4t[i][j] + this->f5t[i][j] + this->f6t[i][j] + this->f7t[i][j] + this->f8t[i][j];
                T tmpu = (this->f1t[i][j] - this->f3t[i][j] + this->f5t[i][j] - this->f6t[i][j] - this->f7t[i][j] + this->f8t[i][j])/tmprho;
                T tmpv = (this->f2t[i][j] - this->f4t[i][j] + this->f5t[i][j] + this->f6t[i][j] - this->f7t[i][j] - this->f8t[i][j])/tmprho;

                T dxt1alpha = -3.0*this->dx*this->t1*this->permeation[i][j];
                T dxt2alpha = -3.0*this->dx*this->t2*this->permeation[i][j];

                this->f1t[i][j] += dxt1alpha*tmpu;
                this->f2t[i][j] += dxt1alpha*tmpv;
                this->f3t[i][j] += dxt1alpha*(-tmpu);
                this->f4t[i][j] += dxt1alpha*(-tmpv);
                this->f5t[i][j] += dxt2alpha*(tmpu + tmpv);
                this->f6t[i][j] += dxt2alpha*(-tmpu + tmpv);
                this->f7t[i][j] += dxt2alpha*(-tmpu - tmpv);
                this->f8t[i][j] += dxt2alpha*(tmpu - tmpv);
            }
        }
    }
}