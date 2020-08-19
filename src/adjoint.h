//*****************************************************************************
//  Title       :   src/adjoint.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/08/18
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <cassert>


#include "lbm.h"


namespace PANSLBM2 {
    template<class T>
    class AdjointLBM : public LBM<T> {
public:
        AdjointLBM() = delete;
        AdjointLBM(const LBM<T>& _lbm);
        ~AdjointLBM() {};

        void Inlet(T _u, T _v);
        void Stream();
        void UpdateMacro();
        void Collision();
        void ExternalForce();
    };


    template<class T>
    AdjointLBM<T>::AdjointLBM(const LBM<T>& _lbm) : LBM<T>::LBM(_lbm) {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->f0t[i] = T(); this->f0tp1[i] = T();
            this->f1t[i] = T(); this->f1tp1[i] = T();
            this->f2t[i] = T(); this->f2tp1[i] = T();
            this->f3t[i] = T(); this->f3tp1[i] = T();
            this->f4t[i] = T(); this->f4tp1[i] = T();
            this->f5t[i] = T(); this->f5tp1[i] = T();
            this->f6t[i] = T(); this->f6tp1[i] = T();
            this->f7t[i] = T(); this->f7tp1[i] = T();
            this->f8t[i] = T(); this->f8tp1[i] = T();
        }
    }


    template<class T>
    void AdjointLBM<T>::Inlet(T _u, T _v) {

    }


    template<class T>
    void AdjointLBM<T>::Stream() {

    }


    template<class T>
    void AdjointLBM<T>::UpdateMacro() {

    }


    template<class T>
    void AdjointLBM<T>::Collision() {

    }


    template<class T>
    void AdjointLBM<T>::ExternalForce() {
        
    }
}