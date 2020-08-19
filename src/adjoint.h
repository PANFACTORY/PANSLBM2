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
        ~AdjointLBM();

private:
        
    };


    template<class T>
    AdjointLBM<T>::AdjointLBM(const LBM<T>& _lbm) : LBM<T>::LBM(_lbm) {
std::cout << _lbm.dt << std::endl;
    }


    template<class T>
    AdjointLBM<T>::~AdjointLBM() {

    }
}