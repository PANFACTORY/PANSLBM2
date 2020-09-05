//*****************************************************************************
//  Title       :   src/thermal.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/09/03
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include "lbm.h"


namespace PANSLBM2 {
    template<class T>
    class ThermalLBM : public LBM<T> {
public:
        ThermalLBM() = delete;
        ThermalLBM(int _nx, int _ny, T _viscosity);
        ~ThermalLBM();

        void Inlet(T _temperature);
        void UpdateMacro();
        void Collision();
        void ExternalForce();

        void CaptureMacro(const LBM<T>& _lbm);

        T GetTemperature(int _i, int _j) const;

private:
        T *temperature;
    };


    template<class T>
    ThermalLBM<T>::ThermalLBM(int _nx, int _ny, T _viscosity) : LBM<T>::LBM(_nx, _ny, _viscosity) {
        this->temperature = new T[this->nx*this->ny]; 

        for (int i = 0; i < this->nx*this->ny; i++) {
            this->temperature[i] = T();
        }
    }


    template<class T>
    ThermalLBM<T>::~ThermalLBM() {
        delete[] this->temperature;
    }


    template<class T>
    void ThermalLBM<T>::Inlet(T _temperature) {
        for (int i = 0; i < this->nx; i++) {
            T temperature0 = 6.0*(_temperature - this->f0t[this->ny*i] - this->f1t[this->ny*i] - this->f3t[this->ny*i] - this->f4t[this->ny*i] - this->f7t[this->ny*i] - this->f8t[this->ny*i]);
            this->f2t[this->ny*i] = temperature0/9.0;
            this->f5t[this->ny*i] = temperature0/36.0;
            this->f6t[this->ny*i] = temperature0/36.0;
        }
    }


    template<class T>
    void ThermalLBM<T>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->temperature[i] = this->f0t[i] + this->f1t[i] + this->f2t[i] + this->f3t[i] + this->f4t[i] + this->f5t[i] + this->f6t[i] + this->f7t[i] + this->f8t[i];
        }
    }


    template<class T>
    void ThermalLBM<T>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            T g0eq = this->t0*this->temperature[i];
            T g1eq = this->t1*this->temperature[i]*(1.0 + 3.0*this->u[i]);
            T g2eq = this->t1*this->temperature[i]*(1.0 + 3.0*this->v[i]);
            T g3eq = this->t1*this->temperature[i]*(1.0 - 3.0*this->u[i]);
            T g4eq = this->t1*this->temperature[i]*(1.0 - 3.0*this->v[i]);
            T g5eq = this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] + 3.0*this->v[i]);
            T g6eq = this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] + 3.0*this->v[i]);
            T g7eq = this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] - 3.0*this->v[i]);
            T g8eq = this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] - 3.0*this->v[i]);

            this->f0tp1[i] = (1.0 - this->omega)*this->f0t[i] + this->omega*this->t0*this->temperature[i];
            this->f1tp1[i] = (1.0 - this->omega)*this->f1t[i] + this->omega*this->t1*this->temperature[i]*(1.0 + 3.0*this->u[i]);
            this->f2tp1[i] = (1.0 - this->omega)*this->f2t[i] + this->omega*this->t1*this->temperature[i]*(1.0 + 3.0*this->v[i]);
            this->f3tp1[i] = (1.0 - this->omega)*this->f3t[i] + this->omega*this->t1*this->temperature[i]*(1.0 - 3.0*this->u[i]);
            this->f4tp1[i] = (1.0 - this->omega)*this->f4t[i] + this->omega*this->t1*this->temperature[i]*(1.0 - 3.0*this->v[i]);
            this->f5tp1[i] = (1.0 - this->omega)*this->f5t[i] + this->omega*this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] + 3.0*this->v[i]);
            this->f6tp1[i] = (1.0 - this->omega)*this->f6t[i] + this->omega*this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] + 3.0*this->v[i]);
            this->f7tp1[i] = (1.0 - this->omega)*this->f7t[i] + this->omega*this->t2*this->temperature[i]*(1.0 - 3.0*this->u[i] - 3.0*this->v[i]);
            this->f8tp1[i] = (1.0 - this->omega)*this->f8t[i] + this->omega*this->t2*this->temperature[i]*(1.0 + 3.0*this->u[i] - 3.0*this->v[i]);
        }
    }


    template<class T>
    void ThermalLBM<T>::ExternalForce() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            
        }
    }


    template<class T>
    void ThermalLBM<T>::CaptureMacro(const LBM<T>& _lbm) {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->rho[this->ny*i + j] = _lbm.GetRho(i, j);  this->u[this->ny*i + j] = _lbm.GetU(i, j);  this->v[this->ny*i + j] = _lbm.GetV(i, j);
            }
        }
    }


    template<class T>
    T ThermalLBM<T>::GetTemperature(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->temperature[this->ny*_i + _j];
    }
}