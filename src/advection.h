//*****************************************************************************
//  Title       :   src/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P>
    class AD {
public:
        AD() = delete;
        AD(P<T>* _f, T _diffusivity);
        AD(const AD<T>& _e);
        virtual ~AD();

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void SetTemperature(int _i, int _j, int _temperature);
        virtual void SetFlux(int _i, int _j, T _ux, T _uy, T _q);
        virtual void ExternalForce();

        virtual T GetTemperature(int _i, int _j) const;

        const int nx, ny;

protected:
        T omega;
        P<T>* f;
        T *temperature;
    };


    template<class T, template<class>class P>
    AD<T, P>::AD(P<T>* _f, T _diffusivity) : nx(_f->nx), ny(_f->ny) {
        assert(0 < _f->nx && 0 < _f->ny);
        this->f = _f;
        this->omega = 1.0/(3.0*_diffusivity*this->dt/pow(this->dx, 2.0) + 0.5);

        this->temperature = new T[this->nx*this->ny];
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->temperature[i] = T();
        }
    }


    template<class T, template<class>class P>
    AD<T, P>::AD(const AD<T>& _e) : nx(_e.nx), ny(_e.ny) {
        this->f = _e.f;
        this->omega = _e.omega;

        this->temperature = new T[this->nx*this->ny];
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->temperature[i] = _e.temperature[i];
        }
    }


    template<class T, template<class>class P>
    AD<T, P>::~AD() {
        delete[] this->temperature;
    }


    template<class T, template<class>class P>
    void AD<T, P>::UpdateMacro() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            this->temperature[i] = T();
            for (int j = 0; j < P<T>::nc; j++) {
                this->temperature[i] += this->f->ft[j][i];
            }
        }
    }


    template<class T, template<class>class P>
    void AD<T, P>::Collision() {
        for (int i = 0; i < this->nx*this->ny; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][i];
                }
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*P<T>::ei[j]*this->temperature[i]*(1.0 + 3.0*ciu);
            }
        }
    }


    template<class T, template<class>class P>
    void AD<T, P>::SetTemperature(int _i, int _j, int _temperature) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T temperature0 = 6.0*(_temperature - this->f->ft[0][ij] - this->f->ft[2][ij] - this->f->ft[3][ij] - this->f->ft[4][ij] - this->f->ft[6][ij] - this->f->ft[7][ij]);
            this->f->ft[1][ij] = temperature0/9.0;
            this->f->ft[5][ij] = temperature0/36.0;
            this->f->ft[8][ij] = temperature0/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(_temperature - this->f->ft[0][ij] - this->f->ft[1][ij] - this->f->ft[2][ij] - this->f->ft[4][ij] - this->f->ft[5][ij] - this->f->ft[8][ij]);
            this->f->ft[3][ij] = temperature0/9.0;
            this->f->ft[6][ij] = temperature0/36.0;
            this->f->ft[7][ij] = temperature0/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_temperature - this->f->ft[0][ij] - this->f->ft[1][ij] - this->f->ft[3][ij] - this->f->ft[4][ij] - this->f->ft[7][ij] - this->f->ft[8][ij]);
            this->f->ft[2][ij] = temperature0/9.0;
            this->f->ft[5][ij] = temperature0/36.0;
            this->f->ft[6][ij] = temperature0/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(_temperature - this->f->ft[0][ij] - this->f->ft[1][ij] - this->f->ft[2][ij] - this->f->ft[3][ij] - this->f->ft[5][ij] - this->f->ft[6][ij]);
            this->f->ft[4][ij] = temperature0/9.0;
            this->f->ft[7][ij] = temperature0/36.0;
            this->f->ft[8][ij] = temperature0/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }
    
    
    template<class T, template<class>class P>
    void AD<T, P>::SetFlux(int _i, int _j, T _ux, T _uy, T _q) {
        int ij = this->ny*_i + _j;
        if (_i == 0) {
            T temperature0 = 6.0*(_q + this->f->ft[3][ij] + this->f->ft[6][ij] + this->f->ft[7][ij])/(1.0 - 3.0*_ux);
            this->f->ft[1][ij] = temperature0*(1.0 + 3.0*_ux)/9.0;
            this->f->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->f->ft[8][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
        } else if (_i == this->nx - 1) {
            T temperature0 = 6.0*(-_q + this->f->ft[1][ij] + this->f->ft[5][ij] + this->f->ft[8][ij])/(1.0 + 3.0*_ux);
            this->f->ft[3][ij] = temperature0*(1.0 - 3.0*_ux)/9.0;
            this->f->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
            this->f->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy/36.0;
        } else if (_j == 0) {
            T temperature0 = 6.0*(_q + this->f->ft[4][ij] + this->f->ft[7][ij] + this->f->ft[8][ij])/(1.0 - 3.0*_uy);
            this->f->ft[2][ij] = temperature0*(1.0 + 3.0*_uy)/9.0;
            this->f->ft[5][ij] = temperature0*(1.0 + 3.0*_ux + 3.0*_uy)/36.0;
            this->f->ft[6][ij] = temperature0*(1.0 - 3.0*_ux + 3.0*_uy)/36.0;
        } else if (_j == this->ny - 1) {
            T temperature0 = 6.0*(-_q + this->f->ft[2][ij] + this->f->ft[5][ij] + this->f->ft[6][ij])/(1.0 + 3.0*_uy);
            this->f->ft[4][ij] = temperature0*(1.0 - 3.0*_uy)/9.0;
            this->f->ft[7][ij] = temperature0*(1.0 - 3.0*_ux - 3.0*_uy)/36.0;
            this->f->ft[8][ij] = temperature0*(1.0 + 3.0*_ux - 3.0*_uy)/36.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T, template<class>class P>
    void AD<T, P>::ExternalForce() {
        /*for (int i = 0; i < this->nx*this->ny; i++) {
            T temperature0 = this->g.f0t[i] + this->g.f1t[i] + this->g.f2t[i] + this->g.f3t[i] + this->g.f4t[i] + this->g.f5t[i] + this->g.f6t[i] + this->g.f7t[i] + this->g.f8t[i];
            
            T rhog = 1.6e-4*(temperature0 - 1.5);

            T dxt0alpha = 3.0*this->dx*this->t0;
            T dxt1alpha = 3.0*this->dx*this->t1;
            T dxt2alpha = 3.0*this->dx*this->t2;
    
            this->f.f2t[i] += dxt1alpha*rhog;
            this->f.f4t[i] += dxt1alpha*-rhog;
            this->f.f5t[i] += dxt2alpha*rhog;
            this->f.f6t[i] += dxt2alpha*rhog;
            this->f.f7t[i] += dxt2alpha*-rhog;
            this->f.f8t[i] += dxt2alpha*-rhog;
        }*/
    }


    template<class T, template<class>class P>
    T AD<T, P>::GetTemperature(int _i, int _j) const {
        assert(0 <= _i && _i < this->nx && 0 <= _j && _j < this->ny);
        return this->temperature[this->ny*_i + _j];
    }
}