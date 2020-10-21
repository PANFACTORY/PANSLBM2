//*****************************************************************************
//  Title       :   src/navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P>
    class NS {
public:
        NS() = delete;
        NS(P<T>* _f, T _viscosity);
        NS(const NS<T, P>& _e);
        virtual ~NS();

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _d, int _i, int _j) const;
        
protected:
        T omega;
        P<T>* f;
        T *rho, *u[P<T>::nd];
    };


    template<class T, template<class>class P>
    NS<T, P>::NS(P<T>* _f, T _viscosity) {
        assert(0 < _f->nx && 0 < _f->ny);
        this->f = _f;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/(this->f->dx*this->f->dx) + 0.5);

        this->rho = new T[this->f->np];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->f->np];
        }

        for (int i = 0; i < this->f->np; i++) {
            this->rho[i] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = T();
            }
        }
    }


    template<class T, template<class>class P>
    NS<T, P>::NS(const NS<T, P>& _e) {
        this->f = _e.f;
        this->omega = _e.omega;

        this->rho = new T[this->f->np];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->f->np];
        }

        for (int i = 0; i < this->f->np; i++) {
            this->rho[i] = _e.rho[i];
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = _e.u[j][i];
            }
        }
    }


    template<class T, template<class>class P>
    NS<T, P>::~NS() {
        delete[] this->rho;
        for (int i = 0; i < P<T>::nd; i++) {
            delete[] this->u[i];
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::UpdateMacro() {
        for (int i = 0; i < this->f->np; i++) {
            this->rho[i] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                this->rho[i] += this->f->ft[j][i];
                for (int k = 0; k < P<T>::nd; k++) {
                    this->u[k][i] += P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] /= this->rho[i];
            }
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::Collision() {
        for (int i = 0; i < this->f->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][i];
                    uu += this->u[k][i]*this->u[k][i];
                }
                T fieq = P<T>::ei[j]*this->rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*fieq;
            }
        }
    }


    template<class T, template<class>class P>
    void NS<T, P>::ExternalForce() {
        for (int i = 0; i < this->f->np; i++) {
            
        }
    }


    template<class T, template<class>class P>
    T NS<T, P>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->rho[this->f->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NS<T, P>::GetU(int _d, int _i, int _j) const {
        assert(0 < _d < P<T>::nd && 0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->u[_d][this->f->ny*_i + _j];
    }
}