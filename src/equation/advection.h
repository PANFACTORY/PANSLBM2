//*****************************************************************************
//  Title       :   src/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P, template<class>class Q>
    class AD {
public:
        AD() = delete;
        AD(P<T>* _f, Q<T>* _g, T _viscosity, T _diffusivity);
        AD(const AD<T, P, Q>& _e);
        virtual ~AD();

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        virtual T GetRho(int _i, int _j) const;
        virtual T GetU(int _d, int _i, int _j) const;
        virtual T GetTemperature(int _i, int _j) const;

protected:
        const int np;           //  np : number of particle
        T omegaf, omegag;
        P<T> *f;
        Q<T> *g;
        T *rho, *u[P<T>::nd], *temperature;
    };


    template<class T, template<class>class P, template<class>class Q>
    AD<T, P, Q>::AD(P<T>* _f, Q<T>* _g, T _viscosity, T _diffusivity) : np(_f->np) {
        // ここでfとgの整合性を確認すること
        this->f = _f;
        this->omegaf = 1.0/(3.0*_viscosity*this->f->dt/(this->f->dx*this->f->dx) + 0.5);
        this->g = _g;
        this->omegag = 1.0/(3.0*_diffusivity*this->g->dt/(this->g->dx*this->g->dx) + 0.5);
        
        this->rho = new T[this->np];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->np];
        }
        this->temperature = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = T();
            }
            this->temperature[i] = T();
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    AD<T, P, Q>::AD(const AD<T, P, Q>& _e) : np(_e.np) {
        this->f = _e.f;
        this->omegaf = _e.omegaf;
        this->g = _e.g;
        this->omegag = _e.omegag;

        this->rho = new T[this->np];
        for (int i = 0; i < P<T>::nd; i++) {
            this->u[i] = new T[this->np];
        }
        this->temperature = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = _e.rho[i];
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][i] = _e.u[j][i];
            }
            this->temperature[i] = _e.temperature[i];
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    AD<T, P, Q>::~AD() {
        delete[] this->rho;
        for (int i = 0; i < P<T>::nd; i++) {
            delete[] this->u[i];
        }
        delete[] this->temperature;
    }


    template<class T, template<class>class P, template<class>class Q>
    void AD<T, P, Q>::UpdateMacro() {
        for (int i = 0; i < this->np; i++) {
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
            this->temperature[i] = T();
            for (int j = 0; j < Q<T>::nc; j++) {
                this->temperature[i] += this->g->ft[j][i];
            }
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    void AD<T, P, Q>::Collision() {
        for (int i = 0; i < this->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][i];
                    uu += this->u[k][i]*this->u[k][i];
                }
                T feq = P<T>::ei[j]*this->rho[i]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                this->f->ftp1[j][i] = (1.0 - this->omegaf)*this->f->ft[j][i] + this->omegaf*feq;
            }
            for (int j = 0; j < Q<T>::nc; j++) {
                T ciu = T();
                for (int k = 0; k < Q<T>::nd; k++) {
                    ciu += Q<T>::ci[j][k]*this->u[k][i];
                }
                T geq = Q<T>::ei[j]*this->temperature[i]*(1.0 + 3.0*ciu);
                this->g->ftp1[j][i] = (1.0 - this->omegag)*this->g->ft[j][i] + this->omegag*geq;
            }
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    void AD<T, P, Q>::ExternalForce() {
        for (int i = 0; i < this->np; i++) {
            T temperature0 = T();
            for (int j = 0; j < Q<T>::nc; j++) {
                temperature0 += this->g->ft[j][i];
            }            
            T rhog = 1.6e-5*(temperature0 - 1.5);
            for (int j = 0; j < P<T>::nc; j++) {
                this->f->ft[j][i] += 3.0*this->f->dx*P<T>::ei[j]*Q<T>::ci[j][1]*rhog;
            }
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetRho(int _i, int _j) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->rho[this->f->ny*_i + _j];
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetU(int _d, int _i, int _j) const {
        assert(0 < _d < P<T>::nd && 0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->u[_d][this->f->ny*_i + _j];
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetTemperature(int _i, int _j) const {
        assert(0 <= _i && _i < this->g->nx && 0 <= _j && _j < this->g->ny);
        return this->temperature[this->g->ny*_i + _j];
    }
}