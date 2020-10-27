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
        ~AD();

        void UpdateMacro();
        void Collision();
        void ExternalForce();

        template<class ...Ts>
        void SetFt(int _i, T _rho, Ts ..._u);
        void SetGt(int _i, T _Temperature);

        T GetRho(int _i) const;
        T GetU(int _d, int _i) const;
        T GetTemperature(int _i) const;

private:
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
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np];
        }
        this->temperature = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = T();
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
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np];
        }
        this->temperature = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = _e.rho[i];
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = _e.u[k][i];
            }
            this->temperature[i] = _e.temperature[i];
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    AD<T, P, Q>::~AD() {
        delete[] this->rho;
        for (int k = 0; k < P<T>::nd; k++) {
            delete[] this->u[k];
        }
        delete[] this->temperature;
    }


    template<class T, template<class>class P, template<class>class Q>
    void AD<T, P, Q>::UpdateMacro() {
        for (int i = 0; i < this->np; i++) {
            this->rho[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = T();
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
                T ciu = T(), uu = T();
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
    template<class ...Ts>
    void AD<T, P, Q>::SetFt(int _i, T _rho, Ts ..._u) {
        assert(P<T>::nd == sizeof...(Ts) && 0 <= _i && _i < this->np);
        T us[P<T>::nd] = { _u... };
        for (int j = 0; j < P<T>::nc; j++) {
            T ciu = T(), uu = T();
            for (int k = 0; k < P<T>::nd; k++) {
                ciu += P<T>::ci[j][k]*us[k];
                uu += us[k]*us[k];
            }
            this->f->ft[j][_i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    void AD<T, P, Q>::SetGt(int _i, T _Temperature) {
        assert(0 <= _i && _i < this->np);
        for (int j = 0; j < Q<T>::nc; j++) {
            T ciu = T();
            for (int k = 0; k < Q<T>::nd; k++) {
                ciu += Q<T>::ci[j][k]*this->u[k][_i];
            }
            this->g->ft[j][_i] = Q<T>::ei[j]*_Temperature*(1.0 + 3.0*ciu);
        }
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetRho(int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->rho[_i];
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetU(int _d, int _i) const {
        assert(0 <= _d && _d < P<T>::nd && 0 <= _i && _i < this->np);
        return this->u[_d][_i];
    }


    template<class T, template<class>class P, template<class>class Q>
    T AD<T, P, Q>::GetTemperature(int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->temperature[_i];
    }
}