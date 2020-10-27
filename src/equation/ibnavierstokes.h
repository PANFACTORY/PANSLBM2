//*****************************************************************************
//  Title       :   src/ibnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P>
    class IBNS {
public:
        IBNS() = delete;
        IBNS(P<T>* _f, T _viscosity);
        IBNS(const NS<T, P>& _e);
        ~IBNS();

        void UpdateMacro();
        void Collision();
        void ExternalForce();

        template<class ...Ts>
        void SetFt(int _i, T _rho, Ts ..._u);
        template<class ...Ts>
        void SetXk(int _k, Ts ..._xk);
        template<class ...Ts>
        void SetUk(int _k, Ts ..._uk);

        T GetRho(int _i) const;
        T GetU(int _d, int _i) const;
        
private:
        const int np;                   //  np : number of particle
        T omega;
        P<T>* f;
        T *rho, *u[P<T>::nd];           //  Macroscopic value
        T *Xk[P<T>::nd], *Uk[P<T>::nd]; //  Boundary point coordinate and velocity
    };


    template<class T, template<class>class P>
    IBNS<T, P>::IBNS(P<T>* _f, T _viscosity) : np(_f->np) {
        assert(0 < _f->np);
        this->f = _f;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/(this->f->dx*this->f->dx) + 0.5);

        this->rho = new T[this->np];
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np];
        }

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = T();
            }
        }
    }


    template<class T, template<class>class P>
    IBNS<T, P>::IBNS(const NS<T, P>& _e) : np(_e.np) {
        this->f = _e.f;
        this->omega = _e.omega;

        this->rho = new T[this->np];
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np];
        }

        for (int i = 0; i < this->np; i++) {
            this->rho[i] = _e.rho[i];
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = _e.u[k][i];
            }
        }
    }


    template<class T, template<class>class P>
    IBNS<T, P>::~IBNS() {
        delete[] this->rho;
        for (int k = 0; k < P<T>::nd; k++) {
            delete[] this->u[k];
        }
    }


    template<class T, template<class>class P>
    void IBNS<T, P>::UpdateMacro() {
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
        }
    }


    template<class T, template<class>class P>
    void IBNS<T, P>::Collision() {
        for (int i = 0; i < this->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T(), uu = T();
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
    void IBNS<T, P>::ExternalForce() {
        for (int i = 0; i < this->np; i++) {
            
        }
    }


    template<class T, template<class>class P>
    template<class ...Ts>
    void IBNS<T, P>::SetFt(int _i, T _rho, Ts ..._u) {
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


    template<class T, template<class>class P>
    template<class ...Ts>
    void IBNS<T, P>::SetXk(int _k, Ts ..._xk) {
        assert(P<T>::nd == sizeof...(Ts));
        
    }


    template<class T, template<class>class P>
    template<class ...Ts>
    void IBNS<T, P>::SetUk(int _k, Ts ..._uk) {
        assert(P<T>::nd == sizeof...(Ts));
        
    }


    template<class T, template<class>class P>
    T IBNS<T, P>::GetRho(int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->rho[_i];
    }


    template<class T, template<class>class P>
    T IBNS<T, P>::GetU(int _d, int _i) const {
        assert(0 <= _d && _d < P<T>::nd && 0 <= _i && _i < this->np);
        return this->u[_d][_i];
    }
}