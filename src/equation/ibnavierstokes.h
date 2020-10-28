//*****************************************************************************
//  Title       :   src/ibnavierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/20
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
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
        void SetXb(int _k, Ts ..._xb);
        template<class ...Ts>
        void SetUb(int _k, Ts ..._ub);

        T GetRho(int _i) const;
        T GetU(int _d, int _i) const;
        
private:
        const int np, nb;                               //  np : number of particle, nb : number of boundary point
        T omega;
        P<T>* f;
        T *rho, *u[P<T>::nd], *G[P<T>::nd];             //  Macroscopic value
        T *Xb[P<T>::nd], *Ub[P<T>::nd], *Gb[P<T>::nd];  //  Boundary point coordinate and velocity
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
        //  Step0：境界点上の体積力を初期化する
        for (int k = 0; k < this->nb; k++) {
            for (int d = 0; d < P<T>::nd; d++) {
                this->Gb[d][k] = T();
            }
        }

        for (int l = 0; l < 5; l++) {
            //  格子点 → 境界点
            for (int k = 0; k < this->nb; k ++) {
                //  Step3：格子点上の流速 → 境界点上の流速
                T vb[P<T>::nd];
                for (int d = 0; d < P<T>::nd; d++) {
                    vb[d] = T();
                }
                for (int i = 0; i < this->f->nx; i++) {
                    for (int j = 0; j < this->f->ny; j++) {
                        int ij = this->f->ny*i + j;
                        T w = 1.0, rx = fabs(this->Xb[0][k] - i), ry = fabs(this->Xb[1][k] - j);
                        w *= rx < 1.0 ? (3.0 - 2.0*rx + sqrt(1.0 + 4.0*rx - 4.0*rx*rx))/8.0 : (1.0 <= rx && rx < 2.0 ? (5.0 - 2.0*rx - sqrt(-7.0 + 12.0*rx - 4.0*rx*rx))/8.0 : T());
                        w *= ry < 1.0 ? (3.0 - 2.0*ry + sqrt(1.0 + 4.0*ry - 4.0*ry*ry))/8.0 : (1.0 <= ry && ry < 2.0 ? (5.0 - 2.0*ry - sqrt(-7.0 + 12.0*ry - 4.0*ry*ry))/8.0 : T());
                        vb0 += this->u[0][ij]*w*pow(this->f->dx, P<T>::nd);
                        vb1 += this->u[1][ij]*w*pow(this->f->dx, P<T>::nd);
                    }
                }

                //  Step4：境界点上の流速 → 境界点上の体積力
                for (int d = 0; d < P<T>::nd; d++) {
                    this->Gb[d][k] += this->f->dx*(this->Ub[d][k] - vb[d]);
                }
            }

            //  境界点 → 格子点
            for (int i = 0; i < this->f->nx; i++) {
                for (int j = 0; j < this->f->ny; j++) {
                    int ij = this->f->ny*i + j;

                    //  Step1：境界点上の体積力 → 格子点上の体積力
                    for (int d = 0; d < P<T>::nd; d++) {
                        this->G[d][ij] = T();
                    }
                    for (int k = 0; k < this->nb; k++) {
                        T w = 1.0, rx = fabs(this->Xb[0][k] - i), ry = fabs(this->Xb[1][k] - j);
                        w *= rx < 1.0 ? (3.0 - 2.0*rx + sqrt(1.0 + 4.0*rx - 4.0*rx*rx))/8.0 : (1.0 <= rx && rx < 2.0 ? (5.0 - 2.0*rx - sqrt(-7.0 + 12.0*rx - 4.0*rx*rx))/8.0 : T());
                        w *= ry < 1.0 ? (3.0 - 2.0*ry + sqrt(1.0 + 4.0*ry - 4.0*ry*ry))/8.0 : (1.0 <= ry && ry < 2.0 ? (5.0 - 2.0*ry - sqrt(-7.0 + 12.0*ry - 4.0*ry*ry))/8.0 : T());
                        this->G[0][ij] += this->Gb[0][k]*w*pow(this->f->dx, P<T>::nd);
                        this->G[1][ij] += this->Gb[1][k]*w*pow(this->f->dx, P<T>::nd);
                    }

                    //  Step2：格子点上の体積力 → 格子点上の流速
                    for (int d = 0; d < P<T>::nd; d++) {
                        this->u[d][ij] += this->f->dx*this->G[d][ij];
                    }
                }
            }
        }
        
        //  Step5：速度分布関数を更新する
        for (int i = 0; i < this->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciG = T();
                for (int d = 0; d < P<T>::nd; d++) {
                    ciG += P<T>::ci[j][d]*this->G[d][i];
                }
                this->f->ft[j] += 3.0*this->f->dx*P<T>::ei[j]*ciG;
            }
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
    void IBNS<T, P>::SetXb(int _k, Ts ..._xb) {
        assert(P<T>::nd == sizeof...(Ts));
        
    }


    template<class T, template<class>class P>
    template<class ...Ts>
    void IBNS<T, P>::SetUb(int _k, Ts ..._ub) {
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