//*****************************************************************************
//  Title       :   src/nsadjoint.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/17
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#pragma once
#include <cmath>
#include <cassert>


namespace PANSLBM2 {
    template<class T, template<class>class P>
    class NSAdjoint {
public:
        NSAdjoint() = delete;
        NSAdjoint(P<T>* _f, T _viscosity, int _nt);
        NSAdjoint(const NSAdjoint<T, P>& _e);
        virtual ~NSAdjoint();

        void SetAlpha(int _i, T _alpha);

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        virtual void iUpdateMacro();
        virtual void iCollision();
        virtual void iExternalForce();

        bool CheckConvergence(T _eps);
        void SetFt(T _rho, T _ux, T _uy);

        virtual T GetRho(int _i, int _t) const;
        virtual T GetU(int _d, int _i, int _t) const;
        virtual T GetQ(int _i) const;
        virtual T GetV(int _d, int _i) const;
        virtual T GetR(int _d, int _i) const;
        virtual T GetSensitivity(int _i) const;
        
        const int nt;
        int t = 0, tmax;

protected:
        const int np;           //  np : number of particle
        P<T>* f;
        T omega;
        T *rho, *u[P<T>::nd], *q, *v[P<T>::nd], *r[P<T>::nd], *alpha, *sensitivity;
    };


    template<class T, template<class>class P>
    NSAdjoint<T, P>::NSAdjoint(P<T>* _f, T _viscosity, int _nt) : np(_f->np), nt(_nt) {
        assert(0 < _f->np && 0 < _nt);
        
        this->t = 0;
        this->f = _f;
        this->tmax = 0;
        this->omega = 1.0/(3.0*_viscosity*this->f->dt/(this->f->dx*this->f->dx) + 0.5);
        
        this->rho = new T[this->np*this->nt];
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np*this->nt];
        }
        
        for (int i = 0; i < this->np*this->nt; i++) {
            this->rho[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = T();
            }
        }

        this->q = new T[this->np];
         for (int k = 0; k < P<T>::nd; k++) {
            this->v[k] = new T[this->np];
            this->r[k] = new T[this->np];
        }
        this->alpha = new T[this->np];
        this->sensitivity = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->q[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = T();
                this->r[k][i] = T();
            }
            this->alpha[i] = T();
            this->sensitivity[i] = T();
        }
    }


    template<class T, template<class>class P>
    NSAdjoint<T, P>::NSAdjoint(const NSAdjoint<T, P>& _e) : np(_e.np), nt(_e.nt) {
        this->t = _e.t;
        this->f = _e.f;
        this->tmax = _e.tmax;
        this->omega = _e.omega;
        
        this->rho = new T[this->np];
        for (int k = 0; k < P<T>::nd; k++) {
            this->u[k] = new T[this->np*this->nt];
        }

        for (int i = 0; i < this->np*this->nt; i++) {
            this->rho[i] = _e.rho[i];
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][i] = _e.u[k][i];
            }
        }

        this->q = new T[this->np];
        for (int k = 0; k < P<T>::nd; k++) {
            this->v[k] = new T[this->np];
            this->r[k] = new T[this->np];
        }
        this->alpha = new T[this->np];
        this->sensitivity = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->q[i] = _e.q[i];
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = _e.v[k][i];
                this->r[k][i] = _e.r[k][i];
            }
            this->alpha[i] = _e.alpha[i];
            this->sensitivity[i] = _e.sensitivity[i];
        }
    }


    template<class T, template<class>class P>
    NSAdjoint<T, P>::~NSAdjoint() {
        delete[] this->rho;
        delete[] this->q;
        for (int k = 0; k < P<T>::nd; k++) {
            delete[] this->u[k];
            delete[] this->v[k];
            delete[] this->r[k];
        }
        delete[] this->alpha;
        delete[] this->sensitivity;
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::SetAlpha(int _i, T _alpha) {
        assert(0 <= _i && _i < this->np);
        this->alpha[_i] = _alpha;
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::UpdateMacro() {
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;
            this->rho[ti] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][ti] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                this->rho[ti] += this->f->ft[j][i];
                for (int k = 0; k < P<T>::nd; k++) {
                    this->u[k][ti] += P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][ti] /= this->rho[ti];
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::Collision() {
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][ti];
                    uu += this->u[k][ti]*this->u[k][ti];
                }
                T fieq = P<T>::ei[j]*this->rho[ti]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*fieq;
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::ExternalForce() {
        for (int i = 0; i < this->np; i++) {
            T tmprho = T();
            T tmpu[P<T>::nd];
            for (int k = 0; k < P<T>::nd; k++) {
                tmpu[k] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                tmprho += this->f->ft[j][i];
                for (int k = 0; k < P<T>::nd; k++) {
                    tmpu[k] += P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }
            for (int k = 0; k < P<T>::nd; k++) {
                tmpu[k] /= tmprho;
            }

            for (int j = 0; j < P<T>::nc; j++) {
                for (int k = 0; k < P<T>::nd; k++) {
                    this->f->ft[j][i] += -3.0*this->f->dx*P<T>::ei[j]*P<T>::ci[j][k]*this->alpha[i]*tmpu[k];
                }
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::iUpdateMacro() {
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;

            this->q[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = T();
                this->r[k][i] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][ti];
                    uu += this->u[k][ti]*this->u[k][ti];
                }
                this->q[i] += this->f->ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                for (int k = 0; k < P<T>::nd; k++) {
                    this->v[k][i] += this->f->ft[j][i]*P<T>::ei[j]*(P<T>::ci[j][k] + 3.0*ciu*P<T>::ci[j][k] - this->u[k][ti]);
                    this->r[k][i] += this->f->ft[j][i]*P<T>::ei[j]*P<T>::ci[j][k];
                }
            }
            for (int k = 0; k < P<T>::nd; k++) {
                this->sensitivity[i] += 3.0*this->f->dx*r[k][i]*this->u[k][ti];
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::iCollision() {
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;
            for (int j = 0; j < P<T>::nc; j++) {
                T feq = this->q[i];
                for (int k = 0; k < P<T>::nd; k++) {
                    feq += 3.0*this->v[k][i]*(P<T>::ci[j][k] - this->u[k][ti]);
                }
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*feq;
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::iExternalForce() {
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;

            T tmpm[P<T>::nd];
            for (int k = 0; k < P<T>::nd; k++) {
                tmpm[k] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    tmpm[k] += P<T>::ei[j]*P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }

            for (int j = 0; j < P<T>::nc; j++) {
                for (int k = 0; k < P<T>::nd; k++) {
                    this->f->ft[j][i] -= 3.0*this->alpha[i]*tmpm[k]*(P<T>::ci[j][k] - this->u[k][ti])/this->rho[ti];
                }
            }
        }
    }


    template<class T, template<class>class P>
    bool NSAdjoint<T, P>::CheckConvergence(T _eps) {
        T unorm = T();
        T dunorm = T();
        for (int i = 0; i < this->np; i++) {
            int ti = this->np*this->t + i;
            int tm1i = this->np*(this->t - 1) + i;
            for (int k = 0; k < P<T>::nd; k++) {
                unorm += this->u[k][ti]*this->u[k][ti];
                dunorm += (this->u[k][ti] - this->u[k][tm1i])*(this->u[k][ti] - this->u[k][tm1i]);
            }
        }
        this->tmax = this->t;
        return sqrt(dunorm/unorm) < _eps ? true : false;
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::SetFt(T _rho, T _ux, T _uy) {
        for (int i = 0; i < this->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = P<T>::ci[j][0]*_ux + P<T>::ci[j][0]*_uy;
                T uu = _ux*_ux + _uy*_uy;
                this->f->ft[j][i] = P<T>::ei[j]*_rho*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
            }
        }
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetRho(int _i, int _t) const {
        assert(0 <= _i && _i < this->np && 0 <= _t && _t < this->nt);
        return this->rho[this->np*_t + _i];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetU(int _d, int _i, int _t) const {
        assert(0 <= _i && _i < this->np && 0 <= _t && _t < this->nt);
        return this->u[_d][this->np*_t + _i];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetQ(int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->q[_i];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetV(int _d, int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->v[_d][_i];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetR(int _d, int _i) const {
        assert(0 <= _i && _i < this->np);
        return this->r[_d][_i];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetSensitivity(int _i) const {
        assert(0 <= _i && _i < this->np);
        T sens = T();
        for (int k = 0; k < P<T>::nd; k++) {
            sens += 3.0*this->r[k][_i]*this->u[k][this->np*this->tmax + _i];
        }
        return sens;
    }
}