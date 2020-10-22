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

        template<class F>
        void SetAlpha(F _f);

        virtual void UpdateMacro();
        virtual void Collision();
        virtual void ExternalForce();

        virtual void iUpdateMacro();
        virtual void iCollision();
        virtual void iExternalForce();

        bool CheckConvergence(T _eps);
        void SwitchDirection();

        virtual T GetRho(int _i, int _j, int _t) const;
        virtual T GetU(int _d, int _i, int _j, int _t) const;
        virtual T GetQ(int _i, int _j) const;
        virtual T GetV(int _d, int _i, int _j) const;
        virtual T GetSensitivity(int _i, int _j) const;
        
        const int nt;
        int t = 0;

protected:
        const int np;           //  np : number of particle
        P<T>* f;
        int tmax;
        T omega;
        T *rho, *u[P<T>::nd], *q, *v[P<T>::nd], *alpha, *sensitivity;
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
        }
        this->alpha = new T[this->np];
        this->sensitivity = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->q[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = T();
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
        }
        this->alpha = new T[this->np];
        this->sensitivity = new T[this->np];

        for (int i = 0; i < this->np; i++) {
            this->q[i] = _e.q[i];
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = _e.v[k][i];
            }
            this->alpha[i] = _e.alpha[i];
            this->sensitivity[i] = _e.sensitivity[i];
        }
    }


    template<class T, template<class>class P>
    NSAdjoint<T, P>::~NSAdjoint() {
        delete[] this->rho;
        for (int i = 0; i < P<T>::nd; i++) {
            delete[] this->u[i];
        }
        delete[] this->q;
        for (int i = 0; i < P<T>::nd; i++) {
            delete[] this->v[i];
        }
        delete[] this->alpha;
        delete[] this->sensitivity;
    }


    template<class T, template<class>class P>
    template<class F>
    void NSAdjoint<T, P>::SetAlpha(F _f) {
        for (int i = 0; i < this->f->nx; i++) {
            for (int j = 0; j < this->f->ny; j++) {
                this->alpha[this->f->ny*i + j] = _f(i, j);
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::UpdateMacro() {
        for (int i = 0; i < this->np; i++) {
            int ii = this->np*this->t + i;
            this->rho[ii] = T();
            for (int j = 0; j < P<T>::nd; j++) {
                this->u[j][ii] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                this->rho[ii] += this->f->ft[j][i];
                for (int k = 0; k < P<T>::nd; k++) {
                    this->u[k][ii] += P<T>::ci[j][k]*this->f->ft[j][i];
                }
            }
            for (int k = 0; k < P<T>::nd; k++) {
                this->u[k][ii] /= this->rho[ii];
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::Collision() {
        for (int i = 0; i < this->np; i++) {
            int ii = this->np*this->t + i;
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += P<T>::ci[j][k]*this->u[k][ii];
                    uu += this->u[k][ii]*this->u[k][ii];
                }
                T fieq = P<T>::ei[j]*this->rho[ii]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
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
            int ii = this->np*this->t + i;

            this->q[i] = T();
            for (int k = 0; k < P<T>::nd; k++) {
                this->v[k][i] = T();
            }
            for (int j = 0; j < P<T>::nc; j++) {
                T ciu = T();
                T uu = T();
                for (int k = 0; k < P<T>::nd; k++) {
                    ciu += -P<T>::ci[j][k]*this->u[k][ii];
                    uu += this->u[k][ii]*this->u[k][ii];
                }
                this->q[i] += this->f->ft[j][i]*P<T>::ei[j]*(1.0 + 3.0*ciu + 4.5*ciu*ciu - 1.5*uu);
                for (int k = 0; k < P<T>::nd; k++) {
                    this->v[k][i] += this->f->ft[j][i]*P<T>::ei[j]*(-P<T>::ci[j][k] - 3.0*ciu*P<T>::ci[j][k] - this->u[k][ii]);
                }
                this->sensitivity[i] += 3.0*this->f->dx*this->f->ft[j][i]*P<T>::ei[j]*ciu;
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::iCollision() {
        for (int i = 0; i < this->np; i++) {
            int ii = this->np*this->t + i;
            for (int j = 0; j < P<T>::nc; j++) {
                T feq = this->q[i];
                for (int k = 0; k < P<T>::nd; k++) {
                    feq += 3.0*this->v[k][i]*(-P<T>::ci[j][k] - this->u[k][ii]);
                }
                this->f->ftp1[j][i] = (1.0 - this->omega)*this->f->ft[j][i] + this->omega*feq;
            }
        }
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::iExternalForce() {
        for (int i = 0; i < this->np; i++) {
            T tmpm[P<T>::nd];
            for (int k = 0; k < P<T>::nd; k++) {
                tmpm[k] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    tmpm[k] += -P<T>::ei[j]*P<T>::ci[j][k]*this->f->ft[j][i];
                }
                tmpm[k] *= 3.0*this->f->dx*this->alpha[i];
            }

            int ii = this->np*this->t + i;
            for (int j = 0; j < P<T>::nc; j++) {
                for (int k = 0; k < P<T>::nd; k++) {
                    this->f->ft[j][i] -= tmpm[k]*(-P<T>::ci[j][k] - this->u[k][ii])/this->rho[ii];
                }
            }
        }
    }


    template<class T, template<class>class P>
    bool NSAdjoint<T, P>::CheckConvergence(T _eps) {
        T unorm = T();
        T dunorm = T();
        for (int i = 0; i < this->np; i++) {
            int ii = this->np*this->t + i;
            int jj = this->np*(this->t - 1) + i;
            for (int k = 0; k < P<T>::nd; k++) {
                unorm += this->u[k][ii]*this->u[k][ii];
                dunorm += (this->u[k][ii] - this->u[k][jj])*(this->u[k][ii] - this->u[k][jj]);
            }
        }
        this->tmax = this->t;
        return sqrt(dunorm/unorm) < _eps ? true : false;
    }


    template<class T, template<class>class P>
    void NSAdjoint<T, P>::SwitchDirection() {
        for (int i = 0; i < this->np; i++) {
            for (int j = 0; j < P<T>::nc; j++) {
                this->f->ft[j][i] = T();
                this->f->ftp1[j][i] = T();
            }
        }
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetRho(int _i, int _j, int _t) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->rho[this->np*_t + this->f->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetU(int _d, int _i, int _j, int _t) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->u[_d][this->np*_t + this->f->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetQ(int _i, int _j) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->q[this->f->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetV(int _d, int _i, int _j) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->v[_d][this->f->ny*_i + _j];
    }


    template<class T, template<class>class P>
    T NSAdjoint<T, P>::GetSensitivity(int _i, int _j) const {
        assert(0 <= _i && _i < this->f->nx && 0 <= _j && _j < this->f->ny);
        return this->sensitivity[this->f->ny*_i + _j];
    }
}