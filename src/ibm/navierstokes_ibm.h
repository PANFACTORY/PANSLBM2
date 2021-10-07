#pragma once
#include "src/equation/navierstokes.h"

namespace PANSLBM2{
    template<class T, template<class>class P>
    class IBNS {
public:
        IBNS() = delete;
        IBNS(const P<T>& _p, int _nb, T _dv, int _lmax = 5, T _eps = T(1e-5)) : nx(_p.nx), ny(_p.ny), nxyz(_p.nxyz), nb(_nb) {
            this->uxl = new T[this->nxyz];  this->uyl = new T[this->nxyz];
            this->gxl = new T[this->nxyz];  this->gyl = new T[this->nxyz];

            this->bpx = new T[this->nb];    this->bpy = new T[this->nb];
            this->bvx = new T[this->nb];    this->bvy = new T[this->nb];
            this->bux = new T[this->nb];    this->buy = new T[this->nb];
            this->bgx = new T[this->nb];    this->bgy = new T[this->nb];
            
            this->dv = _dv;
            this->eps = _eps;
            this->lmax = _lmax;
            this->width = 2;
            /*
(-width+1=)-1 +0 +1 +2(= width)
            |  | *|  |
                    ↑
            (int)bpx
            */
        }
        IBNS(const IBNS<T, P>&) = delete;
        ~IBNS() {
            delete[] this->uxl, this->uyl, this->gxl, this->gyl;
            delete[] this->bpx, this->bpy, this->bvx, this->bvy, this->bux, this->buy, this->bgx, this->bgy;
        }

        void Update(P<T>& _p, T *_rho, T *_ux, T *_uy);

        void SetBP(int _n, T _bpx, T _bpy) {    this->bpx[_n] = _bpx;   this->bpy[_n] = _bpy;   }
        void SetBV(int _n, T _bvx, T _bvy) {    this->bvx[_n] = _bvx;   this->bvy[_n] = _bvy;   }
        void GetBU(int _n, T& _bux, T& _buy) {  _bux = this->bux[_n];   _buy = this->buy[_n];   }
        void GetBG(int _n, T& _bgx, T& _bgy) {  _bgx = -this->bgx[_n];  _bgy = -this->bgy[_n];   }

private:
        const int nx, ny, nxyz, nb;
        T dv, eps;
        int lmax, width;
        T *uxl, *uyl, *gxl, *gyl;
        T *bpx, *bpy, *bvx, *bvy, *bux, *buy, *bgx, *bgy;

        static T W(T _r) {
            T absr = fabs(_r);
            if (absr < 1.0) {
                return (3.0 - 2.0*absr + sqrt(1.0 + 4.0*absr - 4.0*pow(absr, 2.0)))/8.0;
            } else if (absr < 2.0) {
                return (5.0 - 2.0*absr - sqrt(-7.0 + 12.0*absr - 4.0*pow(absr, 2.0)))/8.0;
            } else {
                return T();
            }
        }
    };

    template<class T, template<class>class P>
    void IBNS<T, P>::Update(P<T>& _p, T *_rho, T *_ux, T *_uy) {
        //  Update macro
        for (int idx = 0; idx < _p.nxyz; ++idx) {
            NS::Macro<T, P>(_rho[idx], _ux[idx], _uy[idx], _p.f0, _p.f, idx);
        }

        //**********STEP0**********
        for (int n = 0; n < this->nb; ++n) {
            this->bux[n] = T();
            this->buy[n] = T();

            for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                int i = (int)this->bpx[n] + ii;
                if (0 <= i && i < this->nx) {
                    for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                        int j = (int)this->bpy[n] + jj;
                        if (0 <= j && j < this->ny) {
                            this->bux[n] += _ux[_p.Index(i, j)]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n]);
                            this->buy[n] += _uy[_p.Index(i, j)]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n]);
                        }
                    }
                }
            }
            
            this->bgx[n] = this->bvx[n] - this->bux[n];
            this->bgy[n] = this->bvy[n] - this->buy[n];
        }

        for (int l = 0; l < this->lmax; l++) {
            //**********STEP1**********
            for (int idx = 0; idx < this->nxyz; ++idx) {
                this->gxl[idx] = T();
                this->gyl[idx] = T();
            } 

            for (int n = 0; n < this->nb; ++n) {
                for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                    int i = (int)this->bpx[n] + ii;
                    if (0 <= i && i < this->nx) {
                        for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                            int j = (int)this->bpy[n] + jj;
                            if (0 <= j && j < this->ny) {
                                this->gxl[_p.Index(i, j)] += this->bgx[n]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n])*this->dv;
                                this->gyl[_p.Index(i, j)] += this->bgy[n]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n])*this->dv;
                            }
                        }
                    }
                }
            }

            //**********STEP2**********
            for (int idx = 0; idx < this->nxyz; ++idx) {
                this->uxl[idx] = _ux[idx] + this->gxl[idx];
                this->uyl[idx] = _uy[idx] + this->gyl[idx];
            } 

            //**********STEP3**********
            T unorm = T();
            for (int n = 0; n < this->nb; ++n) {
                this->bux[n] = T();
                this->buy[n] = T();

                for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                    int i = (int)this->bpx[n] + ii;
                    if (0 <= i && i < this->nx) {
                        for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                            int j = (int)this->bpy[n] + jj;
                            if (0 <= j && j < this->ny) {
                                this->bux[n] += this->uxl[_p.Index(i, j)]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n]);
                                this->buy[n] += this->uyl[_p.Index(i, j)]*IBNS<T, P>::W(i - this->bpx[n])*IBNS<T, P>::W(j - this->bpy[n]);
                            }
                        }
                    }
                }

                T tmpnorm = sqrt(pow(this->bvx[n] - this->bux[n], 2.0) + pow(this->bvy[n] - this->buy[n], 2.0));
                unorm = unorm < tmpnorm ? tmpnorm : unorm;
            }
            
            //**********STEP4**********
            if (unorm < this->eps) {
                break;
            }

            for (int n = 0; n < this->nb; ++n) {
                this->bgx[n] += this->bvx[n] - this->bux[n];
                this->bgy[n] += this->bvy[n] - this->buy[n]; 
            }
        }
    
        //**********STEP5**********
        for (int idx = 0; idx < _p.nxyz; ++idx) {
            for (int c = 1; c < P<T>::nc; ++c) {
                _p.f[P<T>::IndexF(idx, c)] += 3.0*P<T>::ei[c]*(P<T>::cx[c]*this->gxl[idx] + P<T>::cy[c]*this->gyl[idx]);
            } 
        }
    }
}