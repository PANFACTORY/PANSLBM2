#pragma once
#include "src/equation/advection.h"

namespace PANSLBM2{
    template<class T, template<class>class Q>
    class IBAD {
public:
        IBAD() = delete;
        IBAD(const Q<T>& _q, int _nb, T _dv, int _lmax = 5, T _eps = T(1e-5)) : nx(_q.nx), ny(_q.ny), nxyz(_q.nxyz), nb(_nb) {
            this->teml = new T[this->nxyz];
            this->gl = new T[this->nxyz];

            this->bpx = new T[this->nb];    
            this->bpy = new T[this->nb];
            this->btem = new T[this->nb];
            this->bg = new T[this->nb];
            this->btemd = new T[this->nb];

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
        IBAD(const IBAD<T, Q>&) = delete;
        ~IBAD() {
            delete[] this->teml, this->gl;
            delete[] this->bpx, this->bpy, this->btem, this->bg, this->btemd;
        }

        void Update(const Q<T>& _q, const T *_tem);
        void ExternalForceIB(T *_g0, T *_g, int _idx) const;

        void SetBP(int _n, T _bpx, T _bpy) {    this->bpx[_n] = _bpx;   this->bpy[_n] = _bpy;   }
        void SetBT(int _n, T _btem) {  this->btem[_n] = _btem;   }
        void GetBG(int _n, T& _bg) {  _bg = -this->bg[_n];   }

private:
        const int nx, ny, nxyz, nb;
        T dv, eps;
        int lmax, width;
        T *teml, *gl;
        T *bpx, *bpy, *btem, *bg, *btemd;

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

    template<class T, template<class>class Q>
    void IBAD<T, Q>::Update(const Q<T>& _q, const T *_tem) {
        //**********STEP0**********
        for (int n = 0; n < this->nb; ++n) {
            this->btem[n] = T();

            for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                int i = (int)this->bpx[n] + ii;
                if (0 <= i && i < this->nx) {
                    for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                        int j = (int)this->bpy[n] + jj;
                        if (0 <= j && j < this->ny) {
                            this->btem[n] += _tem[_q.Index(i, j)]*IBAD<T, Q>::W(i - this->bpx[n])*IBAD<T, Q>::W(j - this->bpy[n]);
                        }
                    }
                }
            }
            
            this->bg[n] = this->btemd[n] - this->btem[n];
        }

        for (int l = 0; l < this->lmax; l++) {
            //**********STEP1**********
            for (int idx = 0; idx < this->nxyz; ++idx) {
                this->gl[idx] = T();
            } 

            for (int n = 0; n < this->nb; ++n) {
                for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                    int i = (int)this->bpx[n] + ii;
                    if (0 <= i && i < this->nx) {
                        for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                            int j = (int)this->bpy[n] + jj;
                            if (0 <= j && j < this->ny) {
                                this->gl[_q.Index(i, j)] += this->bg[n]*IBAD<T, Q>::W(i - this->bpx[n])*IBAD<T, Q>::W(j - this->bpy[n])*this->dv;
                            }
                        }
                    }
                }
            }

            //**********STEP2**********
            for (int idx = 0; idx < this->nxyz; ++idx) {
                this->teml[idx] = _tem[idx] + this->gl[idx];
            } 

            //**********STEP3**********
            T unorm = T();
            for (int n = 0; n < this->nb; ++n) {
                this->btem[n] = T();
                
                for (int ii = -this->width + 1; ii <= this->width; ++ii) {
                    int i = (int)this->bpx[n] + ii;
                    if (0 <= i && i < this->nx) {
                        for (int jj = -this->width + 1; jj <= this->width; ++jj) {
                            int j = (int)this->bpy[n] + jj;
                            if (0 <= j && j < this->ny) {
                                this->btem[n] += this->teml[_q.Index(i, j)]*IBAD<T, Q>::W(i - this->bpx[n])*IBAD<T, Q>::W(j - this->bpy[n]);
                            }
                        }
                    }
                }

                T tmpnorm = fabs(this->btemd[n] - this->btem[n]);
                unorm = unorm < tmpnorm ? tmpnorm : unorm;
            }
            
            //**********STEP4**********
            if (unorm < this->eps) {
                break;
            }

            for (int n = 0; n < this->nb; ++n) {
                this->bg[n] += this->btemd[n] - this->btem[n];
            }
        }
    }

    template<class T, template<class>class Q>
    void IBAD<T, Q>::ExternalForceIB(T *_g0, T *_g, int _idx) const {
        //**********STEP5**********
        _g0[_idx] += Q<T>::ei[0]*this->gl[_idx];
        for (int c = 1; c < Q<T>::nc; ++c) {
            _g[Q<T>::IndexF(_idx, c)] += Q<T>::ei[c]*this->gl[_idx];
        } 
    }

    namespace AD {
        //  Function of Update macro and Collide of force convection for 2D
        template<class T, template<class>class P, template<class>class Q, template<class, template<class>class>class B>
        void MacroCollideForceConvectionIB(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            const B<T, Q>& _b, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                _b.ExternalForceIB(_q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of natural convection for 2D
        template<class T, template<class>class P, template<class>class Q, template<class, template<class>class>class B>
        void MacroCollideNaturalConvectionIB(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            T _gx, T _gy, T _tem0, const B<T, Q>& _b, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc]; 
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
            #pragma omp parallel for private(feq, geq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                _b.ExternalForceIB(_q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  Save macro if need
                if (_issave) {
                    _rho[idx] = rho;
                    _ux[idx] = ux;
                    _uy[idx] = uy;
                    _tem[idx] = tem;
                    _qx[idx] = qx;
                    _qy[idx] = qy;
                }

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    }
}