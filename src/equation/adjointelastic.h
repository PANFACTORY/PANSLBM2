//*****************************************************************************
//  Title       :   src/equation/adjointelastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#ifdef _USE_AVX_DEFINES
    #include "../equation_avx/adjointelastic_avx.h"
#endif

namespace PANSLBM2 {
    namespace AEL {
        //  Function of updating macroscopic values of AEL for 2D
        template<class T, template<class>class P>
        void Macro(T &_irho, T &_imx, T &_imy, T &_isxx, T &_isxy, T &_isyx, T &_isyy, const T *_f0, const T *_f, int _idx) {
            _irho = P<T>::ei[0]*_f0[_idx];
            _imx = T();
            _imy = T();
            _isxx = T();
            _isxy = T();
            _isyx = T();
            _isyy = T();
            for (int c = 1; c < P<T>::nc; ++c) {
                T f = _f[P<T>::IndexF(_idx, c)];
                _irho += P<T>::ei[c]*f;
                _imx += P<T>::ei[c]*P<T>::cx[c]*f;
                _imy += P<T>::ei[c]*P<T>::cy[c]*f;
                _isxx += P<T>::ei[c]*P<T>::cx[c]*P<T>::cx[c]*f;
                _isxy += P<T>::ei[c]*P<T>::cx[c]*P<T>::cy[c]*f;
                _isyx += P<T>::ei[c]*P<T>::cy[c]*P<T>::cx[c]*f;
                _isyy += P<T>::ei[c]*P<T>::cy[c]*P<T>::cy[c]*f;
            }
        }

        //  Function of getting equilibrium of AEL for 2D
        template<class T, template<class>class P>
        void Equilibrium(T *_feq, T _irho, T _imx, T _imy, T _isxx, T _isxy, T _isyx, T _isyy, T _gamma) {
            for (int c = 0; c < P<T>::nc; ++c) {
                T imc = _imx*P<T>::cx[c] + _imy*P<T>::cy[c];
                T cisc = P<T>::cx[c]*_isxx*P<T>::cx[c] + P<T>::cx[c]*_isxy*P<T>::cy[c] + P<T>::cy[c]*_isyx*P<T>::cx[c] + P<T>::cy[c]*_isyy*P<T>::cy[c];
                T irhocc = _irho*(P<T>::cx[c]*P<T>::cx[c] + P<T>::cy[c]*P<T>::cy[c]);
                _feq[c] = 3.0*imc + 4.5*_gamma*cisc - 1.5*_gamma*irhocc;
            }
        }

        //  Function of setting boundary condition of AEL set iStress for D2Q9 along x edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetStressAlongXEdge(P<T>& _p, int _i, int _directionx, Fv0 _txbc, Fv1 _tybc, const T *_rho, Ff _bctype) {
            int i = _i - _p.offsetx;
            if (0 <= i && i < _p.nx) {
                for (int j = 0; j < _p.ny; ++j) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T tx = _txbc(i + _p.offsetx, j + _p.offsety), ty = _tybc(i + _p.offsetx, j + _p.offsety);
                        if (_directionx == -1) {
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] - (4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + 2.0*tx/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 5)] - (4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + 2.0*(tx - ty)/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 8)] - (4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 + 2.0*(tx + ty)/_rho[idx];
                        } else if (_directionx == 1) {
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] - (4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 - 2.0*tx/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 6)] - (4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 - 2.0*(tx + ty)/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 7)] - (4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0 - 2.0*(tx - ty)/_rho[idx];
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition of AEL set iStress for D2Q9 along y edge
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetStressAlongYEdge(P<T>& _p, int _j, int _directiony, Fv0 _txbc, Fv1 _tybc, const T *_rho, Ff _bctype) {
            int j = _j - _p.offsety;
            if (0 <= j && j < _p.ny) {
                for (int i = 0; i < _p.nx; ++i) {
                    if (_bctype(i + _p.offsetx, j + _p.offsety)) {
                        int idx = _p.Index(i, j);
                        T tx = _txbc(i + _p.offsetx, j + _p.offsety), ty = _tybc(i + _p.offsetx, j + _p.offsety);
                        if (_directiony == -1) {
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] - (4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 + 2.0*ty/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 6)] - (4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 + 2.0*(ty + tx)/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 5)] - (4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0 + 2.0*(ty - tx)/_rho[idx];
                        } else if (_directiony == 1) {
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] - (4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - 2.0*ty/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 8)] - (4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - 2.0*(ty + tx)/_rho[idx];
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 7)] - (4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0 - 2.0*(ty - tx)/_rho[idx];
                        }
                    }
                }
            }
        }
    }

    namespace AEL {
        //  Function of Update macro and Collide of AEL for 2D
        template<class T, template<class>class P>
        void MacroCollide(P<T>& _p, T *_irho, T *_imx, T *_imy, T *_isxx, T *_isxy, T *_isyx, T *_isyy, T _tau, const T *_gamma, bool _issave = false) {
            T omega = 1.0/_tau, iomega = 1.0 - omega, feq[P<T>::nc];
            #pragma omp parallel for private(feq)
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T irho, imx, imy, isxx, isxy, isyx, isyy;
                Macro<T, P>(irho, imx, imy, isxx, isxy, isyx, isyy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _irho[idx] = irho;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _isxx[idx] = isxx;
                    _isxy[idx] = isxy;
                    _isyx[idx] = isyx;
                    _isyy[idx] = isyy;
                }

                //  Collide
                Equilibrium<T, P>(feq, irho, imx, imy, isxx, isxy, isyx, isyy, _gamma[idx]);
                _p.f0[idx] = iomega*_p.f0[idx] + omega*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomega*_p.f[idxf] + omega*feq[c];
                }
            }
        }

        //  Function of setting initial condition of AEL for 2D
        template<class T, template<class>class P>
        void InitialCondition(P<T>& _p, const T *_irho, const T *_imx, const T *_imy, const T *_isxx, const T *_isxy, const T *_isyx, const T *_isyy, const T *_gamma) {
            T feq[P<T>::nc];
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                Equilibrium<T, P>(feq, _irho[idx], _imx[idx], _imy[idx], _isxx[idx], _isxy[idx], _isyx[idx], _isyy[idx], _gamma[idx]);
                _p.f0[idx] = feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = feq[c];
                }
            }
        }

        //  Function of setting boundary condition of AEL set iStress for 2D
        template<class T, template<class>class P, class Fv0, class Fv1, class Ff>
        void iBoundaryConditionSetStress(P<T>& _p, Fv0 _txbc, Fv1 _tybc, const T *_rho, Ff _bctype) {
            iBoundaryConditionSetStressAlongXEdge(_p, 0, -1, _txbc, _tybc, _rho, _bctype);          //  On xmin
            iBoundaryConditionSetStressAlongXEdge(_p, _p.lx - 1, 1, _txbc, _tybc, _rho, _bctype);   //  On xmax
            iBoundaryConditionSetStressAlongYEdge(_p, 0, -1, _txbc, _tybc, _rho, _bctype);          //  On ymin
            iBoundaryConditionSetStressAlongYEdge(_p, _p.ly - 1, 1, _txbc, _tybc, _rho, _bctype);   //  On ymax
        } 
    }
}