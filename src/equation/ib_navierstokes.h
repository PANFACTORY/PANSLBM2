//*****************************************************************************
//  Title       :   src/equation/ib_navierstokes.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/15
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <cmath>

namespace PANSLBM2 {
    namespace NS {
        //*********************************************************************
        //  Navier-Stokes 2D    :   External Force with Immersed boundary method
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceIB(
            P<T>& _p, T *_ux, T *_uy, T *_uxl, T *_uyl, T *_gxl, T *_gyl, 
            int _nb, T *_bpx, T *_bpy, T *_bux, T *_buy, T *_bgx, T *_bgy, T *_bvx, T *_bvy,
            T _dv, int _lmax = 5, T _eps = T(1.0e-5)
        ) {
            assert(P<T>::nd == 2);

            auto W = [](T _r) {
                T absr = fabs(_r);
                if (absr < 1.0) {
                    return (3.0 - 2.0*absr + sqrt(1.0 + 4.0*absr - 4.0*pow(absr, 2.0)))/8.0;
                } else if (absr < 2.0) {
                    return (5.0 - 2.0*absr - sqrt(-7.0 + 12.0*absr - 4.0*pow(absr, 2.0)))/8.0;
                } else {
                    return T();
                }
            };

            //**********STEP0**********
            for (int n = 0; n < _nb; n++) {
                _bux[n] = T();
                _buy[n] = T();

                /*
                -1 +0 +1 +2
                |  | *|  |
                   ↑
                (int)_bpx
                */
                for (int i = -1; i <= 2; i++) {
                    int ii = (int)_bpx[n] + i;
                    if (0 <= ii && ii < _p.nx) {
                        for (int j = -1; j <= 2; j++) {
                            int jj = (int)_bpy[n] + j;
                            if (0 <= jj && jj < _p.ny) {
                                _bux[n] += _ux[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                                _buy[n] += _uy[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                            }
                        }
                    }
                }
                
                _bgx[n] = (_bvx[n] - _bux[n])/_p.dx;
                _bgy[n] = (_bvy[n] - _buy[n])/_p.dx;
            }

            for (int l = 0; l < _lmax; l++) {
                //**********STEP1**********
                for (int ij = 0; ij < _p.np; ij++) {
                    _gxl[ij] = T();
                    _gyl[ij] = T();
                } 

                for (int n = 0; n < _nb; n++) {
                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _gxl[_p.GetIndex(ii, jj)] += _bgx[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*_dv;
                                    _gyl[_p.GetIndex(ii, jj)] += _bgy[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*_dv;
                                }
                            }
                        }
                    }
                }

                //**********STEP2**********
                for (int ij = 0; ij < _p.np; ij++) {
                    _uxl[ij] = _ux[ij] + _gxl[ij]*_p.dx;
                    _uyl[ij] = _uy[ij] + _gyl[ij]*_p.dx;
                } 

                //**********STEP3**********
                T unorm = T();
                for (int n = 0; n < _nb; n++) {
                    _bux[n] = T();
                    _buy[n] = T();

                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _bux[n] += _uxl[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                                    _buy[n] += _uyl[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                                }
                            }
                        }
                    }

                    unorm = unorm < sqrt(pow(_bvx[n] - _bux[n], 2.0) + pow(_bvy[n] - _buy[n], 2.0)) ? sqrt(pow(_bvx[n] - _bux[n], 2.0) + pow(_bvy[n] - _buy[n], 2.0)) : unorm;
                }
                
                //**********STEP4**********
                if (unorm < _eps) {
                    break;
                }

                for (int n = 0; n < _nb; n++) {
                    _bgx[n] += (_bvx[n] - _bux[n])*_p.dx;
                    _bgy[n] += (_bvy[n] - _buy[n])*_p.dx; 
                }
            }

            //**********STEP5**********
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _p.ft[j][i] += 3.0*_p.dx*P<T>::ei[j]*(P<T>::cx[j]*_gxl[i] + P<T>::cy[j]*_gyl[i]);
                }
            } 
        }


        //*********************************************************************
        //  Navier-Stokes 3D    :   External Force with Immersed boundary method
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceIB(
            P<T>& _p, T *_ux, T *_uy, T *_uz, T *_uxl, T *_uyl, T *_uzl, T *_gxl, T *_gyl, T *_gzl, 
            int _nb, T *_bpx, T *_bpy, T *_bpz, T *_bux, T *_buy, T *_buz, T *_bgx, T *_bgy, T *_bgz, T *_bvx, T *_bvy, T *_bvz,
            T _dv, int _lmax = 5, T _eps = T(1.0e-5)
        ) {
            assert(P<T>::nd == 3);

            auto W = [](T _r) {
                T absr = fabs(_r);
                if (absr < 1.0) {
                    return (3.0 - 2.0*absr + sqrt(1.0 + 4.0*absr - 4.0*pow(absr, 2.0)))/8.0;
                } else if (absr < 2.0) {
                    return (5.0 - 2.0*absr - sqrt(-7.0 + 12.0*absr - 4.0*pow(absr, 2.0)))/8.0;
                } else {
                    return T();
                }
            };

            //**********STEP0**********
            for (int n = 0; n < _nb; n++) {
                _bux[n] = T();
                _buy[n] = T();
                _buz[n] = T();

                /*
                -1 +0 +1 +2
                |  | *|  |
                   ↑
                (int)_bpx
                */
                for (int i = -1; i <= 2; i++) {
                    int ii = (int)_bpx[n] + i;
                    if (0 <= ii && ii < _p.nx) {
                        for (int j = -1; j <= 2; j++) {
                            int jj = (int)_bpy[n] + j;
                            if (0 <= jj && jj < _p.ny) {
                                for (int k = 0; k <= 2; k++) {
                                    int kk = (int)_bpz[n] + k;
                                    if (0 <= kk && kk < _p.nz) {
                                        _bux[n] += _ux[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                        _buy[n] += _uy[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                        _buz[n] += _uz[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                    }
                                }
                            }
                        }
                    }
                }
                
                _bgx[n] = (_bvx[n] - _bux[n])/_p.dx;
                _bgy[n] = (_bvy[n] - _buy[n])/_p.dx;
                _bgz[n] = (_bvz[n] - _buz[n])/_p.dx;
            }

            for (int l = 0; l < _lmax; l++) {
                //**********STEP1**********
                for (int i = 0; i < _p.np; i++) {
                    _gxl[i] = T();
                    _gyl[i] = T();
                    _gzl[i] = T();
                } 

                for (int n = 0; n < _nb; n++) {
                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    for (int k = -1; k <= 2; k++) {
                                        int kk = (int)_bpz[n] + k;
                                        if (0 <= kk && kk < _p.nz) {
                                            _gxl[_p.GetIndex(ii, jj, kk)] += _bgx[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n])*_dv;
                                            _gyl[_p.GetIndex(ii, jj, kk)] += _bgy[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n])*_dv;
                                            _gzl[_p.GetIndex(ii, jj, kk)] += _bgz[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n])*_dv;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                //**********STEP2**********
                for (int i = 0; i < _p.np; i++) {
                    _uxl[i] = _ux[i] + _gxl[i]*_p.dx;
                    _uyl[i] = _uy[i] + _gyl[i]*_p.dx;
                    _uzl[i] = _uz[i] + _gzl[i]*_p.dx;
                } 

                //**********STEP3**********
                T unorm = T();
                for (int n = 0; n < _nb; n++) {
                    _bux[n] = T();
                    _buy[n] = T();
                    _buz[n] = T();

                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    for (int k = -1; k <= 2; k++) {
                                        int kk = (int)_bpz[n] + k;
                                        if (0 <= kk && kk < _p.nz) {
                                            _bux[n] += _uxl[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                            _buy[n] += _uyl[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                            _buz[n] += _uzl[_p.GetIndex(ii, jj, kk)]*W(ii - _bpx[n])*W(jj - _bpy[n])*W(kk - _bpz[n]);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    unorm = unorm < sqrt(pow(_bvx[n] - _bux[n], 2.0) + pow(_bvy[n] - _buy[n], 2.0) + pow(_bvz[n] - _buz[n], 2.0)) ? sqrt(pow(_bvx[n] - _bux[n], 2.0) + pow(_bvy[n] - _buy[n], 2.0) + pow(_bvz[n] - _buz[n], 2.0)) : unorm;
                }
                
                //**********STEP4**********
                if (unorm < _eps) {
                    break;
                }

                for (int n = 0; n < _nb; n++) {
                    _bgx[n] += (_bvx[n] - _bux[n])*_p.dx;
                    _bgy[n] += (_bvy[n] - _buy[n])*_p.dx; 
                    _bgz[n] += (_bvz[n] - _buz[n])*_p.dx; 
                }
            }

            //**********STEP5**********
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _p.ft[j][i] += 3.0*_p.dx*P<T>::ei[j]*(P<T>::cx[j]*_gxl[i] + P<T>::cy[j]*_gyl[i] + P<T>::cz[j]*_gzl[i]);
                }
            } 
        }
    }
}