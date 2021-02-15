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
        void ExternalForceIB(P<T>& _p, T *_ux, T *_uy, T *_uxl, T *_uyl, T *_gxl, T *_gyl, int _nb, T *_bpx, T *_bpy, T *_bux, T *_buy, T *_bgx, T *_bgy, T _dv, int _lmax = 5, T _eps = T(1.0e-5)) {
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
            for (int k = 0; k < _nb; k++) {
                _bux[k] = T();
                _buy[k] = T();

                /*
                -1 +0 +1 +2
                |  | *|  |
                   ↑
                (int)_bpx
                */
                for (int i = -1; i <= 2; i++) {
                    int ii = (int)_bpx[k] + i;
                    if (0 <= ii && ii < _p.nx) {
                        for (int j = -1; j <= 2; j++) {
                            int jj = (int)_bpy[k] + j;
                            if (0 <= jj && jj < _p.ny) {
                                _bux[k] += _ux[_p.GetIndex(ii, jj)]*W(ii - _bpx[k])*W(jj - _bpy[k]);
                                _buy[k] += _uy[_p.GetIndex(ii, jj)]*W(ii - _bpx[k])*W(jj - _bpy[k]);
                            }
                        }
                    }
                }
                
                _bgx[k] = -_bux[k]/_p.dx;
                _bgy[k] = -_buy[k]/_p.dx;
            }

            for (int l = 0; l < _lmax; l++) {
                //**********STEP1**********
                for (int ij = 0; ij < _p.np; ij++) {
                    _gxl[ij] = T();
                    _gyl[ij] = T();
                } 

                for (int k = 0; k < _nb; k++) {
                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[k] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[k] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _gxl[_p.GetIndex(ii, jj)] += _bgx[k]*W(ii - _bpx[k])*W(jj - _bpy[k])*_dv;
                                    _gyl[_p.GetIndex(ii, jj)] += _bgy[k]*W(ii - _bpx[k])*W(jj - _bpy[k])*_dv;
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
                for (int k = 0; k < _nb; k++) {
                    _bux[k] = T();
                    _buy[k] = T();

                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[k] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[k] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _bux[k] += _uxl[_p.GetIndex(ii, jj)]*W(ii - _bpx[k])*W(jj - _bpy[k]);
                                    _buy[k] += _uyl[_p.GetIndex(ii, jj)]*W(ii - _bpx[k])*W(jj - _bpy[k]);
                                }
                            }
                        }
                    }

                    unorm = unorm < sqrt(pow(_bux[k], 2.0) + pow(_buy[k], 2.0)) ? sqrt(pow(_bux[k], 2.0) + pow(_buy[k], 2.0)) : unorm;
                }
                
                //**********STEP4**********
                if (unorm < _eps) {
                    break;
                }

                for (int k = 0; k < _nb; k++) {
                    _bgx[k] += -_bux[k]*_p.dx;
                    _bgy[k] += -_buy[k]*_p.dx; 
                }
            }

            //**********STEP5**********
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _p.ft[j][i] += 3.0*_p.dx*P<T>::ei[j]*(P<T>::cx[j]*_gxl[i] + P<T>::cy[j]*_gyl[i]);
                }
            } 
        }
    }
}