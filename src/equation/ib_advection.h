//*****************************************************************************
//  Title       :   src/equation/ib_advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/02/17
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include <cmath>

namespace PANSLBM2 {
    namespace AD {
        //*********************************************************************
        //  Advection 2D    :   External Force with Immersed boundary method (Fix thermal)
        //*********************************************************************
        template<class T, template<class>class P>
        void ExternalForceIBThermal(
            P<T>& _p, T *_tem, T *_teml, T *_gl, 
            int _nb, T *_bpx, T *_bpy, T *_btem, T *_bg, T *_btemd,
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
                _btem[n] = T();
            
                for (int i = -1; i <= 2; i++) {
                    int ii = (int)_bpx[n] + i;
                    if (0 <= ii && ii < _p.nx) {
                        for (int j = -1; j <= 2; j++) {
                            int jj = (int)_bpy[n] + j;
                            if (0 <= jj && jj < _p.ny) {
                                _btem[n] += _tem[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                            }
                        }
                    }
                }
                
                _bg[n] = (_btemd[n] - _btem[n])/_p.dx;
            }

            for (int l = 0; l < _lmax; l++) {
                //**********STEP1**********
                for (int ij = 0; ij < _p.np; ij++) {
                    _gl[ij] = T();
                } 

                for (int n = 0; n < _nb; n++) {
                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _gl[_p.GetIndex(ii, jj)] += _bg[n]*W(ii - _bpx[n])*W(jj - _bpy[n])*_dv;
                                }
                            }
                        }
                    }
                }

                //**********STEP2**********
                for (int ij = 0; ij < _p.np; ij++) {
                    _teml[ij] = _tem[ij] + _gl[ij]*_p.dx;
                } 

                //**********STEP3**********
                T unorm = T();
                for (int n = 0; n < _nb; n++) {
                    _btem[n] = T();

                    for (int i = -1; i <= 2; i++) {
                        int ii = (int)_bpx[n] + i;
                        if (0 <= ii && ii < _p.nx) {
                            for (int j = -1; j <= 2; j++) {
                                int jj = (int)_bpy[n] + j;
                                if (0 <= jj && jj < _p.ny) {
                                    _btem[n] += _teml[_p.GetIndex(ii, jj)]*W(ii - _bpx[n])*W(jj - _bpy[n]);
                                }
                            }
                        }
                    }

                    unorm = unorm < fabs(_btemd[n] - _btem[n]) ? fabs(_btemd[n] - _btem[n]) : unorm;
                }
                
                //**********STEP4**********
                if (unorm < _eps) {
                    break;
                }

                for (int n = 0; n < _nb; n++) {
                    _bg[n] += (_btemd[n] - _btem[n])*_p.dx;
                }
            }

            //**********STEP5**********
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    _p.ft[j][i] += _p.dx*P<T>::ei[j]*_gl[i];
                }
            } 
        }
    }
}