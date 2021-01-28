//*****************************************************************************
//  Title       :   src/equation/elastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/01/27
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace EL {
        //*********************************************************************
        //  Elastic 2D  :   Update macroscopic values
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _p, T* _rho, T* _ux, T* _uy, T* _sxx, T* _sxy, T* _syx, T* _syy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _p.np; i++) {
                _ux[i] = T();
                _uy[i] = T();
                _sxx[i] = T();
                _sxy[i] = T();
                _syx[i] = T();
                _syy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _ux[i] += P<T>::cx[j]*_p.ft[j][i];
                    _uy[i] += P<T>::cy[j]*_p.ft[j][i];
                    _sxx[i] -= P<T>::cx[j]*P<T>::cx[j]*_p.ft[j][i];
                    _sxy[i] -= P<T>::cx[j]*P<T>::cy[j]*_p.ft[j][i];
                    _syx[i] -= P<T>::cy[j]*P<T>::cx[j]*_p.ft[j][i];
                    _syy[i] -= P<T>::cy[j]*P<T>::cy[j]*_p.ft[j][i]; 
                }
                _ux[i] /= _rho[i];
                _uy[i] /= _rho[i];
            }
        }

        //*********************************************************************
        //  Elastic 2D  :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _elasticy, P<T>& _p, T* _rho, T* _ux, T* _uy, T* _sxx, T* _sxy, T* _syx, T* _syy) {
            assert(P<T>::nd == 2);
//  要修正：緩和時間はelasticyの関数になるはずだがその具体的形が不明のため
            T omega = 1.0/0.8;
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T cu = P<T>::cx[j]*_ux[i] + P<T>::cy[j]*_uy[i];
                    T csc = P<T>::cx[j]*_sxx[i]*P<T>::cx[j] + P<T>::cx[j]*_sxy[i]*P<T>::cy[j] + P<T>::cy[j]*_syx[i]*P<T>::cx[j] + P<T>::cy[j]*_syy[i]*P<T>::cy[j];
                    T trs = _sxx[i] + _syy[i];
                    T feq = P<T>::ei[j]*(3.0*_rho[i]*cu - 4.5*csc + 1.5*trs);
                    _p.ftp1[j][i] = (1.0 - omega)*_p.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Elastic 2D  :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _p, T _rho, T _ux, T _uy, T _sxx, T _sxy, T _syx, T _syy) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _p.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T cu = P<T>::cx[j]*_ux + P<T>::cy[j]*_uy;
                T csc = P<T>::cx[j]*_sxx*P<T>::cx[j] + P<T>::cx[j]*_sxy*P<T>::cy[j] + P<T>::cy[j]*_syx*P<T>::cx[j] + P<T>::cy[j]*_syy*P<T>::cy[j];
                T trs = _sxx + _syy;
                _p.ft[j][_i] = P<T>::ei[j]*(3.0*_rho*cu - 4.5*csc + 1.5*trs);
            }
        }
    }
}