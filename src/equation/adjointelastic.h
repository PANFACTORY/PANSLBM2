//*****************************************************************************
//  Title       :   src/equation/adjointelastic.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/01/28
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>

namespace PANSLBM2 {
    namespace EL {
        //*********************************************************************
        //  Elastic 2D  :   Expand macroscopic values
        //*********************************************************************
        template<class T, template<class>class P>
        void ExpandMacro(P<T>& _p, T* _sxx, T* _sxy, T* _syx, T* _syy, T* _gamma) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _p.np; i++) {
                _sxx[i] *= _gamma[i];
                _sxy[i] *= _gamma[i];
                _syx[i] *= _gamma[i];
                _syy[i] *= _gamma[i];
            }
        }
    }

    namespace AEL {
        //*********************************************************************
        //  Adjoint elastic 2D  :   Update macroscopic values
        //*********************************************************************
        template<class T, template<class>class P>
        void UpdateMacro(P<T>& _p, T* _irho, T* _imx, T* _imy, T* _isxx, T* _isxy, T* _isyx, T* _isyy) {
            assert(P<T>::nd == 2);
            for (int i = 0; i < _p.np; i++) {
                _irho[i] = T();
                _imx[i] = T();
                _imy[i] = T();
                _isxx[i] = T();
                _isxy[i] = T();
                _isyx[i] = T();
                _isyy[i] = T();
                for (int j = 0; j < P<T>::nc; j++) {
                    _irho[i] += P<T>::ei[j]*_p.ft[j][i];
                    _imx[i] += P<T>::ei[j]*P<T>::cx[j]*_p.ft[j][i];
                    _imy[i] += P<T>::ei[j]*P<T>::cy[j]*_p.ft[j][i];
                    _isxx[i] += P<T>::ei[j]*P<T>::cx[j]*P<T>::cx[j]*_p.ft[j][i];
                    _isxy[i] += P<T>::ei[j]*P<T>::cx[j]*P<T>::cy[j]*_p.ft[j][i];
                    _isyx[i] += P<T>::ei[j]*P<T>::cy[j]*P<T>::cx[j]*_p.ft[j][i];
                    _isyy[i] += P<T>::ei[j]*P<T>::cy[j]*P<T>::cy[j]*_p.ft[j][i];
                }
            }
        }

        //*********************************************************************
        //  Adjoint elastic 2D  :   Collision term
        //*********************************************************************
        template<class T, template<class>class P>
        void Collision(T _elasticy, P<T>& _p, T* _irho, T* _imx, T* _imy, T* _isxx, T* _isxy, T* _isyx, T* _isyy, T* _gamma) {
            assert(P<T>::nd == 2);
//  要修正：緩和時間はelasticyの関数になるはずだがその具体的形が不明のため
            T omega = 1.0/(3.0*_elasticy*_p.dt/(_p.dx*_p.dx) + 0.5);
            for (int i = 0; i < _p.np; i++) {
                for (int j = 0; j < P<T>::nc; j++) {
                    T imc = _imx[i]*P<T>::cx[j] + _imy[i]*P<T>::cy[j];
                    T cisc = P<T>::cx[j]*_isxx[i]*P<T>::cx[j] + P<T>::cx[j]*_isxy[i]*P<T>::cy[j] + P<T>::cy[j]*_isyx[i]*P<T>::cx[j] + P<T>::cy[j]*_isyy[i]*P<T>::cy[j];
                    T irhocc = _irho[i]*(P<T>::cx[j]*P<T>::cx[j] + P<T>::cy[j]*P<T>::cy[j]);
                    T feq = 3.0*imc + 4.5*_gamma[i]*cisc - 1.5*_gamma[i]*irhocc;
                    _p.ftp1[j][i] = (1.0 - omega)*_p.ft[j][i] + omega*feq;
                }
            }
        }

        //*********************************************************************
        //  Adjoint elastic 2D  :   Initial condition
        //*********************************************************************
        template<class T, template<class>class P>
        void InitialCondition(int _i, P<T>& _p, T _irho, T _imx, T _imy, T _isxx, T _isxy, T _isyx, T _isyy, T _gamma) {
            assert(P<T>::nd == 2 && 0 <= _i && _i < _p.np);
            for (int j = 0; j < P<T>::nc; j++) {
                T imc = _imx*P<T>::cx[j] + _imy*P<T>::cy[j];
                T cisc = P<T>::cx[j]*_isxx*P<T>::cx[j] + P<T>::cx[j]*_isxy*P<T>::cy[j] + P<T>::cy[j]*_isyx*P<T>::cx[j] + P<T>::cy[j]*_isyy*P<T>::cy[j];
                T irhocc = _irho*(P<T>::cx[j]*P<T>::cx[j] + P<T>::cy[j]*P<T>::cy[j]);
                _p.ft[j][_i] = 3.0*imc + 4.5*_gamma*cisc -1.5*_gamma*irhocc;
            }
        }
    }
}