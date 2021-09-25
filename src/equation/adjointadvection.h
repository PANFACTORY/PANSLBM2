//*****************************************************************************
//  Title       :   src/equation/adjointadvection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/03
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>
#include "adjointnavierstokes.h"

namespace {
    const int SetT = 1;
    const int SetQ = 2;
}

namespace PANSLBM2 {
    namespace AAD {
        //  Function of updating macroscopic values of AAD for 2D
        template<class T, template<class>class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, const T *_g0, const T *_g, int _idx) {
            _item = Q<T>::ei[0]*_g0[_idx];
            _iqx = T();
            _iqy = T();
            for (int c = 1; c <Q<T>::nc; ++c) {
                _item += Q<T>::ei[c]*_g[Q<T>::IndexF(_idx, c)];
                _iqx += Q<T>::ei[c]*Q<T>::cx[c]*_g[Q<T>::IndexF(_idx, c)];
                _iqy += Q<T>::ei[c]*Q<T>::cy[c]*_g[Q<T>::IndexF(_idx, c)];
            }
        }

        //  Function of updating macroscopic values of AAD for 3D
        template<class T, template<class>class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, T &_iqz, const T *_g0, const T *_g, int _idx) {
            _item = Q<T>::ei[0]*_g0[_idx];
            _iqx = T();
            _iqy = T();
            _iqz = T();
            for (int c = 1; c <Q<T>::nc; ++c) {
                _item += Q<T>::ei[c]*_g[Q<T>::IndexF(_idx, c)];
                _iqx += Q<T>::ei[c]*Q<T>::cx[c]*_g[Q<T>::IndexF(_idx, c)];
                _iqy += Q<T>::ei[c]*Q<T>::cy[c]*_g[Q<T>::IndexF(_idx, c)];
                _iqz += Q<T>::ei[c]*Q<T>::cz[c]*_g[Q<T>::IndexF(_idx, c)];
            }
        }

        //  Function of getting equilibrium of AAD for 2D
        template<class T, template<class>class Q>
        T Equilibrium(T _item, T _iqx, T _iqy, T _ux, T _uy, int _c) {
            return _item + 3.0*(_ux*_iqx + _uy*_iqy);
        }

        //  Function of getting equilibrium of AAD for 3D
        template<class T, template<class>class Q>
        T Equilibrium(T _item, T _iqx, T _iqy, T _iqz, T _ux, T _uy, T _uz, int _c) {
            return _item + 3.0*(_ux*_iqx + _uy*_iqy + _uz*_iqz);
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(const T *_rho, const T *_ux, const T *_uy, T _imx, T _imy, const T *_tem, T _iqx, T _iqy, T _omegag, T *_f0, T *_f, const T *_alpha, int _idx) {
            _f0[_idx] += 3.0*(
                -_ux[_idx]*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                -_uy[_idx]*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy)
            )/(_rho[_idx] + _alpha[_idx]);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*(
                    (P<T>::cx[c] - _ux[_idx])*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                    (P<T>::cy[c] - _uy[_idx])*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy)
                )/(_rho[_idx] + _alpha[_idx]);
            }
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 3D
        template<class T, template<class>class P>
        void ExternalForceBrinkman(const T *_rho, const T *_ux, const T *_uy, const T *_uz, T _imx, T _imy, T _imz, const T *_tem, T _iqx, T _iqy, T _iqz, T _omegag, T *_f0, T *_f, const T *_alpha, int _idx) {
            _f0[_idx] += 3.0*(
                -_ux[_idx]*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                -_uy[_idx]*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy) +
                -_uz[_idx]*(_tem[_idx]*_iqz*_omegag - _alpha[_idx]*_imz)
            )/(_rho[_idx] + _alpha[_idx]);
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*(
                    (P<T>::cx[c] - _ux[_idx])*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                    (P<T>::cy[c] - _uy[_idx])*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy) +
                    (P<T>::cz[c] - _uz[_idx])*(_tem[_idx]*_iqz*_omegag - _alpha[_idx]*_imz)
                )/(_rho[_idx] + _alpha[_idx]);
            }
        }

        //  Function of applying external force with heat exchange of AAD for 2D/3D
        template<class T, template<class>class Q>
        void ExternalForceHeatExchange(T _item, T *_g0, T *_g, const T *_beta, int _idx) {
            _g0[_idx] -= _beta[_idx]*(1.0 + _item)/(1.0 + _beta[_idx]);
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] -= _beta[_idx]*(1.0 + _item)/(1.0 + _beta[_idx]);
            }
        }

        //  Function of applying external force with natural convection of AAD for 2D
        template<class T, template<class>class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _gx, T _gy, T *_g0, T *_g, int _idx) {
            _g0[_idx] += 3.0*(_imx*_gx + _imy*_gy);
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += 3.0*(_imx*_gx + _imy*_gy);
            }
        }

        //  Function of applying external force with natural convection of AAD for 3D
        template<class T, template<class>class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _imz, T _gx, T _gy, T _gz, T *_g0, T *_g, int _idx) {
            _g0[_idx] += 3.0*(_imx*_gx + _imy*_gy + _imz*_gz);
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += 3.0*(_imx*_gx + _imy*_gy + _imz*_gz);
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange), Collide and Stream of AAD for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamHeatExchange(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(item, _q.f0, _q.f, _beta, idx);
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                }
            }
        }
    
        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity constant)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity constant)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity heterogenious)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity heterogenious)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamForceConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity constant)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T _diffusivity, T _gx, T _gy, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, _gx, _gy, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity constant)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, T _diffusivity, T _gx, T _gy, T _gz, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity heterogenious)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, T _gx, T _gy, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                //  Update macro
                T ip, iux, iuy, imx, imy;
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                T item, iqx, iqy;
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, _gx, _gy, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity heterogenious)
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideStreamNaturalConvection(
            P<T>& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q<T>& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, T _gx, T _gy, T _gz, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                //  Update macro
                T ip, iux, iuy, iuz, imx, imy, imz;
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                T item, iqx, iqy, iqz;
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  External force with Brinkman model
                ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f0, _p.f, _alpha, idx);
                ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f0, _p.f, idx);
                ExternalForceNaturalConvection<T, Q>(imx, imy, imz, _gx, _gy, _gz, _q.f0, _q.f, idx);
                Macro<T, Q>(item, iqx, iqy, iqz, _q.f0, _q.f, idx);

                //  Save macro if need
                if (_issave) {
                    _ip[idx] = ip;
                    _iux[idx] = iux;
                    _iuy[idx] = iuy;
                    _iuz[idx] = iuz;
                    _imx[idx] = imx;
                    _imy[idx] = imy;
                    _imz[idx] = imz;
                    _item[idx] = item;
                    _iqx[idx] = iqx;
                    _iqy[idx] = iqy;
                    _iqz[idx] = iqz;
                }

                //  Collide
                _p.f0[idx] = (1.0 - omegaf)*_p.f0[idx] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, 0);
                for (int c = 1; c < P<T>::nc; ++c) {
                    _p.f[P<T>::IndexF(idx, c)] = (1.0 - omegaf)*_p.f[P<T>::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                }
                _q.f0[idx] = (1.0 - omegag)*_q.f0[idx] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = (1.0 - omegag)*_q.f[Q<T>::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of setting initial condition of AAD for 2D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_ux, const T *_uy, const T *_item, const T *_iqx, const T *_iqy) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _q.f0[idx] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _ux[idx], _uy[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of setting initial condition of AAD for 3D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _q.f0[idx] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _iqz[idx], _ux[idx], _uy[idx], _uz[idx], 0);
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _iqz[idx], _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D2Q9
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetT(Q<T>& _q, const T *_ux, const T *_uy, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        T rho0 = -(4.0*(1.0 + 3.0*_ux[idx])*_q.f[Q<T>::IndexF(idx, 1)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 5)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 8)])/(6.0*(1.0 + 3.0*_ux[idx]));
                        _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {                    
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        T rho0 = -(4.0*(1.0 - 3.0*_ux[idx])*_q.f[Q<T>::IndexF(idx, 3)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 6)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 7)])/(6.0*(1.0 - 3.0*_ux[idx]));
                        _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        T rho0 = -(4.0*(1.0 + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 2)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 5)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 6)])/(6.0*(1.0 + 3.0*_uy[idx]));
                        _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        T rho0 = -(4.0*(1.0 - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 4)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 7)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q<T>::IndexF(idx, 8)])/(6.0*(1.0 - 3.0*_uy[idx]));
                        _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D3Q15
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetT(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/12.0
                                -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]))
                                -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {            
                        if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(_q.nx - 1, j, k);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]))
                                -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                        }
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, 0, k);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/12.0
                                -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]))
                                -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, _q.ny - 1, k);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                -_uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]))
                                -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                        }
                    }
                }
            }
            //  On zmin
            if (_q.PEz == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz)) {
                            int idx = _q.Index(i, j, 0);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/12.0
                                -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]))
                                -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]));
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                        }
                    }
                }
            }
            //  On zmax
            if (_q.PEz == _q.mz - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz)) {
                            int idx = _q.Index(i, j, _q.nz - 1);
                            T rho0 = -(8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/12.0
                                -_ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]))
                                -_uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]));
                            _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }   
        }
    
        //  Function of setting boundary condition set iQ of AAD for D2Q9
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQ(Q<T>& _q, const T *_ux, const T *_uy, Ff _bctype, T _eps = T()) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        T rho0 = (
                            (1.0 + 3.0*_ux[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])
                            + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])
                            - 12.0*_eps
                        )/(6.0*(1.0 - 3.0*_ux[idx]));
                        _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        T rho0 = (
                            (1.0 - 3.0*_ux[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])
                            + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])
                            - 12.0*_eps
                        )/(6.0*(1.0 + 3.0*_ux[idx]));
                        _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        T rho0 = (
                            (1.0 + 3.0*_uy[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])
                            + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])
                            - 12.0*_eps
                        )/(6.0*(1.0 - 3.0*_uy[idx]));
                        _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        T rho0 = (
                            (1.0 - 3.0*_uy[idx])*(4.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])
                            + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])
                            - 12.0*_eps
                        )/(6.0*(1.0 + 3.0*_uy[idx]));
                        _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                        _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set iQ of AAD for D3Q15
        template<class T, template<class>class Q, class Ff>
        void iBoundaryConditionSetQ(Q<T>& _q, const T *_ux, const T *_uy, const T *_uz, Ff _bctype, T _eps = T()) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            T rho0 = (
                                (1.0 + 3.0*_ux[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])
                                + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])
                                + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])
                                - 24.0*_eps
                            )/(12.0*(1.0 - 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 4)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(_q.nx - 1, j, k);
                            T rho0 = (
                                (1.0 - 3.0*_ux[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])
                                + 3.0*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])
                                + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])
                                - 24.0*_eps
                            )/(12.0*(1.0 + 3.0*_ux[idx]));
                            _q.f[Q<T>::IndexF(idx, 1)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                        }
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, 0, k);
                            T rho0 = (
                                (1.0 + 3.0*_uy[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])
                                + 3.0*_uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])
                                + 3.0*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])
                                - 24.0*_eps
                            )/(12.0*(1.0 - 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 5)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, _q.ny - 1, k);
                            T rho0 = (
                                (1.0 - 3.0*_uy[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])
                                + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])
                                + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])
                                - 24.0*_eps
                            )/(12.0*(1.0 + 3.0*_uy[idx]));
                            _q.f[Q<T>::IndexF(idx, 2)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                        }
                    }
                }
            }
            //  On zmin
            if (_q.PEz == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz)) {
                            int idx = _q.Index(i, j, 0);
                            T rho0 = (
                                (1.0 + 3.0*_uz[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])
                                + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])
                                + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])
                                - 24.0*_eps
                            )/(12.0*(1.0 - 3.0*_uz[idx]));
                            _q.f[Q<T>::IndexF(idx, 6)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 10)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 11)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 12)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 13)] = rho0;
                        }
                    }
                }
            }
            //  On zmax
            if (_q.PEz == _q.mz - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz)) {
                            int idx = _q.Index(i, j, _q.nz - 1);
                            T rho0 = (
                                (1.0 - 3.0*_uz[idx])*(8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])
                                + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])
                                + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])
                                - 24.0*_eps
                            )/(12.0*(1.0 + 3.0*_uz[idx]));
                            _q.f[Q<T>::IndexF(idx, 3)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 7)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 8)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 9)] = rho0;
                            _q.f[Q<T>::IndexF(idx, 14)] = rho0;
                        }
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D2Q9
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRho(P<T>& _p, Q<T>& _q, const T *_rho, const T *_ux, const T *_uy, const T *_tem, Ff _bctype, T _eps = T()) { 
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _p.offsetx, j + _p.offsety)) {
                        int idx = _q.Index(0, j);
                        T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                        T flux0 = T();
                        if (_bctype(0 + _p.offsetx, j + _p.offsety) == SetT) {
                            flux0 = _tem[idx]*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/(2.0*(1.0 + 3.0*_ux[idx])*_rho[idx]);
                        } else if (_bctype(0 + _p.offsetx, j + _p.offsety) == SetQ) {
                            flux0 = -_tem[idx]*(
                                (4.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/3.0
                                + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/2.0
                            )/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                        _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety)) { 
                        int idx = _q.Index(_q.nx - 1, j);    
                        T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 7)])/3.0;
                        T flux0 = T();
                        if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety) == SetT) {
                            flux0 = _tem[idx]*_uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_ux[idx])*_rho[idx]);
                        } else if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety) == SetQ) {
                            flux0 = -_tem[idx]*(
                                (4.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/3.0
                                + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/2.0
                            )/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                        _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _p.offsetx, 0 + _p.offsety)) { 
                        int idx = _q.Index(i, 0);
                        T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 6)])/3.0;
                        T flux0 = T();
                        if (_bctype(i + _p.offsetx, 0 + _p.offsety) == SetT) {
                            flux0 = _tem[idx]*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/(2.0*(1.0 + 3.0*_uy[idx])*_rho[idx]);
                        } else if (_bctype(i + _p.offsetx, 0 + _p.offsety) == SetQ) {
                            flux0 = -_tem[idx]*(
                                (4.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/3.0
                                + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/2.0
                            )/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                        _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety)) {  
                        int idx = _q.Index(i, _q.ny - 1);  
                        T rho0 = -(4.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)])/3.0;
                        T flux0 = T();
                        if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety) == SetT) {
                            flux0 = _tem[idx]*_ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_uy[idx])*_rho[idx]);
                        } else if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety) == SetQ) {
                            flux0 = -_tem[idx]*(
                                (4.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/3.0 
                                + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 7)])/2.0
                            )/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                        _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                        _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15
        template<class T, template<class>class P, template<class>class Q, class Ff>
        void iBoundaryConditionSetRho(P<T>& _p, Q<T>& _q, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, Ff _bctypeg, T _eps = T()) {
            //  On xmin
            if (_p.PEx == 0) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(0, j, k);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 1)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 12)])/6.0;
                            T flux0 = T();
                            if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                            } else if (_bctype(0 + _p.offsetx, j + _p.offsety, k + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/6.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 4)] = _p.f[P<T>::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
            //  On xmax
            if (_p.PEx == _p.mx - 1) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(_p.nx - 1, j, k);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 4)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 13)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            T flux0 = T();
                            if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                            } else if (_bctype((_p.nx - 1) + _p.offsetx, j + _p.offsety, k + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 1)] = _p.f[P<T>::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
            //  On ymin
            if (_p.PEy == 0) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, 0, k);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 2)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                            T flux0 = T();
                            if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                            } else if (_bctype(i + _p.offsetx, 0 + _p.offsety, k + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/6.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 5)] = _p.f[P<T>::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
            //  On ymax
            if (_p.PEy == _p.my - 1) {
                for (int k = 0; k < _p.nz; ++k) {
                    for (int i = 0; i < _p.nx; ++i) {
                        if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz)) {
                            int idx = _p.Index(i, _p.ny - 1, k);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 5)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            T flux0 = T();
                            if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                            } else if (_bctype(i + _p.offsetx, (_p.ny - 1) + _p.offsety, k + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                    + _uz[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 2)] = _p.f[P<T>::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
            //  On zmin
            if (_p.PEz == 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz)) {
                            int idx = _p.Index(i, j, 0);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 3)] + _p.f[P<T>::IndexF(idx, 7)] + _p.f[P<T>::IndexF(idx, 8)] + _p.f[P<T>::IndexF(idx, 9)] + _p.f[P<T>::IndexF(idx, 14)])/6.0;
                            T flux0 = T();
                            if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                            } else if (_bctype(i + _p.offsetx, j + _p.offsety, 0 + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/6.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uz[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 6)] = _p.f[P<T>::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 10)] = _p.f[P<T>::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 11)] = _p.f[P<T>::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 12)] = _p.f[P<T>::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 13)] = _p.f[P<T>::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
            //  On zmax
            if (_p.PEz == _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz)) {
                            int idx = _p.Index(i, j, _p.nz - 1);
                            T rho0 = -(8.0*_p.f[P<T>::IndexF(idx, 6)] + _p.f[P<T>::IndexF(idx, 10)] + _p.f[P<T>::IndexF(idx, 11)] + _p.f[P<T>::IndexF(idx, 12)] + _p.f[P<T>::IndexF(idx, 13)])/6.0;
                            T flux0 = T();
                            if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz) == SetT) {
                                flux0 = _tem[idx]*(
                                    _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                            } else if (_bctype(i + _p.offsetx, j + _p.offsety, (_p.nz - 1) + _p.offsetz) == SetQ) {
                                flux0 = -_tem[idx]*(
                                    (8.0*_q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/6.0
                                    + _ux[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                    + _uy[idx]*(_q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/4.0
                                )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                            }
                            T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uz[idx])*_rho[idx]);
                            _p.f[P<T>::IndexF(idx, 3)] = _p.f[P<T>::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 7)] = _p.f[P<T>::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 8)] = _p.f[P<T>::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 9)] = _p.f[P<T>::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                            _p.f[P<T>::IndexF(idx, 14)] = _p.f[P<T>::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                        }
                    }
                }
            }
        }
    
        //  Function of getting sensitivity of temperature at heat source for D2Q9
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            const T *_ux, const T *_uy, const T *_imx, const T *_imy,
            Q<T>& _q, const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_g0, const T *_g, const *_ig0, const T *_ig,
            T *_dfds, const T *_diffusivity, const T *_dads, const T *_dkds, Fv _qnbc, Ff _bctype
        ) {
            //  Brinkman term and diffusivity term
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                T sumg = _g0[idx]*_ig0[idx];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    sumg += _ig[Q<T>::IndexF(idx, c)]*_g[Q<T>::IndexF(idx, c)];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
            }

            //  Boundary term along xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        _dfds[idx] += _qnbc(0 + _q.offsetx, j + _q.offsety)*_dkds[idx]*(
                            (1.0 + 3.0*_ux[idx])*(-6.0 + 4.0*_ig[Q<T>::IndexF(idx, 1)] + _ig[Q<T>::IndexF(idx, 5)] + _ig[Q<T>::IndexF(idx, 8)])
                            + 3.0*_uy[idx]*(_ig[Q<T>::IndexF(idx, 5)] - _ig[Q<T>::IndexF(idx, 8)])
                        )/(36.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
            //  Boundary term along xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        _dfds[idx] += _qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety)*_dkds[idx]*(
                            (1.0 - 3.0*_ux[idx])*(-6.0 + 4.0*_ig[Q<T>::IndexF(idx, 3)] + _ig[Q<T>::IndexF(idx, 6)] + _ig[Q<T>::IndexF(idx, 7)])
                            + 3.0*_uy[idx]*(_ig[Q<T>::IndexF(idx, 6)] - _ig[Q<T>::IndexF(idx, 7)])
                        )/(36.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
            //  Boundary term along ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        _dfds[idx] += _qnbc(i + _q.offsetx, 0 + _q.offsety)*_dkds[idx]*(
                            (1.0 + 3.0*_uy[idx])*(-6.0 + 4.0*_ig[Q<T>::IndexF(idx, 2)] + _ig[Q<T>::IndexF(idx, 5)] + _ig[Q<T>::IndexF(idx, 6)])
                            + 3.0*_ux[idx]*(_ig[Q<T>::IndexF(idx, 5)] - _ig[Q<T>::IndexF(idx, 6)])
                        )/(36.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
            //  Boundary term along ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        _dfds[idx] += _qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety)*_dkds[idx]*(
                            (1.0 - 3.0*_uy[idx])*(-6.0 + 4.0*_ig[Q<T>::IndexF(idx, 4)] + _ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 8)])
                            + 3.0*_ux[idx]*(_ig[Q<T>::IndexF(idx, 8)] - _ig[Q<T>::IndexF(idx, 7)])
                        )/(36.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15
        template<class T, template<class>class Q, class Fv, class Ff>
        void SensitivityTemperatureAtHeatSource(
            const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz,
            Q<T>& _q, const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz, const T *_g0, const T *_g, const T *_ig0, const T *_ig,
            T *_dfds, const T *_diffusivity, const T *_dads, const T *_dkds, Fv _qnbc, Ff _bctype
        ) {
            //  Brinkman term and diffusivity term
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                T sumg = _g0[idx]*_ig0[idx];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    sumg += _ig[Q<T>::IndexF(idx, c)]*_g[Q<T>::IndexF(idx, c)];
                }
                _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
            }

            //  Boundary term along xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            _dfds[idx] += _qnbc(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)*_dkds[idx]*(
                                (1.0 + 3.0*_ux[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 1)] + _ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 9)] + _ig[Q<T>::IndexF(idx, 10)] + _ig[Q<T>::IndexF(idx, 12)])
                                + 3.0*_uy[idx]*(_ig[Q<T>::IndexF(idx, 7)] - _ig[Q<T>::IndexF(idx, 9)] + _ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 12)])
                                + 3.0*_uz[idx]*(_ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 9)] - _ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 12)])
                            )/(72.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
            //  Boundary term along xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(_q.nx - 1, j, k);
                            _dfds[idx] += _qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz)*_dkds[idx]*(
                                (1.0 - 3.0*_ux[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 4)] + _ig[Q<T>::IndexF(idx, 8)] + _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 13)] + _ig[Q<T>::IndexF(idx, 14)])
                                + 3.0*_uy[idx]*(_ig[Q<T>::IndexF(idx, 8)] - _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 13)] - _ig[Q<T>::IndexF(idx, 14)])
                                + 3.0*_uz[idx]*(_ig[Q<T>::IndexF(idx, 8)] - _ig[Q<T>::IndexF(idx, 11)] - _ig[Q<T>::IndexF(idx, 13)] + _ig[Q<T>::IndexF(idx, 14)])
                            )/(72.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
            //  Boundary term along ymin
            if (_q.PEy == 0) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, 0, k);
                            _dfds[idx] += _qnbc(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz)*_dkds[idx]*(
                                (1.0 + 3.0*_uy[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 2)] + _ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 8)] + _ig[Q<T>::IndexF(idx, 10)] + _ig[Q<T>::IndexF(idx, 13)])
                                + 3.0*_uz[idx]*(_ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 8)] - _ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 13)])
                                + 3.0*_ux[idx]*(_ig[Q<T>::IndexF(idx, 7)] - _ig[Q<T>::IndexF(idx, 8)] + _ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 13)])
                            )/(72.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
            //  Boundary term along ymax
            if (_q.PEy == _q.my - 1) {
                for (int k = 0; k < _q.nz; ++k) {
                    for (int i = 0; i < _q.nx; ++i) {
                        if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(i, _q.ny - 1, k);
                            _dfds[idx] += _qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz)*_dkds[idx]*(
                                (1.0 - 3.0*_uy[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 5)] + _ig[Q<T>::IndexF(idx, 9)] + _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 12)] + _ig[Q<T>::IndexF(idx, 14)])
                                + _uz[idx]*(_ig[Q<T>::IndexF(idx, 9)] - _ig[Q<T>::IndexF(idx, 11)] - _ig[Q<T>::IndexF(idx, 12)] + _ig[Q<T>::IndexF(idx, 14)])
                                + _ux[idx]*(_ig[Q<T>::IndexF(idx, 9)] - _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 12)] - _ig[Q<T>::IndexF(idx, 14)])
                            )/(72.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
            //  Boundary term along zmin
            if (_q.PEz == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz)) {
                            int idx = _q.Index(i, j, 0);
                            _dfds[idx] += _qnbc(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz)*_dkds[idx]*(
                                (1.0 + 3.0*_uz[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 3)] + _ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 8)] + _ig[Q<T>::IndexF(idx, 9)] + _ig[Q<T>::IndexF(idx, 14)])
                                + _ux[idx]*(_ig[Q<T>::IndexF(idx, 7)] - _ig[Q<T>::IndexF(idx, 8)] + _ig[Q<T>::IndexF(idx, 9)] - _ig[Q<T>::IndexF(idx, 14)])
                                + _uy[idx]*(_ig[Q<T>::IndexF(idx, 7)] + _ig[Q<T>::IndexF(idx, 8)] - _ig[Q<T>::IndexF(idx, 9)] - _ig[Q<T>::IndexF(idx, 14)])
                            )/(72.0*(1.0 - 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
            //  Boundary term along zmax
            if (_q.PEz == _q.mz - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    for (int j = 0; j < _q.ny; ++j) {
                        if (_bctype(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz)) {
                            int idx = _q.Index(i, j, _q.nz - 1);
                            _dfds[idx] += _qnbc(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz)*_dkds[idx]*(
                                (1.0 - 3.0*_uz[idx])*(-12.0 + 8.0*_ig[Q<T>::IndexF(idx, 6)] + _ig[Q<T>::IndexF(idx, 10)] + _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 12)] + _ig[Q<T>::IndexF(idx, 13)])
                                + _ux[idx]*(_ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 11)] + _ig[Q<T>::IndexF(idx, 12)] - _ig[Q<T>::IndexF(idx, 13)])
                                + _uy[idx]*(_ig[Q<T>::IndexF(idx, 10)] - _ig[Q<T>::IndexF(idx, 11)] - _ig[Q<T>::IndexF(idx, 12)] + _ig[Q<T>::IndexF(idx, 13)])
                            )/(72.0*(1.0 + 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                        }
                    }
                }
            }
        }
    }
}