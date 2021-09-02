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

namespace PANSLBM2 {
    namespace {
        const int SetiT = 1;
        const int SetiQ = 2;
    }

    namespace AAD {
        //  Function of updating macroscopic values of AAD for 2D
        template<class T, class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, const T *_g, int _idx) {
            _item = T();
            _iqx = T();
            _iqy = T();
            for (int c = 0; c <Q::nc; ++c) {
                _item += Q::ei[c]*_g[Q::IndexF(_idx, c)];
                _iqx += Q::ei[c]*Q::cx[c]*_g[Q::IndexF(_idx, c)];
                _iqy += Q::ei[c]*Q::cy[c]*_g[Q::IndexF(_idx, c)];
            }
        }

        //  Function of updating macroscopic values of AAD for 3D
        template<class T, class Q>
        void Macro(T &_item, T &_iqx, T &_iqy, T &_iqz, const T *_g, int _idx) {
            _item = T();
            _iqx = T();
            _iqy = T();
            _iqz = T();
            for (int c = 0; c <Q::nc; ++c) {
                _item += Q::ei[c]*_g[Q::IndexF(_idx, c)];
                _iqx += Q::ei[c]*Q::cx[c]*_g[Q::IndexF(_idx, c)];
                _iqy += Q::ei[c]*Q::cy[c]*_g[Q::IndexF(_idx, c)];
                _iqz += Q::ei[c]*Q::cz[c]*_g[Q::IndexF(_idx, c)];
            }
        }

        //  Function of getting equilibrium of AAD for 2D
        template<class T, class Q>
        T Equilibrium(T _item, T _iqx, T _iqy, T _ux, T _uy, int _c) {
            return _item + 3.0*(_ux*_iqx + _uy*_iqy);
        }

        //  Function of getting equilibrium of AAD for 3D
        template<class T, class Q>
        T Equilibrium(T _item, T _iqx, T _iqy, T _iqz, T _ux, T _uy, T _uz, int _c) {
            return _item + 3.0*(_ux*_iqx + _uy*_iqy + _uz*_iqz);
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 2D
        template<class T, class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, T _imx, T _imy, const T *_tem, T _iqx, T _iqy, T _omegag, T *_f, const T *_alpha, int _idx
        ) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] += 3.0*(
                    (P::cx[c] - _ux[_idx])*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                    (P::cy[c] - _uy[_idx])*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy)
                )/(_rho[_idx] + _alpha[_idx]);
            }
        }

        //  Function of applying external force with Brinkman model and advection of AAD for 3D
        template<class T, class P>
        void ExternalForceBrinkman(
            const T *_rho, const T *_ux, const T *_uy, const T *_uz, T _imx, T _imy, T _imz, const T *_tem, T _iqx, T _iqy, T _iqz, T _omegag, T *_f, const T *_alpha, int _idx
        ) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] += 3.0*(
                    (P::cx[c] - _ux[_idx])*(_tem[_idx]*_iqx*_omegag - _alpha[_idx]*_imx) + 
                    (P::cy[c] - _uy[_idx])*(_tem[_idx]*_iqy*_omegag - _alpha[_idx]*_imy) +
                    (P::cz[c] - _uz[_idx])*(_tem[_idx]*_iqz*_omegag - _alpha[_idx]*_imz)
                )/(_rho[_idx] + _alpha[_idx]);
            }
        }

        //  Function of applying external force with heat exchange of AAD for 2D/3D
        template<class T, class Q>
        void ExternalForceHeatExchange(T _item, T *_g, const T *_beta, int _idx) {
            for (int c = 0; c < Q::nc; ++c) {
                _g[Q::IndexF(_idx, c)] -= _beta[_idx]*(1.0 + _item)/(1.0 + _beta[_idx]);
            }
        }

        //  Function of applying external force with natural convection of AAD for 2D
        template<class T, class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _gx, T _gy, T *_g, int _idx) {
            for (int c = 0; c < Q::nc; ++c) {
                _g[Q::IndexF(_idx, c)] += 3.0*(_imx*_gx + _imy*_gy);
            }
        }

        //  Function of applying external force with natural convection of AAD for 3D
        template<class T, class Q>
        void ExternalForceNaturalConvection(T _imx, T _imy, T _imz, T _gx, T _gy, T gz, T *_g, int _idx) {
            for (int c = 0; c < Q::nc; ++c) {
                _g[Q::IndexF(_idx, c)] += 3.0*(_imx*_gx + _imy*_gy + _imz*_gz);
            }
        }

        //  Function of Update macro, External force(Brinkman, Heat exchange), Collide and Stream of AAD for 2D
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_HeatExchange(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    ExternalForceHeatExchange<T, Q>(item, _q.f, _beta, idx);
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

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

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }
    
        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);

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

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T ip, iux, iuy, iuz, imx, imy, imz;
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        T item, iqx, iqy, iqz;
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f, _alpha, idx);
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);

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

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i - P::cx[c], j - P::cy[c], k - P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c], k - Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);

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

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                        //  Update macro
                        T ip, iux, iuy, iuz, imx, imy, imz;
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        T item, iqx, iqy, iqz;
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f, _alpha, idx);
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);

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

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i - P::cx[c], j - P::cy[c], k - P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c], k - Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T _diffusivity, T _gx, T _gy, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    ExternalForceNaturalConvection(imx, imy, _gx, _gy, _q.f, idx);
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

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

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, T _diffusivity, T _gx, T _gy, T _gz, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T ip, iux, iuy, iuz, imx, imy, imz;
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        T item, iqx, iqy, iqz;
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f, _alpha, idx);
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        ExternalForceNaturalConvection(imx, imy, imz, _gx, _gy, _gz, _q.f, idx);
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

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

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i - P::cx[c], j - P::cy[c], k - P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c], k - Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 2D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, T *_ip, T *_iux, T *_iuy, T *_imx, T *_imy, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, const T *_diffusivity, T _gx, T _gy, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                    //  Update macro
                    T ip, iux, iuy, imx, imy;
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    T item, iqx, iqy;
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

                    //  External force with Brinkman model
                    ExternalForceBrinkman<T, P>(_rho, _ux, _uy, imx, imy, _tem, iqx, iqy, omegag, _p.f, _alpha, idx);
                    ANS::Macro<T, P>(ip, iux, iuy, imx, imy, _rho, _ux, _uy, _p.f, idx);
                    ExternalForceNaturalConvection(imx, imy, _gx, _gy, _q.f, idx);
                    Macro<T, Q>(item, iqx, iqy, _q.f, idx);

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

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i - P::cx[c], j - P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], ip, iux, iuy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of AAD for 3D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, const T *_rho, const T *_ux, const T *_uy, const T *_uz, T *_ip, T *_iux, T *_iuy, T *_iuz, T *_imx, T *_imy, T *_imz, const T *_alpha, T _viscosity,
            Q& _q, const T *_tem, T *_item, T *_iqx, T *_iqy, T *_iqz, const T *_diffusivity, T _gx, T _gy, T _gz, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                        //  Update macro
                        T ip, iux, iuy, iuz, imx, imy, imz;
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        T item, iqx, iqy, iqz;
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

                        //  External force with Brinkman model
                        ExternalForceBrinkman<T, P>(_rho, _ux, _uy, _uz, imx, imy, imz, _tem, iqx, iqy, iqz, omegag, _p.f, _alpha, idx);
                        ANS::Macro<T, P>(ip, iux, iuy, iuz, imx, imy, imz, _rho, _ux, _uy, _uz, _p.f, idx);
                        ExternalForceNaturalConvection(imx, imy, imz, _gx, _gy, _gz, _q.f, idx);
                        Macro<T, Q>(item, iqx, iqy, iqz, _q.f, idx);

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

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i - P::cx[c], j - P::cy[c], k - P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*ANS::Equilibrium<T, P>(_ux[idx], _uy[idx], _uz[idx], ip, iux, iuy, iuz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i - Q::cx[c], j - Q::cy[c], k - Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(item, iqx, iqy, iqz, _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of setting initial condition of AAD for 2D
        template<class T, class Q>
        void InitialCondition(Q& _q, const T *_ux, const T *_uy, const T *_item, const T *_iqx, const T *_iqy) {
            for (int idx = 0; idx < _q.nxy; ++idx) {
                for (int c = 0; c < Q::nc; ++c) {
                    _q.f[Q::IndexF(idx, c)] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _ux[idx], _uy[idx], c);
                }
            }
        }

        //  Function of setting initial condition of AAD for 3D
        template<class T, class Q>
        void InitialCondition(Q& _q, const T *_ux, const T *_uy, const T *_uz, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz) {
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                for (int c = 0; c < Q::nc; ++c) {
                    _q.f[Q::IndexF(idx, c)] = Equilibrium<T, Q>(_item[idx], _iqx[idx], _iqy[idx], _iqz[idx], _ux[idx], _uy[idx], _uz[idx], c);
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D2Q9
        template<class T, class Q>
        void BoundaryConditionSetiT(Q& _q, const T *_ux, const T *_uy, const int *_bctype) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == SetiT) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T rho0 = -(4.0*(1.0 + 3.0*_ux[idx])*_q.f[Q::IndexF(idx, 1)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 5)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 8)])/(6.0*(1.0 + 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 3)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == SetiT) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T rho0 = -(4.0*(1.0 - 3.0*_ux[idx])*_q.f[Q::IndexF(idx, 3)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 6)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 7)])/(6.0*(1.0 - 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 1)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == SetiT) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T rho0 = -(4.0*(1.0 + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 2)] + (1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 5)] + (1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 6)])/(6.0*(1.0 + 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 4)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == SetiT) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T rho0 = -(4.0*(1.0 - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 4)] + (1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 7)] + (1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])*_q.f[Q::IndexF(idx, 8)])/(6.0*(1.0 - 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 2)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                }
            }
        }

        //  Function of setting boundary condition set iT of AAD for D3Q15
        template<class T, class Q>
        void BoundaryConditionSetiT(Q& _q, const T *_ux, const T *_uy, const T *_uz, const int *_bctype) {
            int idx, idxbc;

            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  On xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 12)])/12.0
                            -_uy[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]))
                            -_uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/(4.0*(1.0 + 3.0*_ux[idx]));
                        _q.f[Q::IndexF(idx, 4)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }

                    //  On xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/12.0
                            -_uy[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] - _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]))
                            -_uz[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_ux[idx]));
                        _q.f[Q::IndexF(idx, 1)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  On ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 13)])/12.0
                            -_uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]))
                            -_ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/(4.0*(1.0 + 3.0*_uy[idx]));
                        _q.f[Q::IndexF(idx, 5)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }

                    //  On ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/12.0
                            -_uz[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]))
                            -_ux[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 - 3.0*_uy[idx]));
                        _q.f[Q::IndexF(idx, 2)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  On zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 14)])/12.0
                            -_ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]))
                            -_uy[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/(4.0*(1.0 + 3.0*_uz[idx]));
                        _q.f[Q::IndexF(idx, 6)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                    }

                    //  On zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixT) {
                        T rho0 = -(8.0*_q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/12.0
                            -_ux[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]))
                            -_uy[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/(4.0*(1.0 - 3.0*_uz[idx]));
                        _q.f[Q::IndexF(idx, 3)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set iQ of AAD for D2Q9
        template<class T, class Q>
        void BoundaryConditionSetiQ(Q& _q, const T *_ux, const T *_uy, const int *_bctype, T _eps = T()) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == SetiQ) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T rho0 = (
                        (1.0 + 3.0*_ux[idx])*(4.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])
                        + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])
                        - 12.0*_eps
                    )/(6.0*(1.0 - 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 3)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == SetiQ) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T rho0 = (
                        (1.0 - 3.0*_ux[idx])*(4.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])
                        + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])
                        - 12.0*_eps
                    )/(6.0*(1.0 + 3.0*_ux[idx]));
                    _q.f[Q::IndexF(idx, 1)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == SetiQ) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T rho0 = (
                        (1.0 + 3.0*_uy[idx])*(4.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])
                        + 3.0*_ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])
                        - 12.0*_eps
                    )/(6.0*(1.0 - 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 4)] = rho0;
                    _q.f[Q::IndexF(idx, 7)] = rho0;
                    _q.f[Q::IndexF(idx, 8)] = rho0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == SetiQ) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T rho0 = (
                        (1.0 - 3.0*_uy[idx])*(4.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])
                        + 3.0*_ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])
                        - 12.0*_eps
                    )/(6.0*(1.0 + 3.0*_uy[idx]));
                    _q.f[Q::IndexF(idx, 2)] = rho0;
                    _q.f[Q::IndexF(idx, 5)] = rho0;
                    _q.f[Q::IndexF(idx, 6)] = rho0;
                }
            }
        }

        //  Function of setting boundary condition set iQ of AAD for D3Q15
        template<class T, class Q>
        void BoundaryConditionSetiQ(Q& _q, const T *_ux, const T *_uy, const T *_uz, const int *_bctype, T _eps = T()) {
            int idx, idxbc;

            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  On xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 + 3.0*_ux[idx])*(8.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 12)])
                            + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])
                            + 3.0*_uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])
                            - 24.0*_eps
                        )/(12.0*(1.0 - 3.0*_ux[idx]));
                        _q.f[Q::IndexF(idx, 4)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }

                    //  On xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 - 3.0*_ux[idx])*(8.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])
                            + 3.0*_uy[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] - _q.f[Q::IndexF(idx, 14)])
                            + 3.0*_uz[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])
                            - 24.0*_eps
                        )/(12.0*(1.0 + 3.0*_ux[idx]));
                        _q.f[Q::IndexF(idx, 1)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  On ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 + 3.0*_uy[idx])*(8.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 13)])
                            + 3.0*_uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])
                            + 3.0*_ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])
                            - 24.0*_eps
                        )/(12.0*(1.0 - 3.0*_uy[idx]));
                        _q.f[Q::IndexF(idx, 5)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }

                    //  On ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 - 3.0*_uy[idx])*(8.0*_q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])
                            + _uz[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 14)])
                            - 24.0*_eps
                        )/(12.0*(1.0 + 3.0*_uy[idx]));
                        _q.f[Q::IndexF(idx, 2)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  On zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 + 3.0*_uz[idx])*(8.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 14)])
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])
                            + _uy[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])
                            - 24.0*_eps
                        )/(12.0*(1.0 - 3.0*_uz[idx]));
                        _q.f[Q::IndexF(idx, 6)] = rho0;
                        _q.f[Q::IndexF(idx, 10)] = rho0;
                        _q.f[Q::IndexF(idx, 11)] = rho0;
                        _q.f[Q::IndexF(idx, 12)] = rho0;
                        _q.f[Q::IndexF(idx, 13)] = rho0;
                    }

                    //  On zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixQ) {
                        T rho0 = (
                            (1.0 - 3.0*_uz[idx])*(8.0*_q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 13)])
                            + _uy[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])
                            - 24.0*_eps
                        )/(12.0*(1.0 + 3.0*_uz[idx]));
                        _q.f[Q::IndexF(idx, 3)] = rho0;
                        _q.f[Q::IndexF(idx, 7)] = rho0;
                        _q.f[Q::IndexF(idx, 8)] = rho0;
                        _q.f[Q::IndexF(idx, 9)] = rho0;
                        _q.f[Q::IndexF(idx, 14)] = rho0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D2Q9
        template<class T, class P, class Q>
        void BoundaryConditionSetiRho(P& _p, Q& _q, const T *_rho, const T *_ux, const T *_uy, const T *_tem, const int *_bctypef, const int *_bctypeg, T _eps = T()) {
            int idx, idxbc; 
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                idx = _q.Index(0, j);
                idxbc = j + _q.offsetxmin;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])/(2.0*(1.0 + 3.0*_ux[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])/3.0
                            + _uy[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])/2.0
                        )/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                }

                //  On xmax
                idx = _q.Index(_q.nx - 1, j);
                idxbc = j + _q.offsetxmax;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 7)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_ux[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])/3.0
                             + _uy[idx]*(_q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])/2.0
                        )/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                idx = _q.Index(i, 0);
                idxbc = i + _q.offsetymin;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 6)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])/(2.0*(1.0 + 3.0*_uy[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])/3.0
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])/2.0
                        )/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                }

                //  On ymax
                idx = _q.Index(i, _q.ny - 1);
                idxbc = i + _q.offsetymax;
                if (_bctypef[idxbc] == SetiRho && (_bctypeg[idxbc] == SetiT || _bctypeg[idxbc] == SetiQ)) {    
                    T rho0 = -(4.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)])/3.0;
                    T flux0 = T();
                    if (_bctypeg[idxbc] == FixT) {
                        flux0 = _tem[idx]*_ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])/(2.0*(1.0 - 3.0*_uy[idx])*_rho[idx]);
                    } else if (_bctypeg[idxbc] == FixQ) {
                        flux0 = -_tem[idx]*(
                            (4.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])/3.0 
                            + _ux[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 7)])/2.0
                        )/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                    }
                    T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                    _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                    _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                }
            }
        }
    
        //  Function of setting boundary condition set iRho and iT or iQ of AAD for D3Q15
        template<class T, class P, class Q>
        void BoundaryConditionSetiRho(P& _p, Q& _q, const T *_rho, const T *_ux, const T *_uy, const T *_uz, const T *_tem, const int *_bctypef, const int *_bctypeg, T _eps = T()) {
            int idx, idxbc;

            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    //  On xmin
                    idx = _p.Index(0, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmin;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 1)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 9)] + _p.f[P::IndexF(idx, 10)] + _p.f[P::IndexF(idx, 12)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _uy[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/4.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 12)])/6.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/4.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_ux[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 4)] = _p.f[P::IndexF(idx, 1)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 11)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 13)] = _p.f[P::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 14)] = _p.f[P::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                    }

                    //  On xmax
                    idx = _p.Index(_p.nx - 1, j, k);
                    idxbc = _p.IndexBCx(j, k) + _p.offsetxmax;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 4)] + _p.f[P::IndexF(idx, 8)] + _p.f[P::IndexF(idx, 11)] + _p.f[P::IndexF(idx, 13)] + _p.f[P::IndexF(idx, 14)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _uy[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] - _q.f[Q::IndexF(idx, 14)])/4.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_ux[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/6.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] - _q.f[Q::IndexF(idx, 14)])/4.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_ux[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_ux[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 1)] = _p.f[P::IndexF(idx, 4)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 9)] = _p.f[P::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 10)] = _p.f[P::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 12)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                    }
                }
            }
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    //  On ymin
                    idx = _p.Index(i, 0, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymin;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 2)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)] + _p.f[P::IndexF(idx, 10)] + _p.f[P::IndexF(idx, 13)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/4.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 13)])/6.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/4.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uy[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 5)] = _p.f[P::IndexF(idx, 2)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 9)] = _p.f[P::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 11)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 12)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 14)] = _p.f[P::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                    }

                    //  On ymax
                    idx = _p.Index(i, _p.ny - 1, k);
                    idxbc = _p.IndexBCy(k, i) + _p.offsetymax;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 5)] + _p.f[P::IndexF(idx, 9)] + _p.f[P::IndexF(idx, 11)] + _p.f[P::IndexF(idx, 12)] + _p.f[P::IndexF(idx, 14)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _uz[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/4.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_uy[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/6.0
                                + _uz[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/4.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_uy[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uy[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 2)] = _p.f[P::IndexF(idx, 5)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 10)] = _p.f[P::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 13)] = _p.f[P::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                    }
                }
            }
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    //  On zmin
                    idx = _p.Index(i, j, 0);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmin;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 3)] + _p.f[P::IndexF(idx, 7)] + _p.f[P::IndexF(idx, 8)] + _p.f[P::IndexF(idx, 9)] + _p.f[P::IndexF(idx, 14)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/4.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 14)])/6.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/4.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 - 3.0*_uz[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 6)] = _p.f[P::IndexF(idx, 3)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 10)] = _p.f[P::IndexF(idx, 14)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 11)] = _p.f[P::IndexF(idx, 7)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 12)] = _p.f[P::IndexF(idx, 8)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 13)] = _p.f[P::IndexF(idx, 9)] + rho0 + flux0 + obj0;
                    }

                    //  On zmax
                    idx = _p.Index(i, j, _p.nz - 1);
                    idxbc = _p.IndexBCz(i, j) + _p.offsetzmax;
                    if (_bctypef[idxbc] == OUTLET && (_bctypeg[idxbc] == FixT || _bctypeg[idxbc] == FixQ)) {
                        T rho0 = -(8.0*_p.f[P::IndexF(idx, 6)] + _p.f[P::IndexF(idx, 10)] + _p.f[P::IndexF(idx, 11)] + _p.f[P::IndexF(idx, 12)] + _p.f[P::IndexF(idx, 13)])/6.0;
                        T flux0 = T();
                        if (_bctypeg[idxbc] == FixT) {
                            flux0 = _tem[idx]*(
                                _ux[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 13)])/4.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/4.0
                            )/(_rho[idx]*(1.0 - 3.0*_uz[idx]));
                        } else if (_bctypeg[idxbc] == FixQ) {
                            flux0 = -_tem[idx]*(
                                (8.0*_q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/6.0
                                + _ux[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 13)])/4.0
                                + _uy[idx]*(_q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/4.0
                            )/(_rho[idx]*(1.0 + 3.0*_uz[idx]));
                        }
                        T obj0 = _eps*2.0*_tem[idx]/((1.0 + 3.0*_uz[idx])*_rho[idx]);
                        _p.f[P::IndexF(idx, 3)] = _p.f[P::IndexF(idx, 6)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 7)] = _p.f[P::IndexF(idx, 11)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 8)] = _p.f[P::IndexF(idx, 12)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 9)] = _p.f[P::IndexF(idx, 13)] + rho0 + flux0 + obj0;
                        _p.f[P::IndexF(idx, 14)] = _p.f[P::IndexF(idx, 10)] + rho0 + flux0 + obj0;
                    }
                }
            }
        }
    
        //  Function of getting sensitivity of temperature at heat source for D2Q9
        template<class T, class Q>
        void SensitivityTemperatureAtHeatSource(
            const T *_ux, const T *_uy, const T *_imx, const T *_imy,
            Q& _q, const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_g, const T *_ig,
            T *_dfds, const T *_diffusivity, const T *_dads, const T *_dkds, const T *_qnbc, const int *_bctype
        ) {
            int idx, idxbc;

            //  Brinkman term and diffusivity term
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    idx = _q.Index(i, j);
                    _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx]);
                    T sumg = T();
                    for (int c = 0; c < Q::nc; ++c) {
                        sumg += _ig[Q::IndexF(idx, c)]*_g[Q::IndexF(idx, c)];
                    }
                    _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx])));
                }
            }

            //  Boundary term 
            for (int j = 0; j < _q.ny; ++j) {
                //  Along xmin
                idx = _q.Index(0, j);
                idxbc = j + _q.offsetxmin;
                if (_bctype[idxbc] == FixQ) {
                    _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                        (1.0 + 3.0*_ux[idx])*(-6.0 + 4.0*_ig[Q::IndexF(idx, 1)] + _ig[Q::IndexF(idx, 5)] + _ig[Q::IndexF(idx, 8)])
                        + 3.0*_uy[idx]*(_ig[Q::IndexF(idx, 5)] - _ig[Q::IndexF(idx, 8)])
                    )/(36.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                }

                //  Along xmax
                idx = _q.Index(_q.nx - 1, j);
                idxbc = j + _q.offsetxmax;
                if (_bctype[idxbc] == FixQ) {
                    _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                        (1.0 - 3.0*_ux[idx])*(-6.0 + 4.0*_ig[Q::IndexF(idx, 3)] + _ig[Q::IndexF(idx, 6)] + _ig[Q::IndexF(idx, 7)])
                        + 3.0*_uy[idx]*(_ig[Q::IndexF(idx, 6)] - _ig[Q::IndexF(idx, 7)])
                    )/(36.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                //  Along ymin
                idx = _q.Index(i, 0);
                idxbc = i + _q.offsetymin;
                if (_bctype[idxbc] == FixQ) {
                    _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                        (1.0 + 3.0*_uy[idx])*(-6.0 + 4.0*_ig[Q::IndexF(idx, 2)] + _ig[Q::IndexF(idx, 5)] + _ig[Q::IndexF(idx, 6)])
                        + 3.0*_ux[idx]*(_ig[Q::IndexF(idx, 5)] - _ig[Q::IndexF(idx, 6)])
                    )/(36.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                }

                //  Along ymax
                idx = _q.Index(i, _q.ny - 1);
                idxbc = i + _q.offsetymax;
                if (_bctype[idxbc] == FixQ) {
                    _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                        (1.0 - 3.0*_uy[idx])*(-6.0 + 4.0*_ig[Q::IndexF(idx, 4)] + _ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 8)])
                        + 3.0*_ux[idx]*(_ig[Q::IndexF(idx, 8)] - _ig[Q::IndexF(idx, 7)])
                    )/(36.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                }
            }
        }

        //  Function of getting sensitivity of temperature at heat source for D3Q15
        template<class T, class Q>
        void SensitivityTemperatureAtHeatSource(
            const T *_ux, const T *_uy, const T *_uz, const T *_imx, const T *_imy, const T *_imz,
            Q& _q, const T *_tem, const T *_item, const T *_iqx, const T *_iqy, const T *_iqz, const T *_g, const T *_ig,
            T *_dfds, const T *_diffusivity, const T *_dads, const T *_dkds, const T *_qnbc, const int *_bctype
        ) {
            int idx, idxbc;

            //  Brinkman term and diffusivity term
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        idx = _q.Index(i, j, k);
                        _dfds[idx] += 3.0*_dads[idx]*(_ux[idx]*_imx[idx] + _uy[idx]*_imy[idx] + _uz[idx]*_imz[idx]);
                        T sumg = T();
                        for (int c = 0; c < Q::nc; ++c) {
                            sumg += _ig[Q::IndexF(idx, c)]*_g[Q::IndexF(idx, c)];
                        }
                        _dfds[idx] += -3.0/pow(3.0*_diffusivity[idx] + 0.5, 2.0)*_dkds[idx]*(sumg - _tem[idx]*(_item[idx] + 3.0*(_ux[idx]*_iqx[idx] + _uy[idx]*_iqy[idx] + _uz[idx]*_iqz[idx])));
                    }
                }
            }

            //  Boundary term 
            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  Along xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 + 3.0*_ux[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 1)] + _ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 9)] + _ig[Q::IndexF(idx, 10)] + _ig[Q::IndexF(idx, 12)])
                            + 3.0*_uy[idx]*(_ig[Q::IndexF(idx, 7)] - _ig[Q::IndexF(idx, 9)] + _ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 12)])
                            + 3.0*_uz[idx]*(_ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 9)] - _ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 12)])
                        )/(72.0*(1.0 - 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                    }

                    //  Along xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 - 3.0*_ux[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 4)] + _ig[Q::IndexF(idx, 8)] + _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 13)] + _ig[Q::IndexF(idx, 14)])
                            + 3.0*_uy[idx]*(_ig[Q::IndexF(idx, 8)] - _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 13)] - _ig[Q::IndexF(idx, 14)])
                            + 3.0*_uz[idx]*(_ig[Q::IndexF(idx, 8)] - _ig[Q::IndexF(idx, 11)] - _ig[Q::IndexF(idx, 13)] + _ig[Q::IndexF(idx, 14)])
                        )/(72.0*(1.0 + 3.0*_ux[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  Along ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 + 3.0*_uy[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 2)] + _ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 8)] + _ig[Q::IndexF(idx, 10)] + _ig[Q::IndexF(idx, 13)])
                            + 3.0*_uz[idx]*(_ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 8)] - _ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 13)])
                            + 3.0*_ux[idx]*(_ig[Q::IndexF(idx, 7)] - _ig[Q::IndexF(idx, 8)] + _ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 13)])
                        )/(72.0*(1.0 - 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                    }

                    //  Along ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 - 3.0*_uy[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 5)] + _ig[Q::IndexF(idx, 9)] + _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 12)] + _ig[Q::IndexF(idx, 14)])
                            + _uz[idx]*(_ig[Q::IndexF(idx, 9)] - _ig[Q::IndexF(idx, 11)] - _ig[Q::IndexF(idx, 12)] + _ig[Q::IndexF(idx, 14)])
                            + _ux[idx]*(_ig[Q::IndexF(idx, 9)] - _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 12)] - _ig[Q::IndexF(idx, 14)])
                        )/(72.0*(1.0 + 3.0*_uy[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  Along zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 + 3.0*_uz[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 3)] + _ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 8)] + _ig[Q::IndexF(idx, 9)] + _ig[Q::IndexF(idx, 14)])
                            + _ux[idx]*(_ig[Q::IndexF(idx, 7)] - _ig[Q::IndexF(idx, 8)] + _ig[Q::IndexF(idx, 9)] - _ig[Q::IndexF(idx, 14)])
                            + _uy[idx]*(_ig[Q::IndexF(idx, 7)] + _ig[Q::IndexF(idx, 8)] - _ig[Q::IndexF(idx, 9)] - _ig[Q::IndexF(idx, 14)])
                        )/(72.0*(1.0 - 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                    }

                    //  Along zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixQ) {
                        _dfds[idx] += _qnbc[idxbc]*_dkds[idx]*(
                            (1.0 - 3.0*_uz[idx])*(-12.0 + 8.0*_ig[Q::IndexF(idx, 6)] + _ig[Q::IndexF(idx, 10)] + _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 12)] + _ig[Q::IndexF(idx, 13)])
                            + _ux[idx]*(_ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 11)] + _ig[Q::IndexF(idx, 12)] - _ig[Q::IndexF(idx, 13)])
                            + _uy[idx]*(_ig[Q::IndexF(idx, 10)] - _ig[Q::IndexF(idx, 11)] - _ig[Q::IndexF(idx, 12)] + _ig[Q::IndexF(idx, 13)])
                        )/(72.0*(1.0 + 3.0*_uz[idx])*pow(_diffusivity[idx], 2.0));
                    }
                }
            }
        }
    }
}