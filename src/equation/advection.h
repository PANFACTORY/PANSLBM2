//*****************************************************************************
//  Title       :   src/equation/advection.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/08/02
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cassert>
#include "navierstokes.h"

namespace PANSLBM2 {
    namespace {
        const int FixT = 1;
        const int FixQ = 2;
    }

    namespace AD {
        //  Function of updating macroscopic values of AD for 2D
        template<class T, class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T _ux, T _uy, const T *_g, T _omegag, int _idx) {
            _tem = T();
            _qx = T();
            _qy = T();
            for (int c = 0; c < Q::nc; ++c) {
                _tem += _g[Q::IndexF(_idx, c)];
                _qx += Q::cx[c]*_g[Q::IndexF(_idx, c)];
                _qy += Q::cy[c]*_g[Q::IndexF(_idx, c)];
            }
            _qx = (1.0 - 0.5*_omegag)*(_qx - _tem*_ux);
            _qy = (1.0 - 0.5*_omegag)*(_qy - _tem*_uy);
        }

        //  Function of updating macroscopic values of AD for 3D
        template<class T, class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T &_qz, T _ux, T _uy, T _uz, const T *_g, T _omegag, int _idx) {
            _tem = T();
            _qx = T();
            _qy = T();
            _qz = T();
            for (int c = 0; c < Q::nc; ++c) {
                _tem += _g[Q::IndexF(_idx, c)];
                _qx += Q::cx[c]*_g[Q::IndexF(_idx, c)];
                _qy += Q::cy[c]*_g[Q::IndexF(_idx, c)];
                _qz += Q::cz[c]*_g[Q::IndexF(_idx, c)];
            }
            _qx = (1.0 - 0.5*_omegag)*(_qx - _tem*_ux);
            _qy = (1.0 - 0.5*_omegag)*(_qy - _tem*_uy);
            _qz = (1.0 - 0.5*_omegag)*(_qz - _tem*_uz);
        }

        //  Function of getting equilibrium of AD for 2D
        template<class T, class Q>
        T Equilibrium(T _tem, T _ux, T _uy, int _c) {
            T ciu = Q::cx[_c]*_ux + Q::cy[_c]*_uy;
            return Q::ei[_c]*_tem*(1.0 + 3.0*ciu);
        }

        //  Function of getting equilibrium of AD for 3D
        template<class T, class Q>
        T Equilibrium(T _tem, T _ux, T _uy, T _uz, int _c) {
            T ciu = Q::cx[_c]*_ux + Q::cy[_c]*_uy + Q::cz[_c]*_uz;
            return Q::ei[_c]*_tem*(1.0 + 3.0*ciu);
        }

        //  Function of applying external force of AD with natural convection for 2D
        template<class T, class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _tem0, T *_f, int _idx) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] += 3.0*P::ei[c]*(P::cx[c]*_gx + P::cy[c]*_gy)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with natural convection for 3D
        template<class T, class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _gz, T _tem0, T *_f, int _idx) {
            for (int c = 0; c < P::nc; ++c) {
                _f[P::IndexF(_idx, c)] += 3.0*P::ei[c]*(P::cx[c]*_gx + P::cy[c]*_gy + P::cz[c]*_gz)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with heat exchange for 2D/3D
        template<class T, class Q>
        void ExternalForceHeatExchange(T _tem, T *_g, const T *_beta, int _idx) {
            for (int c = 0; c < Q::nc; ++c) {
                _g[Q::IndexF(_idx, c)] += Q::ei[c]*_beta[_idx]*(1.0 - _tem)/(1.0 + _beta[_idx]);
            }
        }

        //  Function of Update macro, Collide and Stream of force convection for 2D
        template<class T, class P, class Q>
        void Macro_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of force convection for 3D
        template<class T, class P, class Q>
        void Macro_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of natural convection for 2D
        template<class T, class P, class Q>
        void Macro_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with natural convection
                    ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of natural convection for 3D
        template<class T, class P, class Q>
        void Macro_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with natural convection
                        ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }
    
        //  Function of Update macro, Collide and Stream of Brinkman and heat exchange for 2D
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_HeatExchange(
            P& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with Brinkman and heat exchange
                    NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    ExternalForceHeatExchange<T, Q>(tem, _q.f, _beta, idx);
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and heat exchange for 3D
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_HeatExchange(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with Brinkman and heat exchange
                        NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        ExternalForceHeatExchange<T, Q>(tem, _q.f, _beta, idx);
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and force convection for 2D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with Brinkman model
                    NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and force convection for 3D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with Brinkman model
                        NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and force convection for 2D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with Brinkman model
                    NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and force convection for 3D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_ForceConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with Brinkman model
                        NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and natural convection for 2D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with Brinkman model and natural convection
                    ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                    NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and natural convection for 3D (diffusivity constant)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), omegag = 1.0/(3.0*_diffusivity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with Brinkman model and natural convection
                        ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                        NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and natural convection for 2D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    int idx = _p.Index(i, j);
                    T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                    //  Update macro
                    T rho, ux, uy;
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    T tem, qx, qy;
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  External force with Brinkman model and natural convection
                    ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                    NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                    NS::Macro<T, P>(rho, ux, uy, _p.f, idx);
                    Macro<T, Q>(tem, qx, qy, ux, uy, _q.f, omegag, idx);

                    //  Save macro if need
                    if (_issave) {
                        _rho[idx] = rho;
                        _ux[idx] = ux;
                        _uy[idx] = uy;
                        _tem[idx] = tem;
                        _qx[idx] = qx;
                        _qy[idx] = qy;
                    }

                    //  Collide and stream
                    for (int c = 0; c < P::nc; ++c) {
                        int idxstream = _p.Index(i + P::cx[c], j + P::cy[c]);
                        _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, c);
                    }
                    for (int c = 0; c < Q::nc; ++c) {
                        int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c]);
                        _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, c);
                    }
                }
            }
        }

        //  Function of Update macro, Collide and Stream of Brinkman and natural convection for 3D (diffusivity heterogenious)
        template<class T, class P, class Q>
        void Macro_Brinkman_Collide_Stream_NaturalConvection(
            P& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5);
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    for (int k = 0; k < _p.nz; ++k) {
                        int idx = _p.Index(i, j, k);
                        T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5);

                        //  Update macro
                        T rho, ux, uy, uz;
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        T tem, qx, qy, qz;
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  External force with Brinkman model and natural convection
                        ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                        NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                        NS::Macro<T, P>(rho, ux, uy, uz, _p.f, idx);
                        Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f, omegag, idx);

                        //  Save macro if need
                        if (_issave) {
                            _rho[idx] = rho;
                            _ux[idx] = ux;
                            _uy[idx] = uy;
                            _uz[idx] = uz;
                            _tem[idx] = tem;
                            _qx[idx] = qx;
                            _qy[idx] = qy;
                            _qz[idx] = qz;
                        }

                        //  Collide and stream
                        for (int c = 0; c < P::nc; ++c) {
                            int idxstream = _p.Index(i + P::cx[c], j + P::cy[c], k + P::cz[c]);
                            _p.fnext[P::IndexF(idxstream, c)] = (1.0 - omegaf)*_p.f[P::IndexF(idx, c)] + omegaf*NS::Equilibrium<T, P>(rho, ux, uy, uz, c);
                        }
                        for (int c = 0; c < Q::nc; ++c) {
                            int idxstream = _q.Index(i + Q::cx[c], j + Q::cy[c], k + Q::cz[c]);
                            _q.fnext[Q::IndexF(idxstream, c)] = (1.0 - omegag)*_q.f[Q::IndexF(idx, c)] + omegag*Equilibrium<T, Q>(tem, ux, uy, uz, c);
                        }
                    }
                }
            }
        }

        //  Function of setting initial condition of AD for 2D
        template<class T, class Q>
        void InitialCondition(Q& _q, const T *_tem, const T *_ux, const T *_uy) {
            for (int i = 0; i < _q.nx; ++ i) {
                for (int j = 0; j < _q.ny; ++j) {
                    int idx = _q.Index(i, j);
                    for (int c = 0; c < Q::nc; ++c) {
                        _q.f[Q::IndexF(idx, c)] = Equilibrium<T, Q>(_tem[idx], _ux[idx], _uy[idx], c);
                    }
                }
            }
        }

        //  Function of setting initial condition of AD for 3D
        template<class T, class Q>
        void InitialCondition(Q& _q, const T *_tem, const T *_ux, const T *_uy, const T *_uz) {
            for (int i = 0; i < _q.nx; ++ i) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        int idx = _q.Index(i, j, k);
                        for (int c = 0; c < Q::nc; ++c) {
                            _q.f[Q::IndexF(idx, c)] = Equilibrium<T, Q>(_tem[idx], _ux[idx], _uy[idx], _uz[idx], c);
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D2Q9
        template<class T, class Q>
        void BoundaryConditionSetT(Q& _q, const T *_tembc, const T *_ux, const T *_uy, const int *_bctype) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == FixT) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)])/(1.0 + 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == FixT) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 8)])/(1.0 - 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == FixT) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)])/(1.0 + 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == FixT) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)])/(1.0 - 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D3Q15
        template<class T, class Q>
        void BoundaryConditionSetT(Q& _q, const T *_tembc, const T *_ux, const T *_uy, const T *_uz, const int *_bctype) {
            int idx, idxbc;

            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  On xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 13)] - _q.f[Q::IndexF(idx, 14)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 12)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  On ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 14)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 13)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  On zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 6)] - _q.f[Q::IndexF(idx, 10)] - _q.f[Q::IndexF(idx, 11)] - _q.f[Q::IndexF(idx, 12)] - _q.f[Q::IndexF(idx, 13)])/(1.0 + 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }

                    //  On zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixT) {
                        T tem0 = 6.0*(_tembc[idxbc] - _q.f[Q::IndexF(idx, 0)] - _q.f[Q::IndexF(idx, 1)] - _q.f[Q::IndexF(idx, 2)] - _q.f[Q::IndexF(idx, 3)] - _q.f[Q::IndexF(idx, 4)] - _q.f[Q::IndexF(idx, 5)] - _q.f[Q::IndexF(idx, 7)] - _q.f[Q::IndexF(idx, 8)] - _q.f[Q::IndexF(idx, 9)] - _q.f[Q::IndexF(idx, 14)])/(1.0 - 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity constant)
        template<class T, class Q>
        void BoundaryConditionSetQ(Q& _q, const T *_qnbc, const T *_ux, const T *_uy, T _diffusivity, const int *_bctype) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == FixQ) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == FixQ) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == FixQ) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == FixQ) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity constant)
        template<class T, class Q>
        void BoundaryConditionSetQ(Q& _q, const T *_qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, const int *_bctype) {
            int idx, idxbc;

            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  On xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  On ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  On zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }

                    //  On zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity heterogenious)
        template<class T, class Q>
        void BoundaryConditionSetQ(Q& _q, const T *_qnbc, const T *_ux, const T *_uy, const T *_diffusivity, const int *_bctype) {
            for (int j = 0; j < _q.ny; ++j) {
                //  On xmin
                if (_bctype[j + _q.offsetxmin] == FixQ) {
                    int idx = _q.Index(0, j), idxbc = j + _q.offsetxmin;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }

                //  On xmax
                if (_bctype[j + _q.offsetxmax] == FixQ) {
                    int idx = _q.Index(_q.nx - 1, j), idxbc = j + _q.offsetxmax;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                    _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }

            for (int i = 0; i < _q.nx; ++i) {
                //  On ymin
                if (_bctype[i + _q.offsetymin] == FixQ) {
                    int idx = _q.Index(i, 0), idxbc = i + _q.offsetymin;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                }

                //  On ymax
                if (_bctype[i + _q.offsetymax] == FixQ) {
                    int idx = _q.Index(i, _q.ny - 1), idxbc = i + _q.offsetymax;
                    T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                    _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                    _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity heterogenious)
        template<class T, class Q>
        void BoundaryConditionSetQ(Q& _q, const T *_qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, const int *_bctype) {
            int idx, idxbc;

            for (int j = 0; j < _q.ny; ++j) {
                for (int k = 0; k < _q.nz; ++k) {
                    //  On xmin
                    idx = _q.Index(0, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 4)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 13)] + _q.f[Q::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On xmax
                    idx = _q.Index(_q.nx - 1, j, k);
                    idxbc = _q.IndexBCx(j, k) + _q.offsetxmax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 1)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int k = 0; k < _q.nz; ++k) {
                for (int i = 0; i < _q.nx; ++i) {
                    //  On ymin
                    idx = _q.Index(i, 0, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 5)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }

                    //  On ymax
                    idx = _q.Index(i, _q.ny - 1, k);
                    idxbc = _q.IndexBCy(k, i) + _q.offsetymax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 2)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }
                }
            }
            for (int i = 0; i < _q.nx; ++i) {
                for (int j = 0; j < _q.ny; ++j) {
                    //  On zmin
                    idx = _q.Index(i, j, 0);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmin;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 6)] + _q.f[Q::IndexF(idx, 10)] + _q.f[Q::IndexF(idx, 11)] + _q.f[Q::IndexF(idx, 12)] + _q.f[Q::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                    }

                    //  On zmax
                    idx = _q.Index(i, j, _q.nz - 1);
                    idxbc = _q.IndexBCz(i, j) + _q.offsetzmax;
                    if (_bctype[idxbc] == FixQ) {
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc[idxbc] + _q.f[Q::IndexF(idx, 3)] + _q.f[Q::IndexF(idx, 7)] + _q.f[Q::IndexF(idx, 8)] + _q.f[Q::IndexF(idx, 9)] + _q.f[Q::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                        _q.f[Q::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                        _q.f[Q::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        _q.f[Q::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                    }
                }
            }
        }
    }
}