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
    namespace AD {
        //  Function of updating macroscopic values of AD for 2D
        template<class T, template<class>class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T _ux, T _uy, const T *_g0, const T *_g, T _omegag, int _idx) {
            _tem = _g0[_idx];
            _qx = T();
            _qy = T();
            for (int c = 1; c < Q<T>::nc; ++c) {
                T g = _g[Q<T>::IndexF(_idx, c)];
                _tem += g;
                _qx += Q<T>::cx[c]*g;
                _qy += Q<T>::cy[c]*g;
            }
            T coef = 1.0 - 0.5*_omegag;
            _qx = coef*(_qx - _tem*_ux);
            _qy = coef*(_qy - _tem*_uy);
        }

        //  Function of updating macroscopic values of AD for 3D
        template<class T, template<class>class Q>
        void Macro(T &_tem, T &_qx, T &_qy, T &_qz, T _ux, T _uy, T _uz, const T *_g0, const T *_g, T _omegag, int _idx) {
            _tem = _g0[_idx];
            _qx = T();
            _qy = T();
            _qz = T();
            for (int c = 1; c < Q<T>::nc; ++c) {
                T g = _g[Q<T>::IndexF(_idx, c)];
                _tem += g;
                _qx += Q<T>::cx[c]*g;
                _qy += Q<T>::cy[c]*g;
                _qz += Q<T>::cz[c]*g;
            }
            T coef = 1.0 - 0.5*_omegag;
            _qx = coef*(_qx - _tem*_ux);
            _qy = coef*(_qy - _tem*_uy);
            _qz = coef*(_qz - _tem*_uz);
        }

        //  Function of getting equilibrium of AD for 2D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _tem, T _ux, T _uy) {
            for (int c = 0; c < Q<T>::nc; ++c) {
                T ciu = Q<T>::cx[c]*_ux + Q<T>::cy[c]*_uy;
                _geq[c] = Q<T>::ei[c]*_tem*(1.0 + 3.0*ciu);
            }
        }

        //  Function of getting equilibrium of AD for 3D
        template<class T, template<class>class Q>
        void Equilibrium(T *_geq, T _tem, T _ux, T _uy, T _uz) {
            for (int c = 0; c < Q<T>::nc; ++c) {
                T ciu = Q<T>::cx[c]*_ux + Q<T>::cy[c]*_uy + Q<T>::cz[c]*_uz;
                _geq[c] = Q<T>::ei[c]*_tem*(1.0 + 3.0*ciu);
            }
        }

        //  Function of applying external force of AD with natural convection for 2D
        template<class T, template<class>class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _tem0, T *_f, int _idx) {
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*P<T>::ei[c]*(P<T>::cx[c]*_gx + P<T>::cy[c]*_gy)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with natural convection for 3D
        template<class T, template<class>class P>
        void ExternalForceNaturalConvection(T _tem, T _gx, T _gy, T _gz, T _tem0, T *_f, int _idx) {
            for (int c = 1; c < P<T>::nc; ++c) {
                _f[P<T>::IndexF(_idx, c)] += 3.0*P<T>::ei[c]*(P<T>::cx[c]*_gx + P<T>::cy[c]*_gy + P<T>::cz[c]*_gz)*(_tem - _tem0);
            }
        }

        //  Function of applying external force of AD with heat exchange for 2D/3D
        template<class T, template<class>class Q>
        void ExternalForceHeatExchange(T _tem, T _beta, T *_g0, T *_g, int _idx) {
            T coef = _beta*(1.0 - _tem)/(1.0 + _beta);
            _g0[_idx] += Q<T>::ei[0]*coef;
            for (int c = 1; c < Q<T>::nc; ++c) {
                _g[Q<T>::IndexF(_idx, c)] += Q<T>::ei[c]*coef;
            }
        }

        //  Function of Update macro and Collide of force convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
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

        //  Function of Update macro and Collide of force convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

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

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of natural convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T _diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc]; 
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);

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

        //  Function of Update macro and Collide of natural convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, T _diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);

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

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c); 
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }
    
        //  Function of Update macro and Collide of Brinkman and heat exchange for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
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

        //  Function of Update macro and Collide of Brinkman and heat exchange for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideHeatExchange(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_beta, T _diffusivity, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc];
            T omegag = 1.0/(3.0*_diffusivity + 0.5), iomegag = 1.0 - omegag, geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman and heat exchange
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                ExternalForceHeatExchange<T, Q>(tem, _beta[idx], _q.f0, _q.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

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

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {            
                    int idxf = Q<T>::IndexF(idx, c);      
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and force convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
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

        //  Function of Update macro and Collide of Brinkman and force convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideForceConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

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

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c);
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of Update macro and Collide of Brinkman and natural convection for 2D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, const T *_diffusivity, 
            T _gx, T _gy, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy;
                NS::Macro<T, P>(rho, ux, uy, _p.f0, _p.f, idx);
                T tem, qx, qy;
                Macro<T, Q>(tem, qx, qy, ux, uy, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, _alpha[idx], _p.f, idx);
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

        //  Function of Update macro and Collide of Brinkman and natural convection for 3D
        template<class T, template<class>class P, template<class>class Q>
        void MacroBrinkmanCollideNaturalConvection(
            P<T>& _p, T *_rho, T *_ux, T *_uy, T *_uz, const T *_alpha, T _viscosity,
            Q<T>& _q, T *_tem, T *_qx, T *_qy, T *_qz, const T *_diffusivity, 
            T _gx, T _gy, T _gz, T _tem0, bool _issave = false
        ) {
            T omegaf = 1.0/(3.0*_viscosity + 0.5), iomegaf = 1.0 - omegaf, feq[P<T>::nc], geq[Q<T>::nc];
#pragma omp parallel for
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                T omegag = 1.0/(3.0*_diffusivity[idx] + 0.5), iomegag = 1.0 - omegag;

                //  Update macro
                T rho, ux, uy, uz;
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                T tem, qx, qy, qz;
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

                //  External force with Brinkman model and natural convection
                ExternalForceNaturalConvection<T, P>(tem, _gx, _gy, _gz, _tem0, _p.f, idx);
                NS::ExternalForceBrinkman<T, P>(rho, ux, uy, uz, _alpha[idx], _p.f, idx);
                NS::Macro<T, P>(rho, ux, uy, uz, _p.f0, _p.f, idx);
                Macro<T, Q>(tem, qx, qy, qz, ux, uy, uz, _q.f0, _q.f, omegag, idx);

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

                //  Collide
                NS::Equilibrium<T, P>(feq, rho, ux, uy, uz);
                _p.f0[idx] = iomegaf*_p.f0[idx] + omegaf*feq[0];
                for (int c = 1; c < P<T>::nc; ++c) {
                    int idxf = P<T>::IndexF(idx, c);
                    _p.f[idxf] = iomegaf*_p.f[idxf] + omegaf*feq[c];
                }
                Equilibrium<T, Q>(geq, tem, ux, uy, uz);
                _q.f0[idx] = iomegag*_q.f0[idx] + omegag*geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    int idxf = Q<T>::IndexF(idx, c); 
                    _q.f[idxf] = iomegag*_q.f[idxf] + omegag*geq[c];
                }
            }
        }

        //  Function of setting initial condition of AD for 2D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_tem, const T *_ux, const T *_uy) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _tem[idx], _ux[idx], _uy[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting initial condition of AD for 3D
        template<class T, template<class>class Q>
        void InitialCondition(Q<T>& _q, const T *_tem, const T *_ux, const T *_uy, const T *_uz) {
            T geq[Q<T>::nc];
            for (int idx = 0; idx < _q.nxyz; ++idx) {
                Equilibrium<T, Q>(geq, _tem[idx], _ux[idx], _uy[idx], _uz[idx]);
                _q.f0[idx] = geq[0];
                for (int c = 1; c < Q<T>::nc; ++c) {
                    _q.f[Q<T>::IndexF(idx, c)] = geq[c];
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D2Q9
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetT(Q<T>& _q, Fv _tembc, const T *_ux, const T *_uy, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        T tem0 = 6.0*(_tembc(0 + _q.offsetx, j + _q.offsety) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        T tem0 = 6.0*(_tembc((_q.nx - 1) + _q.offsetx, j + _q.offsety) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        T tem0 = 6.0*(_tembc(i + _q.offsetx, 0 + _q.offsety) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) { 
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        T tem0 = 6.0*(_tembc(i + _q.offsetx, (_q.ny - 1) + _q.offsety) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set T of AD for D3Q15
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetT(Q<T>& _q, Fv _tembc, const T *_ux, const T *_uy, const T *_uz, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            T tem0 = 6.0*(_tembc(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 13)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*(_tembc((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 12)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*(_tembc(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*(_tembc(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*(_tembc(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 6)] - _q.f[Q<T>::IndexF(idx, 10)] - _q.f[Q<T>::IndexF(idx, 11)] - _q.f[Q<T>::IndexF(idx, 12)] - _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*(_tembc(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz) - _q.f0[idx] - _q.f[Q<T>::IndexF(idx, 1)] - _q.f[Q<T>::IndexF(idx, 2)] - _q.f[Q<T>::IndexF(idx, 3)] - _q.f[Q<T>::IndexF(idx, 4)] - _q.f[Q<T>::IndexF(idx, 5)] - _q.f[Q<T>::IndexF(idx, 7)] - _q.f[Q<T>::IndexF(idx, 8)] - _q.f[Q<T>::IndexF(idx, 9)] - _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        }
                    }
                }
            }
        }
    
        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, T _diffusivity, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(0 + _q.offsetx, j + _q.offsety) + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety) + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, 0 + _q.offsety) + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety) + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity constant)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, T _diffusivity, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity))*_qnbc(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        }
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D2Q9 (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_diffusivity, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype(0 + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(0, j);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(0 + _q.offsetx, j + _q.offsety) + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 7)])/(1.0 - 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On xmax
            if (_q.PEx == _q.mx - 1) {
                for (int j = 0; j < _q.ny; ++j) {
                    if (_bctype((_q.nx - 1) + _q.offsetx, j + _q.offsety)) {
                        int idx = _q.Index(_q.nx - 1, j);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety) + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 + 3.0*_ux[idx]);
                        _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymin
            if (_q.PEy == 0) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, 0 + _q.offsety)) {
                        int idx = _q.Index(i, 0);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, 0 + _q.offsety) + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)])/(1.0 - 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx])/36.0;
                    }
                }
            }
            //  On ymax
            if (_q.PEy == _q.my - 1) {
                for (int i = 0; i < _q.nx; ++i) {
                    if (_bctype(i + _q.offsetx, (_q.ny - 1) + _q.offsety)) {
                        int idx = _q.Index(i, _q.ny - 1);
                        T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety) + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 6)])/(1.0 + 3.0*_uy[idx]);
                        _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                        _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                        _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx])/36.0;
                    }
                }
            }
        }

        //  Function of setting boundary condition set q of AD for D3Q15 (diffusivity heterogenious)
        template<class T, template<class>class Q, class Fv, class Ff>
        void BoundaryConditionSetQ(Q<T>& _q, Fv _qnbc, const T *_ux, const T *_uy, const T *_uz, const T *_diffusivity, Ff _bctype) {
            //  On xmin
            if (_q.PEx == 0) {
                for (int j = 0; j < _q.ny; ++j) {
                    for (int k = 0; k < _q.nz; ++k) {
                        if (_bctype(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz)) {
                            int idx = _q.Index(0, j, k);
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(0 + _q.offsetx, j + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 4)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 13)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 1)] = tem0*(1.0 + 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc((_q.nx - 1) + _q.offsetx, j + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 1)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 12)])/(1.0 + 3.0*_ux[idx]);
                            _q.f[Q<T>::IndexF(idx, 4)] = tem0*(1.0 - 3.0*_ux[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, 0 + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 5)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 - 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 2)] = tem0*(1.0 + 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, (_q.ny - 1) + _q.offsety, k + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 2)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 + 3.0*_uy[idx]);
                            _q.f[Q<T>::IndexF(idx, 5)] = tem0*(1.0 - 3.0*_uy[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, j + _q.offsety, 0 + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 6)] + _q.f[Q<T>::IndexF(idx, 10)] + _q.f[Q<T>::IndexF(idx, 11)] + _q.f[Q<T>::IndexF(idx, 12)] + _q.f[Q<T>::IndexF(idx, 13)])/(1.0 - 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 3)] = tem0*(1.0 + 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 7)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 8)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 9)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 14)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] + 3.0*_uz[idx])/72.0;
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
                            T tem0 = 6.0*((1.0 + 1.0/(6.0*_diffusivity[idx]))*_qnbc(i + _q.offsetx, j + _q.offsety, (_q.nz - 1) + _q.offsetz) + _q.f[Q<T>::IndexF(idx, 3)] + _q.f[Q<T>::IndexF(idx, 7)] + _q.f[Q<T>::IndexF(idx, 8)] + _q.f[Q<T>::IndexF(idx, 9)] + _q.f[Q<T>::IndexF(idx, 14)])/(1.0 + 3.0*_uz[idx]);
                            _q.f[Q<T>::IndexF(idx, 6)] = tem0*(1.0 - 3.0*_uz[idx])/9.0;
                            _q.f[Q<T>::IndexF(idx, 10)] = tem0*(1.0 + 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 11)] = tem0*(1.0 - 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 12)] = tem0*(1.0 + 3.0*_ux[idx] - 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                            _q.f[Q<T>::IndexF(idx, 13)] = tem0*(1.0 - 3.0*_ux[idx] + 3.0*_uy[idx] - 3.0*_uz[idx])/72.0;
                        }
                    }
                }
            } 
        }
    }
}