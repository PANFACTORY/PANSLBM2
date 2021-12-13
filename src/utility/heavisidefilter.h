#pragma once
#include <cmath>
#include <vector>
#include <cassert>

namespace PANSLBM2 {
    namespace HeavisideFilter {
        template<class T, template<class>class P, class F>
        std::vector<T> GetFilteredVariable(P<T>& _p, T _R, T _beta, const std::vector<T> &_s, F _weight) {
            assert(_R > T());
            int nR = (int)_R;

            //  Filter value
            std::vector<T> rho(_p.nxyz, T());
            for(int i1 = 0; i1 < _p.nx; ++i1){
                for(int j1 = 0; j1 < _p.ny; ++j1){
                    for(int k1 = 0; k1 < _p.nz; ++k1){
                        int idx = _p.Index(i1, j1, k1);
                        T wssum = T(), wsum = T();
                        for (int i2 = i1 - nR, i2max = i1 + nR; i2 <= i2max; ++i2) {
                            for (int j2 = j1 - nR, j2max = j1 + nR; j2 <= j2max; ++j2) {
                                for (int k2 = k1 - nR, k2max = k1 + nR; k2 <= k2max; ++k2) {                     
                                    double distance = sqrt(pow(i1 - i2, 2.0) + pow(j1 - j2, 2.0) + pow(k1 - k2, 2.0));
                                    if (distance <= _R) {
                                        T weight = _weight(i1 + _p.offsetx, j1 + _p.offsety, k1 + _p.offsetz, i2 + _p.offsetx, j2 + _p.offsety, k2 + _p.offsetz);
                                        if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                                                //  In self PE
                                            wssum += weight*_s[_p.Index(i2, j2, k2)]
                                            wsum += weight;
                                        } 
                                    }
                                }
                            }
                        }
                        rho[idx] = 0.5*(tanh(0.5*_beta) + tanh(_beta*(wssum/wsum - 0.5)))/tanh(_beta);
                    }
                }
            }

            return rho;
        }
    
        template<class T, template<class>class P>
        std::vector<T> GetFilteredVariable(P<T>& _p, T _R, T _beta, const std::vector<T> &_s) {
            return GetFilteredVariable(_p, _R, _beta, _s, [=](int _i1, int _j1, int _k1, int _i2, int _j2, int _k2){ return (_R - sqrt(pow(_i1 - _i2, 2.0) + pow(_j1 - _j2, 2.0) + pow(_k1 - _k2, 2.0)))/_R; });
        }

        template<class T, template<class>class P, class F>
        std::vector<T> GetFilteredSensitivity(P<T>& _p, T _R, T _beta, const std::vector<T> &_s, const std::vector<T> &_dfdrho, F _weight) {
            assert(_R > T());
            int nR = (int)_R;

            //  Filter value
            std::vector<T> drhodstilde(_p.nxyz, T());
            for(int i1 = 0; i1 < _p.nx; ++i1){
                for(int j1 = 0; j1 < _p.ny; ++j1){
                    for(int k1 = 0; k1 < _p.nz; ++k1){
                        int idx = _p.Index(i1, j1, k1);
                        T wssum = T(), wsum = T();
                        for (int i2 = i1 - nR, i2max = i1 + nR; i2 <= i2max; ++i2) {
                            for (int j2 = j1 - nR, j2max = j1 + nR; j2 <= j2max; ++j2) {
                                for (int k2 = k1 - nR, k2max = k1 + nR; k2 <= k2max; ++k2) {                     
                                    double distance = sqrt(pow(i1 - i2, 2.0) + pow(j1 - j2, 2.0) + pow(k1 - k2, 2.0));
                                    if (distance <= _R) {
                                        T weight = _weight(i1 + _p.offsetx, j1 + _p.offsety, k1 + _p.offsetz, i2 + _p.offsetx, j2 + _p.offsety, k2 + _p.offsetz);
                                        if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                                                //  In self PE
                                            wssum += weight*_v[_p.Index(i2, j2, k2)]
                                            wsum += weight;
                                        } 
                                    }
                                }
                            }
                        }
                        drhodstilde[idx] = 0.5*_beta*(1.0 - pow(tanh(_beta*(wssum/wsum - 0.5)), 2.0))/tanh(0.5*_beta);
                    }
                }
            }

            //  Filter value
            std::vector<T> dfds(_p.nxyz, T());
            for(int i1 = 0; i1 < _p.nx; ++i1){
                for(int j1 = 0; j1 < _p.ny; ++j1){
                    for(int k1 = 0; k1 < _p.nz; ++k1){
                        int idx = _p.Index(i1, j1, k1);
                        T wsum = T();
                        for (int i2 = i1 - nR, i2max = i1 + nR; i2 <= i2max; ++i2) {
                            for (int j2 = j1 - nR, j2max = j1 + nR; j2 <= j2max; ++j2) {
                                for (int k2 = k1 - nR, k2max = k1 + nR; k2 <= k2max; ++k2) {                     
                                    double distance = sqrt(pow(i1 - i2, 2.0) + pow(j1 - j2, 2.0) + pow(k1 - k2, 2.0));
                                    if (distance <= _R) {
                                        T weight = _weight(i1 + _p.offsetx, j1 + _p.offsety, k1 + _p.offsetz, i2 + _p.offsetx, j2 + _p.offsety, k2 + _p.offsetz);
                                        if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                                                //  In self PE
                                            dfds[idx] += weight*_dfdrho[_p.Index(i2, j2, k2)]*drhodstilde[_p.Index(i2, j2, k2)];
                                            wsum += weight;
                                        } 
                                    }
                                }
                            }
                        }
                        dfds[idx] /= wsum;
                    }
                }
            }

            return dfds;
        }

        template<class T, template<class>class P>
        std::vector<T> GetFilteredSensitivity(P<T>& _p, T _R, T _beta, const std::vector<T> &_s, const std::vector<T> &_dfdrho) {
            return GetFilteredSensitivity(_p, _R, _beta, _s, _dfdrho, [=](int _i1, int _j1, int _k1, int _i2, int _j2, int _k2){ return (_R - sqrt(pow(_i1 - _i2, 2.0) + pow(_j1 - _j2, 2.0) + pow(_k1 - _k2, 2.0)))/_R; });
        }
    }
}