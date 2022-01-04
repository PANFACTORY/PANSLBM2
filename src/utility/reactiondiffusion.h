#pragma once
#include <vector>

#include "src/LinearAlgebra/Models/Vector.h"
#include "src/FEM/Controller/ShapeFunction.h"
#include "src/FEM/Controller/GaussIntegration.h"
#include "src/LinearAlgebra/Solvers/CG.h"

using namespace PANSFEM2;

namespace PANSLBM2 {
    template<class T, template<class>class SF, template<class>class IC>
    Matrix<T> ElementalStiffnessMatrix(const std::vector<int>& _element, std::vector<Vector<T> >& _x, T _D, T _dt) {
        Matrix<T> Ke(_element.size(), _element.size());
		
        Matrix<T> X(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);   X(i, 1) = _x[_element[i]](1);
		}

        for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
            Ke += (_D*dNdX.Transpose()*dNdX + N*N.Transpose()/_dt)*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }

        return Ke;
    }

    template<class T, template<class>class SF, template<class>class IC>
    Vector<T> ElementalForceVector(const std::vector<int>& _element, std::vector<Vector<T> >& _x, std::vector<T>& _dFdphi, std::vector<T>& _phi, T _C, T _dt) {
        Vector<T> Fe(_element.size());
		
        Matrix<T> X(_element.size(), 2);
		for(int i = 0; i < _element.size(); i++){
			X(i, 0) = _x[_element[i]](0);   X(i, 1) = _x[_element[i]](1);
		}

        Vector<T> DFDPHI(_element.size()), PHI(_element.size());
        for(int i = 0; i < _element.size(); i++) {
            DFDPHI(i) = _dFdphi[_element[i]];
            PHI(i) = _phi[_element[i]];
        }

        for (int g = 0; g < IC<T>::N; g++) {
            Vector<T> N = SF<T>::N(IC<T>::Points[g]);
            Matrix<T> dNdr = SF<T>::dNdr(IC<T>::Points[g]);
			Matrix<T> dXdr = dNdr*X;
			T J = dXdr.Determinant();
			Matrix<T> dNdX = dXdr.Inverse()*dNdr;
            T dFdphi = N*DFDPHI, phi = N*PHI;
            Fe += N*(-_C*dFdphi + phi/_dt)*J*IC<T>::Weights[g][0]*IC<T>::Weights[g][1];
        }

        return Fe;
    }

    template<class T>
    void Assembling(
        LILCSR<T>& _K, std::vector<T>& _F, std::vector<T>& _u, Matrix<T>& _Ke, Vector<T>& _Fe, 
        const std::vector<int>& _nodetoglobal, const std::vector<int>& _element
    ) {
        for(int i = 0; i < _element.size(); i++) {
            if(_nodetoglobal[_element[i]] != -1) {
                for(int j = 0; j < _element.size(); j++) {
                    //----------Dirichlet condition NOT imposed----------
                    if(_nodetoglobal[_element[j]] != -1) {
                        _K.set(_nodetoglobal[_element[i]], _nodetoglobal[_element[j]], _K.get(_nodetoglobal[_element[i]], _nodetoglobal[_element[j]]) + _Ke(i, j));
                    }
                    //----------Dirichlet condition imposed----------
                    else {
                        _F[_nodetoglobal[_element[i]]] -= _Ke(i, j)*_u[_element[j]];
                    }
                }
                _F[_nodetoglobal[_element[i]]] += _Fe(i);
            }
        }
    }

    template<class T>
    void Disassembling(std::vector<T>& _u, const std::vector<T>& _result, const std::vector<int>& _nodetoglobal) {
        for(int i = 0; i < _nodetoglobal.size(); i++) {
            if(_nodetoglobal[i] != -1) {
                _u[i] = _result[_nodetoglobal[i]]; 
            }
        }
    }

    namespace ReactionDiffusion {
        template<class T, template<class>class P, class Ff, class Fv>
        void UpdateVariables(P<T>& _p, std::vector<T> &_phi, T _f, const std::vector<T> &_dfdphi, T _g, const std::vector<T> &_dgdphi, T _tau, T _dt, Ff _isphifixed, Fv _phivalue) {
            //  Get Lagrange multiplier for constraint function
            T lambda = T();
            if (_g >= T()) {
                T sumdfdphidgdphi = T(), sumdgdphidgdphi = T();
                for (int idx = 0; idx < _p.nxyz; ++idx) {
                    sumdfdphidgdphi += _dfdphi[idx]*_dgdphi[idx];
                    sumdgdphidgdphi += _dgdphi[idx]*_dgdphi[idx];
                }
                lambda = _g == T() ? -sumdfdphidgdphi/sumdgdphidgdphi : -2*sumdfdphidgdphi/sumdgdphidgdphi;
            }

            //  Get topological derivative
            std::vector<T> TDN(_p.nxyz);
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                TDN[idx] = _dfdphi[idx] + lambda*_dgdphi[idx];
            }

            //  Set up mesh where solve reaction diffusion equation
            std::vector<Vector<T> > nodes(_p.nxyz);
            for (int i = 0; i < _p.lx; i++){
                for (int j = 0; j < _p.ly; j++){
                    nodes[i + _p.lx*j] = { T(i), T(j) };
                }
            }
            std::vector<std::vector<int> > elements((_p.lx - 1)*(_p.ly - 1));
            for (int i = 0; i < _p.lx - 1; i++){
                for (int j = 0; j < _p.ly - 1; j++){
                    elements[i + (_p.lx - 1)*j] = { i + _p.lx*j, (i + 1) + _p.lx*j, (i + 1) + _p.lx*(j + 1), i + _p.lx*(j + 1) };
                }
            }
            std::vector<int> nodetoglobal(_p.nxyz, 0);
            int DEGREE = 0;
            for (int i = 0; i < _p.lx; ++i) {
                for (int j = 0; j < _p.ly; ++j) {
                    if (_isphifixed(i, j)) {
                        _phi[i + _p.lx*j] = _phivalue(i, j);
                        nodetoglobal[i + _p.lx*j] = -1;
                    } else {
                        nodetoglobal[i + _p.lx*j] = DEGREE++;
                    }
                }
            }

            //  Get parameter C
            T C = T();
            for(int idx = 0; idx < _p.nxyz; idx++) {
                C += fabs(_dfdphi[idx]);
            }
            C = elements.size()/C;
            
            //  Solve reaction diffusion equation
            LILCSR<T> A(DEGREE, DEGREE);
            std::vector<T> B(DEGREE, T());

            for (int i = 0; i < elements.size(); i++) {
                Matrix<T> Ke = ElementalStiffnessMatrix<T, ShapeFunction4Square, Gauss4Square>(elements[i], nodes, _tau*elements.size(), _dt);
                Vector<T> Fe = ElementalForceVector<T, ShapeFunction4Square, Gauss4Square>(elements[i], nodes, TDN, _phi, C, _dt);
                Assembling(A, B, _phi, Ke, Fe, nodetoglobal, elements[i]);
            }

            CSR<T> Amod(A);	
            Disassembling(_phi, ScalingCG(Amod, B, 100000, 1.0e-10), nodetoglobal);

            //  Modify level set function value
            for(int idx = 0; idx < _p.nxyz; idx++) {
                _phi[idx] = std::max(std::min(1.0, _phi[idx]), -1.0);
            }
        }
    }
}