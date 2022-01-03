#pragma once
#include <iostream>
#include <vector>
#include <numeric>

#include "src/LinearAlgebra/Models/Vector.h"
#include "src/PrePost/Mesher/SquareMesh.h"
#include "src/FEM/Equation/ReactionDiffusion.h"
#include "src/FEM/Equation/General.h"
#include "src/FEM/Controller/ShapeFunction.h"
#include "src/FEM/Controller/GaussIntegration.h"
#include "src/FEM/Controller/BoundaryCondition.h"
#include "src/FEM/Controller/Assembling.h"
#include "src/LinearAlgebra/Solvers/CG.h"

using namespace PANSFEM2;

namespace PANSLBM2 {
    namespace ReactionDiffusion {
        template<class T, template<class>class P>
        void UpdateVariables(P<T>& _p, std::vector<T> &_s, T _f, const std::vector<T> &_dfds, T _g, const std::vector<T> &_dgds, T _tau) {
            T dt = 0.1;

            //  Get Lagrange multiplier for constraint function
            T lambda = T();
            if (_g >= T()) {
                T sumdfdsdgds = T(), sumdgdsdgds = T();
                for (int idx = 0; idx < _p.nxyz; ++idx) {
                    sumdfdsdgds += _dfds[idx]*_dgds[idx];
                    sumdgdsdgds += _dgds[idx]*_dgds[idx];
                }
                lambda = _g == T() ? -sumdfdsdgds/sumdgdsdgds : -2*sumdfdsdgds/sumdgdsdgds;
            }

            //  Set up reaction diffusion equation
            std::vector<Vector<T> > phi(_p.nxyz, Vector<T>(1));
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                phi[idx](0) = _s[idx];
            }

            std::vector<Vector<T> > nodes(_p.nxyz);
            for(int i = 0; i < _p.lx; i++){
                for(int j = 0; j < _p.ly; j++){
                    nodes[i + _p.lx*j] = { T(i), T(j) };
                }
            }
            std::vector<std::vector<int> > elements((_p.lx - 1)*(_p.ly - 1));
            for(int i = 0; i < _p.lx - 1; i++){
                for(int j = 0; j < _p.ly - 1; j++){
                    elements[i + (_p.lx - 1)*j] = { i + _p.lx*j, (i + 1) + _p.lx*j, (i + 1) + _p.lx*(j + 1), i + _p.lx*(j + 1) };
                }
            }
            std::vector<std::pair<std::pair<int, int>, T> > phifixed;
            for(int i = 0; i < _p.lx; i++){
                phifixed.push_back({ { i + _p.lx*0, 0 }, 1.0 });
                phifixed.push_back({ { i + _p.lx*(_p.ly - 1), 0 }, 1.0 });
            }
            for(int j = 0; j < _p.ly; j++){
                phifixed.push_back({ { 0 + _p.lx*j, 0 }, 1.0 });
                phifixed.push_back({ { (_p.lx - 1) + _p.lx*j, 0 }, 1.0 });
            }

            std::vector<Vector<T> > TDN(_p.nxyz, Vector<T>(1));
            for (int idx = 0; idx < _p.nxyz; ++idx) {
                TDN[idx](0) = _dfds[idx];
            }
            
            T C = T();
            for(int idx = 0; idx < _p.nxyz; idx++) {
                C += fabs(TDN[idx](0));
            }
            C = elements.size()/C;

            for (int idx = 0; idx < _p.nxyz; ++idx) {
                TDN[idx](0) = _dfds[idx] + lambda*_dgds[idx];
            }

            std::vector<std::vector<int> > nodetoglobal(_p.nxyz, std::vector<int>(1, 0));
            SetDirichlet(phi, nodetoglobal, phifixed);
            int DEGREE = Renumbering(nodetoglobal);

            LILCSR<T> A(DEGREE, DEGREE);
            std::vector<T> B(DEGREE, T());

            for (int i = 0; i < elements.size(); i++) {
                std::vector<std::vector<std::pair<int, int> > > nodetoelement;
                Matrix<T> Me, Ke;
                Vector<T> Fe;
                ReactionDiffusionConsistentMass<T, ShapeFunction4Square, Gauss4Square>(Me, nodetoelement, elements[i], { 0 }, nodes);
                ReactionDiffusionStiffness<T, ShapeFunction4Square, Gauss4Square>(Ke, nodetoelement, elements[i], { 0 }, nodes, _tau*elements.size());
                ReactionDiffusionReaction<T, ShapeFunction4Square, Gauss4Square>(Fe, nodetoelement, elements[i], { 0 }, nodes, TDN, [&](T _u, Vector<T> _dudX) {
                    return -C*_u;
                });
                Vector<T> phie = ElementVector(phi, nodetoelement, elements[i]);
                Matrix<T> Ae = Me/dt + Ke;
                Vector<T> Be = Me/dt*phie + Fe; 
                Assembling(A, B, phi, Ae, nodetoglobal, nodetoelement, elements[i]);
                Assembling(B, Be, nodetoglobal, nodetoelement, elements[i]);
            }

            CSR<T> Amod(A);	
            std::vector<T> result = ScalingCG(Amod, B, 100000, 1.0e-10);
            Disassembling(phi, result, nodetoglobal);

            for(int idx = 0; idx < _p.nxyz; idx++) {
                _s[idx] = std::max(std::min(1.0, phi[idx](0)), -1.0);
            }
        }
    }
}