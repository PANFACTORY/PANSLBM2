#define _USE_AVX_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/equation/adjointelastic.h"
#include "../src/utility/reactiondiffusion.h"
#include "../src/utility/residual.h"
#include "../src/utility/vtkxmlexport.h"
#include "../src/utility/normalize.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    //********************Set parameters********************
    int lx = 81, ly = 61, nt = 10000, dt = 1000, nk = 50;
    double rho0 = 1e6, stress0 = -1.0, weightlimit = 0.5, dgsi = 0.1, tau = 1e-5;
    D2Q9<double> pf(lx, ly);
    double rho[pf.nxyz], ux[pf.nxyz], uy[pf.nxyz], rx[pf.nxyz], ry[pf.nxyz], sxx[pf.nxyz], sxy[pf.nxyz], syx[pf.nxyz], syy[pf.nxyz];
    double irho[pf.nxyz], imx[pf.nxyz], imy[pf.nxyz], isxx[pf.nxyz], isxy[pf.nxyz], isyx[pf.nxyz], isyy[pf.nxyz];
    double chi[pf.nxyz];

    std::vector<double> s(pf.nxyz, 1.0);
	
    for (int k = 0; k < nk; k++) {
        //********************Set parameter********************
        for (int idx = 0; idx < s.size(); idx++) {
            chi[idx] = s[idx] >= 0.0 ? 1.0 : 0.0;
        }

        //********************Get weight********************
        double g = -1.0;
        std::vector<double> dgds(s.size(), 0.0);
        for(int idx = 0; idx < pf.nxyz; idx++){
            g += chi[idx]/(weightlimit*(double)(lx*ly));
            dgds[idx] = 1.0/(weightlimit*(double)(lx*ly)); 
        }

        //********************Direct Analyse********************
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            rho[idx] = rho0;    
            ux[idx] = 0.0;  uy[idx] = 0.0;
            rx[idx] = 0.0;  ry[idx] = 0.0;
            sxx[idx] = 0.0; sxy[idx] = 0.0; syx[idx] = 0.0; syy[idx] = 0.0;
        }
        EL::InitialCondition(pf, rho, ux, uy, sxx, sxy, syx, syy);
        for (int t = 1; t <= nt; ++t) {
            EL::MacroExtendedCollide(pf, rho, ux, uy, sxx, sxy, syx, syy, 0.8, chi, true);
            pf.Stream();
            pf.BoundaryCondition([=](int _i, int _j) { return _i == 0 ? 1 : 0; });
            EL::BoundaryConditionSetStress(pf, 
                [=](int _i, int _j) { return 0.0; }, 
                [=](int _i, int _j) { return (_i == lx - 1 && fabs(_j - 0.5*ly) < 0.08*ly) ? stress0 : 0.0; }, 
                [=](int _i, int _j) { return _i != 0; }
            );
            pf.SmoothCorner();

            for (int idx = 0; idx < pf.nxyz; ++idx) {
                rx[idx] += ux[idx]; ry[idx] += uy[idx];
            }
        }

        //********************Invert analyze********************
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            irho[idx] = 0.0;
            imx[idx] = 0.0; imy[idx] = 0.0;
            isxx[idx] = 0.0;    isxy[idx] = 0.0;    isyx[idx] = 0.0;    isyy[idx] = 0.0;
        }
        AEL::InitialCondition(pf, irho, imx, imy, isxx, isxy, isyx, isyy, chi);
        for (int t = 1; t <= nt; ++t) {
            AEL::MacroCollide(pf, irho, imx, imy, isxx, isxy, isyx, isyy, 0.8, chi, true);
            pf.iStream();
            pf.iBoundaryCondition([=](int _i, int _j) { return _i == 0 ? 1 : 0; });
            AEL::iBoundaryConditionSetStress(pf, 
                [=](int _i, int _j) { return 0.0; }, 
                [=](int _i, int _j) { return (_i == lx - 1 && fabs(_j - 0.5*ly) < 0.08*ly) ? stress0 : 0.0; }, 
                rho, 
                [=](int _i, int _j) { return _i != 0; }
            );
            pf.SmoothCorner();
        }

        //********************Get sensitivity********************
        double f = 0.0;
        for (int j = 0; j < ly; ++j) {
            if (fabs(j - 0.5*ly) < 0.08*ly) {
                f += stress0*ry[pf.Index(lx - 1, j)];
            }
        }
        std::vector<double> dfds(s.size(), 0.0);
        AEL::SensitivityCompliance(pf, dfds.data(), sxx, sxy, syx, syy, irho, isxx, isxy, isyx, isyy);
        Normalize(dfds.data(), pf.nxyz);
        std::cout << "\r" << std::fixed << std::setprecision(6) << k << " " << f << " " << g << std::endl;

        //********************Check convergence********************
        /*if(g < 0.0 && optimizer.IsConvergence(f)){
            if (MyRank == 0) {
                std::cout << std::endl << "-----Optimized-----" << std::endl;
            }
            break;
        }*/

        //********************Update variable********************
        ReactionDiffusion::UpdateVariables(pf, s, f, dfds, g, dgds, tau, dgsi, [](int _i, int _j){ return false; }, [](int _i, int _j){ return 0.0; });
    }

    //********************Export result********************
    VTKXMLExport file(pf, "result/cantilberLSM");
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "r", 
        [&](int _i, int _j, int _k) { return rx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return ry[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "s", 
        [&](int _i, int _j, int _k) { return sxx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return sxy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [&](int _i, int _j, int _k) { return syx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return syy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "is", 
        [&](int _i, int _j, int _k) { return isxx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return isxy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [&](int _i, int _j, int _k) { return isyx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return isyy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "phi", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "chi", [&](int _i, int _j, int _k) { return chi[pf.Index(_i, _j)]; });

    return 0;
}