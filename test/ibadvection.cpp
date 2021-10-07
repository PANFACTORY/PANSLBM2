#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/ibm/advection_ibm.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int nx = 101, ny = 101, nt = 10000, dt = 100;
    double nu = 0.02, alpha = 0.02, radius = 0.3, v0 = 0.01, Th = 2.0, Tl = 1.0;
    double L = 2.0*radius*nx, surface = L*M_PI;
    D2Q9<double> pf(nx, ny), pg(nx, ny);
    double rho[pf.nxyz], ux[pf.nxyz], uy[pf.nxyz], tem[pg.nxyz], qx[pg.nxyz], qy[pg.nxyz];
    for (int idx = 0; idx < nx*ny; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0;
        tem[idx] = 0.0; qx[idx] = 0.0; qy[idx] = 0.0;
    }

    const int nb = (int)surface;
    double dv = surface/(double)nb;
    IBNS<double, D2Q9> bodyns(pf, nb, dv);
    for (int k = 0; k < nb; k++) {
        bodyns.SetBP(k, 0.5*L*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0, 0.5*L*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0);
        bodyns.SetBV(k, -v0*sin(2.0*M_PI*k/(double)nb), v0*cos(2.0*M_PI*k/(double)nb));
    }
    IBAD<double, D2Q9> bodyad(pg, nb, dv);
    for (int k = 0; k < nb; k++) {
        bodyad.SetBP(k, 0.5*L*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0, 0.5*L*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0);
        bodyad.SetBT(k, Th);
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze loop--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        AD::MacroCollideForceConvectionIB(
            pf, rho, ux, uy, nu, 
            pg, tem, qx, qy, alpha,
            bodyns, bodyad, true
        );
        if (t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
            VTKExport file("result/ibmad_" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("tem", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
            file.AddPointVector("q", 
                [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }
        pf.Stream();
        pg.Stream();
        pf.BoundaryCondition([](int _i, int _j) { return 0; });
        pg.BoundaryCondition([](int _i, int _j) { return 0; });
        pf.SmoothCorner();
        pg.SmoothCorner();

        bodyns.Update(pf, ux, uy);
        bodyad.Update(pg, tem);
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
}