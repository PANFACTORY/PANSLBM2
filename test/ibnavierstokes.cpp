#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/ibm/navierstokes_ibm.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int nx = 101, ny = 101, nt = 10000, dt = 100;
    double nu = 0.1, v0 = 0.01, radius = 0.3;
    double L = 2.0*radius*nx, surface = L*M_PI;
    D2Q9<double> pf(nx, ny);
    double rho[nx*ny], ux[nx*ny], uy[nx*ny];
    for (int idx = 0; idx < nx*ny; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;  uy[idx] = 0.0;
    }

    const int nb = (int)surface;
    double dv = surface/(double)nb;
    IBNS<double, D2Q9> body(pf, nb, dv);
    for (int k = 0; k < nb; k++) {
        body.SetBP(k, 0.5*L*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0, 0.5*L*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0);
        body.SetBV(k, -v0*sin(2.0*M_PI*k/(double)nb), v0*cos(2.0*M_PI*k/(double)nb));
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze loop--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        NS::MacroCollide(pf, rho, ux, uy, nu, true);
        if (t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
            VTKExport file("result/ibm_" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }
        pf.Stream();
        pf.BoundaryCondition([](int _i, int _j) { return 1; });
        pf.SmoothCorner();

        body.Update(pf, rho, ux, uy);
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
}