#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int nx = 201, ny = 21, nt = 100000, dt = 1000;
    double rho0 = 1e6, stress0 = -1.0;
    D2Q9<double> pf(nx, ny);
    double rho[pf.nxyz], ux[pf.nxyz], uy[pf.nxyz], rx[pf.nxyz], ry[pf.nxyz], sxx[pf.nxyz], sxy[pf.nxyz], syx[pf.nxyz], syy[pf.nxyz];
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = rho0;
        ux[idx] = 0.0;  uy[idx] = 0.0;
        rx[idx] = 0.0;  ry[idx] = 0.0;
        sxx[idx] = 0.0; sxy[idx] = 0.0; syx[idx] = 0.0; syy[idx] = 0.0;
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    EL::InitialCondition(pf, rho, ux, uy, sxx, sxy, syx, syy);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt != 0) {
            EL::Macro_Collide_Stream(pf, rho, ux, uy, sxx, sxy, syx, syy, 0.8);
        } else {
            EL::Macro_Collide_Stream(pf, rho, ux, uy, sxx, sxy, syx, syy, 0.8, true);

            std::cout << "t = " << t/dt << std::endl;
            VTKExport file("result/elastic" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointVector("r", 
                [&](int _i, int _j, int _k) { return rx[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return ry[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointTensor("stress", 
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
        }

        pf.Swap();
        pf.BoundaryCondition([=](int _i, int _j) { return _i == 0 ? 1 : 0; });
        EL::BoundaryConditionSetStress(pf, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == nx - 1 ? stress0 : 0.0; }, 
            [=](int _i, int _j) { return _i != 0; }
        );
        pf.SmoothCorner();

        for (int idx = 0; idx < pf.nxyz; ++idx) {
            rx[idx] += ux[idx]; ry[idx] += uy[idx];
        }
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}