#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/equation/adjointelastic.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nx = 81, ny = 61, nt = 10000, dt = 1000;
    double rho0 = 1e6, stress0 = -1.0, p = 0.25, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(nx, ny);
    double rho[pf.nxy], ux[pf.nxy], uy[pf.nxy], rx[pf.nxy], ry[pf.nxy], sxx[pf.nxy], sxy[pf.nxy], syx[pf.nxy], syy[pf.nxy];
    double irho[pf.nxy], imx[pf.nxy], imy[pf.nxy], isxx[pf.nxy], isxy[pf.nxy], isyx[pf.nxy], isyy[pf.nxy];
    double s[pf.nxy], gamma[pf.nxy], sensitivity[pf.nxy];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow(i - 0.5*pf.nx, 2.0) + pow(j, 2.0) < pow(0.15*pf.nx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxy; ++idx) {
        rho[idx] = rho0;    
        ux[idx] = 0.0;  uy[idx] = 0.0;
        rx[idx] = 0.0;  ry[idx] = 0.0;
        sxx[idx] = 0.0; sxy[idx] = 0.0; syx[idx] = 0.0; syy[idx] = 0.0;
        irho[idx] = 0.0;
        imx[idx] = 0.0; imy[idx] = 0.0;
        isxx[idx] = 0.0;    isxy[idx] = 0.0;    isyx[idx] = 0.0;    isyy[idx] = 0.0;
        gamma[idx] = pow(s[idx], p);
    }

    pf.SetBoundary([&](int _i, int _j) {    return _i == 0 ? 1 : 0; });
    int boundarys[pf.nbc];
    pf.SetBoundary(boundarys, [&](int _i, int _j) { return _i == 0 ? 0: 1;  });
    double txbc[pf.nbc] = { 0 }, tybc[pf.nbc] = { 0 };
    for (int j = 0; j < pf.ny; ++j) {
        tybc[j + pf.offsetxmax] = fabs(j - 0.5*pf.ny < 0.08*pf.ny) ? stress0 : 0.0;  
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    
    //--------------------Direct analyze--------------------
    EL::InitialCondition(pf, rho, ux, uy, sxx, sxy, syx, syy);
    for (int t = 1; t <= nt; ++t) {
        EL::MacroExtended_Collide_Stream(pf, rho, ux, uy, sxx, sxy, syx, syy, 0.8, gamma, true);
        pf.Swap();
        pf.BoundaryCondition();
        EL::BoundaryConditionSetStress(pf, txbc, tybc, boundarys);
        pf.SmoothCorner();

        for (int idx = 0; idx < pf.nxy; ++idx) {
            rx[idx] += ux[idx]; ry[idx] += uy[idx];
        }
    }

    //--------------------Invert analyze--------------------
    AEL::InitialCondition(pf, irho, imx, imy, isxx, isxy, isyx, isyy, gamma);
    for (int t = 1; t <= nt; ++t) {
        AEL::Macro_Collide_Stream(pf, irho, imx, imy, isxx, isxy, isyx, isyy, 0.8, gamma, true);
        pf.Swap();
        pf.iBoundaryCondition();
        AEL::BoundaryConditionSetiStress(pf, txbc, tybc, rho, boundarys);
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int idx = 0; idx < pf.nxy; ++idx) {
        sensitivity[idx] = (4.5*(sxx[idx]*isxx[idx] + sxy[idx]*isxy[idx] + syx[idx]*isyx[idx] + syy[idx]*isyy[idx]) - 1.5*irho[idx]*(sxx[idx] + syy[idx]));
        if (sensitivitymax < fabs(sensitivity[idx])) {
            sensitivitymax = fabs(sensitivity[idx]);
        }
    }
    for (int idx = 0; idx < pf.nxy; ++idx) {  
        sensitivity[idx] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointelastic.vtk", nx, ny);
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
    file.AddPointTensor("s", 
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
    file.AddPointScaler("irho", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointTensor("is", 
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
    file.AddPointScaler("sensitivity", [&](int _i, int _j, int _k) { return sensitivity[pf.Index(_i, _j)]; });
    file.AddPointScaler("gamma", [&](int _i, int _j, int _k) { return gamma[pf.Index(_i, _j)]; });
    file.AddPointScaler("ss", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    return 0;
}