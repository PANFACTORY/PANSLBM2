#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/elastic.h"
#include "../src/ibm/ibns.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int nx = 301, ny = 151, nt = 100000, dt = 1000;
    double nu = 0.1, u0 = 0.05;
    D2Q9<double> pf(nx, ny);
    double rho[nx*ny], ux[nx*ny], uy[nx*ny];
    for (int idx = 0; idx < nx*ny; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;  uy[idx] = 0.0;
    }
    pf.SetBoundary([&](int _i, int _j) { return (_j == 0 || _j == pf.ny - 1) ? 2 : 0;  });
    int boundaryup[pf.nbc];
    pf.SetBoundary(boundaryup, [&](int _i, int _j) {    return _i == 0 ? 1 : (_i == pf.nx - 1 ? 2 : 0); });
    double uxbc[pf.nbc], uybc[pf.nbc], rhobc[pf.nbc], usbc[pf.nbc];
    for (int j = 0; j < pf.ny; ++j) {
        uxbc[j + pf.offsetxmin] = u0;   uybc[j + pf.offsetxmin] = 0.0;  rhobc[j + pf.offsetxmin] = 1.0; usbc[j + pf.offsetxmin] = 0.0;
        uxbc[j + pf.offsetxmax] = 0.0;  uybc[j + pf.offsetxmax] = 0.0;  rhobc[j + pf.offsetxmax] = 1.0; usbc[j + pf.offsetxmax] = 0.0;
    }
    for (int i = 0; i < pf.nx; ++i) {
        uxbc[i + pf.offsetymin] = 0.0;  uybc[i + pf.offsetymin] = 0.0;  rhobc[i + pf.offsetymin] = 1.0; usbc[i + pf.offsetymin] = 0.0;
        uxbc[i + pf.offsetymax] = 0.0;  uybc[i + pf.offsetymax] = 0.0;  rhobc[i + pf.offsetymax] = 1.0; usbc[i + pf.offsetymax] = 0.0;
    }

    int mx = 11, my = 101;
    double rho00 = 1e5;
    D2Q9<double> pg(mx, my);
    double rho0[pg.nxy], vx[pg.nxy], vy[pg.nxy], px[pg.nxy], py[pg.nxy], sxx[pg.nxy], sxy[pg.nxy], syx[pg.nxy], syy[pg.nxy];
    for (int idx = 0; idx < pg.nxy; ++idx) {
        rho0[idx] = rho00;
        vx[idx] = 0.0;  vy[idx] = 0.0;
        sxx[idx] = 0.0; sxy[idx] = 0.0; syx[idx] = 0.0; syy[idx] = 0.0;
    }
    for (int i = 0; i < pg.nx; ++i) {
        for (int j = 0; j < pg.ny; ++j) {
            int idx = pg.Index(i, j);
            px[idx] = i + 0.3*(nx - 1); py[idx] = j + 1;
        }
    }
    pg.SetBoundary([&](int _i, int _j) {    return _j == 0 ? 1 : 0; });
    int boundarys[pg.nbc];
    pg.SetBoundary(boundarys, [&](int _i, int _j) { return _j == 0 ? 0: 1;  });
    double txbc[pg.nbc] = { 0 }, tybc[pg.nbc] = { 0 };

    IBNS<double, D2Q9<double> > body(pf, pg.nbc, 1.0);
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze loop--------------------
    EL::InitialCondition(pg, rho0, vx, vy, sxx, sxy, syx, syy);
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        //  Structure
        EL::Macro_Collide_Stream(pg, rho0, vx, vy, sxx, sxy, syx, syy, 0.8, true);
        pg.Swap();
        pg.BoundaryCondition();
        EL::BoundaryConditionSetStress(pg, txbc, tybc, boundarys);
        pg.SmoothCorner();
        for (int idx = 0; idx < pg.nxy; ++idx) {
            px[idx] += vx[idx]; py[idx] += vy[idx];
        }
        if (t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
            VTKExport file("result/fsi_solid_" + std::to_string(t/dt) + ".vtk", mx, my);
            file.AddPointVector("r", 
                [&](int _i, int _j, int _k) { return px[pg.Index(_i, _j)] - (_i + 0.3*(nx - 1)); },
                [&](int _i, int _j, int _k) { return py[pg.Index(_i, _j)] - (_j + 1); },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointVector("v", 
                [&](int _i, int _j, int _k) { return vx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return vy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointTensor("stress", 
                [&](int _i, int _j, int _k) { return sxx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return sxy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; },
                [&](int _i, int _j, int _k) { return syx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return syy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; },
                [](int _i, int _j, int _k) { return 0.0; },
                [](int _i, int _j, int _k) { return 0.0; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }

        //  Immersed Boundary Method
        for (int j = 0; j < pg.ny; ++j) {
            body.SetBP(j + pg.offsetxmin, px[pg.Index(0, j)], py[pg.Index(0, j)]);
            body.SetBV(j + pg.offsetxmin, vx[pg.Index(0, j)], vy[pg.Index(0, j)]);

            body.SetBP(j + pg.offsetxmax, px[pg.Index(pg.nx - 1, j)], py[pg.Index(pg.nx - 1, j)]);
            body.SetBV(j + pg.offsetxmax, vx[pg.Index(pg.nx - 1, j)], vy[pg.Index(pg.nx - 1, j)]);
        }
        for (int i = 0; i < pg.nx; ++i) {
            body.SetBP(i + pg.offsetymin, px[pg.Index(i, 0)], py[pg.Index(i, 0)]);
            body.SetBV(i + pg.offsetymin, vx[pg.Index(i, 0)], vy[pg.Index(i, 0)]);

            body.SetBP(i + pg.offsetymax, px[pg.Index(i, pg.ny - 1)], py[pg.Index(i, pg.ny - 1)]);
            body.SetBV(i + pg.offsetymax, vx[pg.Index(i, pg.ny - 1)], vy[pg.Index(i, pg.ny - 1)]);
        }
        body.Update(pf, ux, uy);
        for (int j = 0; j < pg.ny; ++j) {
            body.GetBG(j + pg.offsetxmin, txbc[j + pg.offsetxmin], tybc[j + pg.offsetxmin]);

            body.GetBG(j + pg.offsetxmax, txbc[j + pg.offsetxmax], tybc[j + pg.offsetxmax]);
        }
        for (int i = 0; i < pg.nx; ++i) {
            body.GetBG(i + pg.offsetymax, txbc[i + pg.offsetymax], tybc[i + pg.offsetymax]);
        }

        //  Fluid
        NS::Macro_Collide_Stream_IBM(pf, rho, ux, uy, nu, body, true);
        if (t%dt == 0) {
            VTKExport file("result/fsi_fluid_" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }
        pf.Swap();
        pf.BoundaryCondition();
        NS::BoundaryConditionSetU(pf, uxbc, uybc, boundaryup);
        NS::BoundaryConditionSetRho(pf, rhobc, usbc, boundaryup);
        pf.SmoothCorner();   
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
}