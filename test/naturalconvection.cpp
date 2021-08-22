#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int nx = 101, ny = 51, nt = 100000, dt = 1000;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    D2Q9<double> pf(nx, ny);
    D2Q9<double> pg(nx, ny);
    double rho[pf.nxy], ux[pf.nxy], uy[pf.nxy], tem[pg.nxy], qx[pg.nxy], qy[pg.nxy];
    for (int idx = 0; idx < pf.nxy; ++ idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        tem[idx] = 0.5*(Th + Tl);   qx[idx] = 0.0;  qy[idx] = 0.0;
    }

    pf.SetBoundary([&](int _i, int _j) {    return (_j == 0 || _j == pf.ny - 1) ? 2 : 0;    });
    pg.SetBoundary([&](int _i, int _j) {    return 0;    });
    int boundaryt[pg.nbc];
    pg.SetBoundary(boundaryt, [&](int _i, int _j) {    return (_j == 0 || _j == pg.ny - 1) ? 1 : 0;    });
    double tembc[pg.nbc] = { 0 };
    for (int i = 0; i < pg.nx; ++i) {
        tembc[i + pg.offsetymin] = Th;
        tembc[i + pg.offsetymax] = Tl;
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        AD::Macro_Collide_Stream_NaturalConvection(
            pf, rho, ux, uy, nu,
            pg, tem, qx, qy, alpha, 
            0.0, 1.6e-5, 0.5*(Th + Tl), true
        );
        if (t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
            VTKExport file("result/thermal" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
            file.AddPointVector("q", 
                [&](int _i, int _j, int _k) { return qx[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }

        pf.Swap();
        pg.Swap();
        pf.BoundaryCondition();
        pg.BoundaryCondition();
        AD::BoundaryConditionSetT(pg, tembc, ux, uy, boundaryt);
        pf.SmoothCorner();
        pg.SmoothCorner();
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}