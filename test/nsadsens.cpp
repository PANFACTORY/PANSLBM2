#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nx = 101, ny = 51, nt = 30000, tmax = nt - 1;
    double nu = 0.1, alp = 0.1/6.0, u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-8, epsdq = 1.0e-8;    
    double qa = 0.01, qb = 0.01, amax = 50.0, bmax = 0.1, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(nx, ny), pg(nx, ny);
    double *rho = new double[nx*ny], *ux = new double[nx*ny], *uy = new double[nx*ny];
    double *tem = new double[nx*ny], *qx = new double[nx*ny], *qy = new double[nx*ny];
    double *irho = new double[nx*ny], *iux = new double[nx*ny], *iuy = new double[nx*ny], *imx = new double[nx*ny], *imy = new double[nx*ny];
    double *item = new double[nx*ny], *iqx = new double[nx*ny], *iqy = new double[nx*ny];   
    double *s = new double[nx*ny], *alpha = new double[nx*ny], *beta = new double[nx*ny], *sensitivity = new double[nx*ny];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow(i - 0.5*pf.nx, 2.0) + pow(j, 2.0) < pow(0.15*pf.nx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        tem[idx] = 1.0; qx[idx] = 0.0;  qy[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0;
        item[idx] = 0.0;    iqx[idx] = 0.0; iqy[idx] = 0.0;
        alpha[idx] = amax/(double)(pf.nx - 1)*qa*(1.0 - s[idx])/(s[idx] + qa);
        beta[idx] = bmax/(double)(pf.nx - 1)*qb*(1.0 - s[idx])/(s[idx] + qb);
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();                                    

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        AD::MacroBrinkmanCollideStreamHeatExchange(
            pf, rho, ux, uy, alpha, nu,
            pg, tem, qx, qy, beta, alp, true
        );
        pf.Stream();
        pg.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ny ? 1 : 0); });
        NS::BoundaryConditionSetU(pf, 
            [=](int _i, int _j) { return -u0*(_j - 0.33*ny)*(_j + 0.33*ny)/(0.33*ny*0.33*ny); }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ny; }
        );
        NS::BoundaryConditionSetRho(pf, 
            [=](int _i, int _j) { return 1.0; }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == nx - 1 && _j < 0.33*ny; }
        );
        pf.SmoothCorner();
        pg.BoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : 0; });
        AD::BoundaryConditionSetT(pg, 
            [=](int _i, int _j) { return tem0; }, 
            ux, uy, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ny; }
        );
        AD::BoundaryConditionSetQ(pg, 
            [=](int _i, int _j) { return q0; }, 
            ux, uy, alp, 
            [=](int _i, int _j) { return _i == nx - 1 || _j >= 0.33*ny; }
        );
        pg.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
    AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
    for (int t = 1; t <= nt; ++t) {
        AAD::MacroBrinkmanCollideStreamHeatExchange(
            pf, rho, ux, uy, irho, iux, iuy, imx, imy, alpha, nu,
            pg, tem, item, iqx, iqy, beta, alp, true
        );
        pg.iStream();
        pf.iStream();
        pg.iBoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : 0; });
        AAD::iBoundaryConditionSetT(pg, ux, uy, [=](int _i, int _j) { return _i == 0 && _j < 0.33*ny; });
        AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return _i == nx - 1 || _j >= 0.33*ny; });
        pg.SmoothCorner();
        pf.iBoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ny ? 1 : 0); });
        ANS::iBoundaryConditionSetU(pf, 
            [=](int _i, int _j) { return -u0*(_j - 0.33*ny)*(_j + 0.33*ny)/(0.33*ny*0.33*ny); }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ny; }
        );
        AAD::iBoundaryConditionSetRho(pf, pg, rho, ux, uy, tem, [=](int _i, int _j) { return _i == nx - 1 && _j < 0.33*ny; });
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        sensitivity[idx] = 3.0*(imx[idx]*ux[idx] + imy[idx]*uy[idx])*(-amax/(double)(pf.nx - 1)*qa*(qa + 1.0)/pow(qa + s[idx], 2.0)) - (1.0 - tem[idx])*(1.0 + item[idx])*(-bmax/(double)(pg.nx - 1)*qb*(qb + 1.0)/pow(qb + s[idx], 2.0));
        if (sensitivitymax < fabs(sensitivity[idx])) {
            sensitivitymax = fabs(sensitivity[idx]);
        }
    }
    for (int idx = 0; idx < nx*ny; ++idx) {  
        sensitivity[idx] /= sensitivitymax;
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    //--------------------Export result--------------------
    VTKExport file("result/nsadsens.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
    file.AddPointVector("q", 
        [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[pf.Index(_i, _j)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointVector("im", 
                [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
    file.AddPointScaler("iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j)]; });
    file.AddPointVector("iq", 
        [&](int _i, int _j, int _k) { return iqx[pg.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iqy[pg.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointScaler("beta", [&](int _i, int _j, int _k) { return beta[pg.Index(_i, _j)]; });
    file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pg.Index(_i, _j)]; });
    
    delete[] rho, ux, uy, tem, qx, qy, irho, iux, iuy, imx, imy, item, iqx, iqy;
    delete[] s, alpha, beta, sensitivity;
 
    return 0;
}