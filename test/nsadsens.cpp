//#define _USE_MPI_DEFINES
#define _USE_AVX_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/vtkxmlexport.h"
#include "../src/utility/normalize.h"

using namespace PANSLBM2;

int main() {
#ifdef _USE_MPI_DEFINES
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int mx = atoi(argv[1]), my = atoi(argv[2]);
    assert(mx*my == PeTot);
#else
    int MyRank = 0, mx = 1, my = 1;
#endif

    //--------------------Setting parameters--------------------
    int lx = 101, ly = 51, nt = 30000, dt = 100;
    double nu = 0.1, alp = 0.1/6.0, u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-8, epsdq = 1.0e-8;    
    double qa = 0.01, qb = 0.01, amax = 50.0, bmax = 0.1, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(lx, ly, MyRank, mx, my), pg(lx, ly, MyRank, mx, my);
    double *rho = new double[pg.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz];   
    double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *beta = new double[pg.nxyz], *dfds = new double[pf.nxyz], *dads = new double[pf.nxyz], *dbds = new double[pg.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow((i + pf.offsetx) - 0.5*pf.nx, 2.0) + pow((j + pf.offsety), 2.0) < pow(0.15*pf.nx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        tem[idx] = 1.0; qx[idx] = 0.0;  qy[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0;
        item[idx] = 0.0;    iqx[idx] = 0.0; iqy[idx] = 0.0;
        alpha[idx] = amax/(double)(pf.nx - 1)*qa*(1.0 - s[idx])/(s[idx] + qa);
        beta[idx] = bmax/(double)(pf.nx - 1)*qb*(1.0 - s[idx])/(s[idx] + qb);
        dads[idx] = -amax/(double)(pf.nx - 1)*qa*(qa + 1.0)/pow(qa + s[idx], 2.0);
        dbds[idx] = -bmax/(double)(pg.nx - 1)*qb*(qb + 1.0)/pow(qb + s[idx], 2.0);
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();                                    

    //--------------------Direct analyze--------------------
    if (MyRank == 0) {
        std::cout << "Direct analyse t = 0";
    }
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        AD::MacroBrinkmanCollideHeatExchange(
            pf, rho, ux, uy, alpha, nu,
            pg, tem, qx, qy, beta, alp, true
        );
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
        }
        pf.Stream();
        pg.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ly ? 1 : 0); });
        NS::BoundaryConditionSetU(pf, 
            [=](int _i, int _j) { return -u0*(_j - 0.33*ly)*(_j + 0.33*ly)/(0.33*ly*0.33*ly); }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ly; }
        );
        NS::BoundaryConditionSetRho(pf, 
            [=](int _i, int _j) { return 1.0; }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == lx - 1 && _j < 0.33*ly; }
        );
        pf.SmoothCorner();
        pg.BoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : 0; });
        AD::BoundaryConditionSetT(pg, 
            [=](int _i, int _j) { return tem0; }, 
            ux, uy, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ly; }
        );
        AD::BoundaryConditionSetQ(pg, 
            [=](int _i, int _j) { return q0; }, 
            ux, uy, alp, 
            [=](int _i, int _j) { return _i == lx - 1 || _j >= 0.33*ly; }
        );
        pg.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    if (MyRank == 0) {
        std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
    }
    ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
    AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
    for (int t = 1; t <= nt; ++t) {
        AAD::MacroBrinkmanCollideHeatExchange(
            pf, rho, ux, uy, irho, iux, iuy, imx, imy, alpha, nu,
            pg, tem, item, iqx, iqy, beta, alp, true
        );
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
        }
        pg.iStream();
        pf.iStream();
        pg.iBoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : 0; });
        AAD::iBoundaryConditionSetT(pg, ux, uy, [=](int _i, int _j) { return _i == 0 && _j < 0.33*ly; });
        AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return _i == lx - 1 || _j >= 0.33*ly; });
        pg.SmoothCorner();
        pf.iBoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ly ? 1 : 0); });
        ANS::iBoundaryConditionSetU(pf, 
            [=](int _i, int _j) { return -u0*(_j - 0.33*ly)*(_j + 0.33*ly)/(0.33*ly*0.33*ly); }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i == 0 && _j < 0.33*ly; }
        );
        AAD::iBoundaryConditionSetRho(pf, pg, rho, ux, uy, tem, [=](int _i, int _j) { return _i == lx - 1 && _j < 0.33*ly; });
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    AAD::SensitivityHeatExchange(pg, dfds, ux, uy, imx, imy, dads, tem, item, dbds);
    Normalize(dfds, pf.nxyz);

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << "\r" << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::string(20, ' ') << std::endl;

    //--------------------Export result--------------------
    VTKXMLExport file(pf, "result/nsadsens");
    file.AddPointData(pf, "p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pg, "T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
    file.AddPointData(pg, "q", 
        [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j)];  });
    file.AddPointData(pf, "ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "im", 
                [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
    file.AddPointData(pg, "iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j)]; });
    file.AddPointData(pg, "iq", 
        [&](int _i, int _j, int _k) { return iqx[pg.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iqy[pg.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointData(pg, "beta", [&](int _i, int _j, int _k) { return beta[pg.Index(_i, _j)]; });
    file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pg.Index(_i, _j)]; });
    
    delete[] rho;  delete[] ux;  delete[] uy;  delete[] tem; delete[] qx;  delete[] qy;
    delete[] irho; delete[] iux; delete[] iuy; delete[] imx; delete[] imy; delete[] item; delete[] iqx; delete[] iqy;
    delete[] s;    delete[] alpha; delete[] beta; delete[] dfds; delete[] dads; delete[] dbds;
 
    return 0;
}