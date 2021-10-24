#define _USE_MATH_DEFINES
//#define _USE_MPI_DEFINES
//#define _USE_AVX_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/vtkxmlexport.h"
#include "../src\utility/normalize.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
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

    //--------------------Set parameters--------------------
    int lx = 51, ly = 101, nt0 = 200000, dt = 500, nt = 50000, period = 5000;
    double viscosity = 0.1, diff_fluid = viscosity/6.0, Th = 1.0, Tl = 0.0, gx = 0.0, gy = 1e-4;
    double alphamax = 1e4, diff_solid = diff_fluid*10.0, qf = 1e-2, qg = 1e0, s_min = 0.1, s_max = 0.9;
    D2Q9<double> pf(lx, ly, MyRank, mx, my), pg(lx, ly, MyRank, mx, my);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pf.nxyz], *qy = new double[pf.nxyz], *qxp = new double[pf.nxyz], *qyp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqxp = new double[pg.nxyz], *iqyp = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0;  ux[idx] = 0.0;  uy[idx] = 0.0;  uxp[idx] = 0.0;  uyp[idx] = 0.0;  tem[idx] = Tl;  qx[idx] = 0.0;  qy[idx] = 0.0;   qxp[idx] = 0.0; qyp[idx] = 0.0;
        irho[idx] = 0.0; iux[idx] = 0.0; iuy[idx] = 0.0; iuxp[idx] = 0.0; iuyp[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0; iqxp[idx] = 0.0; iqyp[idx] = 0.0;
    }
    double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dads = new double[pf.nxyz], *dkds = new double[pf.nxyz], *dfds = new double[pf.nxyz];
    double *gi = new double[pg.nxyz*pg.nc], *igi = new double[pg.nxyz*pg.nc];
    double *directionx = new double[pf.nxyz], *directiony = new double[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = (j + pf.offsety) < lx/2 ? (((i + pf.offsetx) < 0.5*lx && (j + pf.offsety) < 0.25*ly) ? s_min : s_max) : 1.0;
            directionx[idx] = ((i + pf.offsetx) == lx/2 && (j + pf.offsety) > 9*ly/10) ? 1.0 : 0.0;
            directiony[idx] = 0.0;
        }
    }
    for (int idx = 0; idx < pf.nxyz; idx++) {
        diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*s[idx]*(1.0 + qg)/(s[idx] + qg);
        alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - s[idx])/(s[idx] + qf);
        dkds[idx] = (diff_fluid - diff_solid)*qg*(1.0 + qg)/pow(s[idx] + qg, 2.0);
        dads[idx] = -alphamax/(double)(ly - 1)*qf*(1.0 + qf)/pow(s[idx] + qf, 2.0);
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    if (MyRank == 0) {
        std::cout << "Direct analyse t = 0";
    }
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt0 + nt; ++t) {
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
        }
        AD::MacroBrinkmanCollideNaturalConvection(
            pf, rho, ux, uy, alpha, viscosity, 
            pg, tem, qx, qy, diffusivity, gx, gy, Tl, true, gi
        );

        pf.Stream();
        pg.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return 1; });
        pg.BoundaryCondition([=](int _i, int _j) { return 0; });
        AD::BoundaryConditionSetT(pg, 
            [=](int _i, int _j) { return _i == 0 ? Th*(1 - cos(2*M_PI*t/(double)period)) : Tl; }, 
            ux, uy, 
            [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); }
        );
        AD::BoundaryConditionSetQ(pg, 
            [=](int _i, int _j) { return 0.0; },
            ux, uy, alpha, 
            [=](int _i, int _j) { return (_i == 0 && ly/2 <= _j) || (_i == lx - 1 && ly/2 <= _j) || _j == 0 || _j == ly - 1; }
        );
        pf.SmoothCorner();
        pg.SmoothCorner();

        pf.BoundaryConditionAlongXEdge(lx/5, 1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        pf.BoundaryConditionAlongYEdge(ly/2, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pf.BoundaryConditionAlongXEdge(4*lx/5, -1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        pf.BoundaryConditionAlongYEdge(9*ly/10, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pf.SmoothCornerAt(lx/5, ly/2, -1, -1);
        pf.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
        pf.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
        pf.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
        AD::BoundaryConditionSetQAlongXEdge(pg, lx/5, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, alpha, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        AD::BoundaryConditionSetQAlongYEdge(pg, ly/2, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, alpha, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        AD::BoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, [=](int _i, int _j) { return 0.0; }, ux, uy, alpha, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        AD::BoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, alpha, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
        pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
        pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
        pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
    }

    //--------------------Inverse analyze--------------------
    if (MyRank == 0) {
        std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
    }
    ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
    AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
    for (int t = 1; t <= nt0 + nt; ++t) {
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
        }
        AAD::MacroBrinkmanCollideNaturalConvectionMassFlow(
            pf, rho, ux, uy, irho, iux, iuy, imx, imy, alpha, viscosity,
            pg, tem, item, iqx, iqy, diffusivity, gx, gy, 
            directionx, directiony, true, igi
        );

        pf.iStream();
        pg.iStream();
        pf.iBoundaryCondition([=](int _i, int _j) { return 1; });
        pg.iBoundaryCondition([=](int _i, int _j) { return 0; });
        AAD::iBoundaryConditionSetT(pg, ux, uy, [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); });
        AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return (_i == 0 && ly/2 <= _j) || (_i == lx - 1 && ly/2 <= _j) || _j == 0 || _j == ly - 1; });
        pf.SmoothCorner();
        pg.SmoothCorner();

        pf.iBoundaryConditionAlongXEdge(lx/5, 1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        pf.iBoundaryConditionAlongYEdge(ly/2, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pf.iBoundaryConditionAlongXEdge(4*lx/5, -1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        pf.iBoundaryConditionAlongYEdge(9*ly/10, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pf.SmoothCornerAt(lx/5, ly/2, -1, -1);
        pf.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
        pf.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
        pf.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
        AAD::iBoundaryConditionSetQAlongXEdge(pg, lx/5, 1, ux, uy, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        AAD::iBoundaryConditionSetQAlongYEdge(pg, ly/2, 1, ux, uy, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        AAD::iBoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, ux, uy, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
        AAD::iBoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, ux, uy, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
        pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
        pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
        pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
        pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
    }

    //--------------------Get sensitivity--------------------
    AAD::SensitivityBrinkmanDiffusivity(pg, dfds, ux, uy, imx, imy, dads, tem, item, iqx, iqy, gi, igi, diffusivity, dkds);
    Normalize(dfds, pg.nxyz);

    //--------------------Export result--------------------
    VTKXMLExport file(pf, "result/ncpump");
    file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
    file.AddPointData(pf, "q", 
        [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) {   return irho[pf.Index(_i, _j)];  });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) {   return iux[pf.Index(_i, _j)];   },
        [&](int _i, int _j, int _k) {   return iuy[pf.Index(_i, _j)];   },
        [](int _i, int _j, int _k) {   return 0.0;   }
    );
    file.AddPointData(pf, "im", 
        [&](int _i, int _j, int _k) {   return imx[pf.Index(_i, _j)];   },
        [&](int _i, int _j, int _k) {   return imy[pf.Index(_i, _j)];   },
        [](int _i, int _j, int _k) {   return 0.0;   }
    );
    file.AddPointData(pf, "iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j)];  });
    file.AddPointData(pf, "iq", 
        [&](int _i, int _j, int _k) {   return iqx[pg.Index(_i, _j)];   },
        [&](int _i, int _j, int _k) {   return iqy[pg.Index(_i, _j)];   },
        [](int _i, int _j, int _k) {   return 0.0;   }
    );
    file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)];    });
    file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) { return dfds[pf.Index(_i, _j)];    });

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }

    delete[] rho;   delete[] ux;    delete[] uy;    delete[] uxp;   delete[] uyp;
    delete[] tem;   delete[] qx;    delete[] qy;    delete[] qxp;   delete[] qyp;
    delete[] s; delete[] diffusivity;   delete[] alpha; delete[] dkds;  delete[] dads;  delete[] dfds;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] iuxp;  delete[] iuyp;  delete[] imx;   delete[] imy;
    delete[] item;  delete[] iqx;   delete[] iqy;   delete[] iqxp;  delete[] iqyp;
    delete[] gi;    delete[] igi;   delete[] directionx;    delete[] directiony;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}