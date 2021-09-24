#define _USE_MPI_DEFINES
#include <iostream>
#include <cmath>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/utility/vtkxmlexport.h"

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
    int lx = 101, ly = 51, nt = 10000, dt = 100;
    double nu = 0.1, u0 = 0.0109, rho0 = 1.0, epsdu = 1.0e-4, q = 1e-2, amax = 50.0, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz];
    double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *dfds = new double[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow((i + pf.offsetx) - 0.5*pf.lx, 2.0) + pow((j + pf.offsety), 2.0) < pow(0.15*pf.lx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0;
        alpha[idx] = amax/(double)(pf.lx - 1)*q*(1.0 - s[idx])/(s[idx] + q);
    }

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::MacroBrinkmanCollideStream(pf, rho, ux, uy, nu, alpha, true);
        pf.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ly ? 1 : 0); });
        NS::BoundaryConditionSetU(pf, 
            [&](int _i, int _j) { return -u0*(_j - 0.33*ly)*(_j + 0.33*ly)/(0.33*ly*0.33*ly); }, 
            [&](int _i, int _j) { return 0.0; }, 
            [&](int _i, int _j) { return _i == 0 && _j < 0.33*ly; }
        );
        NS::BoundaryConditionSetRho(pf, 
            [&](int _i, int _j) { return 1.0; }, 
            [&](int _i, int _j) { return 0.0; }, 
            [&](int _i, int _j) { return _i == lx - 1 && _j < 0.33*ly; }
        );
        pf.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        ANS::MacroBrinkmanCollideStream(pf, rho, ux, uy, irho, iux, iuy, imx, imy, nu, alpha, true);
        pf.iStream();
        pf.iBoundaryCondition([=](int _i, int _j) { return _j == 0 ? 2 : (_j >= 0.33*ly ? 1 : 0); });
        ANS::iBoundaryConditionSetU(pf, 
            [&](int _i, int _j) { return -u0*(_j - 0.33*ly)*(_j + 0.33*ly)/(0.33*ly*0.33*ly); }, 
            [&](int _i, int _j) { return 0.0; }, 
            [&](int _i, int _j) { return _i == 0 && _j < 0.33*ly; },
            1.0
        );
        ANS::iBoundaryConditionSetRho2D(pf, [&](int _i, int _j) { return _i == lx - 1 && _j < 0.33*ly; });
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double dfdsmax = 0.0, dfdsmaxall;
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        dfds[idx] = 3.0*(imx[idx]*ux[idx] + imy[idx]*uy[idx])*(-amax/(double)(pf.lx - 1)*q*(q + 1.0)/pow(q + s[idx], 2.0));
        if (dfdsmax < fabs(dfds[idx])) {
            dfdsmax = fabs(dfds[idx]);
        }
    }
#ifdef _USE_MPI_DEFINES
    MPI_Allreduce(&dfdsmax, &dfdsmaxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
    dfdsmaxall = dfdsmax;
#endif
    for (int idx = 0; idx < pf.nxyz; ++idx) {  
        dfds[idx] /= dfdsmaxall;
    }

    //--------------------Export result--------------------
    VTKXMLExport file(pf, "result/adjoint");
    file.AddPointData(pf, "p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j)];  });
    file.AddPointData(pf, "ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });
    
    delete[] rho, ux, uy, irho, iux, iuy, imx, imy, s, alpha, dfds;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}