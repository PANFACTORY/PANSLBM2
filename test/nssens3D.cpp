#define _USE_MPI_DEFINES
#include <iostream>
#include <cmath>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d3q15.h"
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

    assert(argc == 4);
    int mx = atoi(argv[1]), my = atoi(argv[2]), mz = atoi(argv[3]);
    assert(mx*my*mz == PeTot);
#else
    int MyRank = 0, mx = 1, my = 1, mz = 1.0;
#endif

    //--------------------Set parameters--------------------
    int lx = 61, ly = 31, lz = 31, nt = 10000, dt = 100;
    double nu = 0.1, u0 = 0.0109, rho0 = 1.0, epsdu = 1.0e-4, q = 1e-2, amax = 50.0, smin = 0.1, smax = 0.9;
    D3Q15<double> pf(lx, ly, lz, MyRank, mx, my, mz);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *iuz = new double[pf.nxyz];
    double *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *imz = new double[pf.nxyz];
    double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *dfds = new double[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            for (int k = 0; k < pf.nz; ++k) {
                int idx = pf.Index(i, j, k);
                s[idx] = pow((i + pf.offsetx) - 0.5*pf.lx, 2.0) + pow((j + pf.offsety), 2.0) + pow((k + pf.offsetz), 2.0) < pow(0.15*pf.lx, 2.0) ? smin : smax;
            }
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;  uz[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; iuz[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; imz[idx] = 0.0;
        alpha[idx] = amax/(double)(pf.lx - 1)*q*(1.0 - s[idx])/(s[idx] + q);
    }

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy, uz);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::MacroBrinkmanCollideStream(pf, rho, ux, uy, uz, nu, alpha, true);
        pf.Swap();
        pf.Synchronize();
        pf.BoundaryCondition([=](int _i, int _j, int _k) { return (_j == 0 || _k == 0) ? 2 : (pow(_j, 2.0) + pow(_k, 2.0) >= pow(0.15*lx, 2.0) ? 1 : 0); });
        NS::BoundaryConditionSetU(pf, 
            [=](int _i, int _j, int _k) { 
                double r = sqrt(pow(_j, 2.0) + pow(_k, 2.0));
                return -u0*(r - 0.15*lx)*(r + 0.15*lx)/(0.15*lx*0.15*lx); 
            }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return _i == 0 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*lx, 2.0); }
        );
        NS::BoundaryConditionSetRho(pf, 
            [=](int _i, int _j, int _k) { return rho0; },  
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return _i == lx - 1 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*lx, 2.0); }
        );
        pf.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    ANS::InitialCondition(pf, ux, uy, uz, irho, iux, iuy, iuz);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        ANS::MacroBrinkmanCollideStream(pf, rho, ux, uy, uz, irho, iux, iuy, iuz, imx, imy, imz, nu, alpha, true);
        pf.Swap();
        pf.iSynchronize();
        pf.iBoundaryCondition([=](int _i, int _j, int _k) { return (_j == 0 || _k == 0) ? 2 : (pow(_j, 2.0) + pow(_k, 2.0) >= pow(0.15*lx, 2.0) ? 1 : 0); });
        ANS::iBoundaryConditionSetU(pf, 
            [=](int _i, int _j, int _k) { 
                double r = sqrt(pow(_j, 2.0) + pow(_k, 2.0));
                return -u0*(r - 0.15*lx)*(r + 0.15*lx)/(0.15*lx*0.15*lx); 
            }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return _i == 0 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*lx, 2.0); },
            1.0
        );
        ANS::iBoundaryConditionSetRho3D(pf, [=](int _i, int _j, int _k) { return _i == lx - 1 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*lx, 2.0); });
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double dfdsmax = 0.0, dfdsmaxall;
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        dfds[idx] = 3.0*(imx[idx]*ux[idx] + imy[idx]*uy[idx] + imz[idx]*uz[idx])*(-amax/(double)(pf.lx - 1)*q*(q + 1.0)/pow(q + s[idx], 2.0));
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
    VTKXMLExport file(pf, "result/adjoint3D");
    file.AddPointData(pf, "p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j, _k)]/3.0; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uz[pf.Index(_i, _j, _k)]; }
    );
    file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j, _k)];  });
    file.AddPointData(pf, "ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j, _k)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return iuz[pf.Index(_i, _j, _k)]; }
    );
    file.AddPointData(pf, "alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j, _k)]; });
    file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j, _k)]; });
    
    delete[] rho, ux, uy, uz, irho, iux, iuy, iuz, imx, imy, imz, s, alpha, dfds;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}