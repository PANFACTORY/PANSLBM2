#include <iostream>
#include <cmath>
#include "mpi.h"

#define _USE_MPI_DEFINES
#include "../src/particle/d3q15.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int mx = atoi(argv[1]), my = atoi(argv[2]), mz = atoi(argv[3]);
    assert(mx*my*mz == PeTot);

    //--------------------Set parameters--------------------
    int lx = 101, ly = 51, lz = 51, nt = 10000, dt = 100;
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

    pf.SetBoundary([&](int _i, int _j, int _k) { return (_j == 0 || _k == 0) ? 2 : (pow(_j, 2.0) + pow(_k, 2.0) >= pow(0.15*pf.lx, 2.0) ? 1 : 0); });
    int *boundaryup = new int[pf.nbc];
    pf.SetBoundary(boundaryup, [&](int _i, int _j, int _k) { 
        if (_i == 0 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*pf.lx, 2.0)) {
            return 1;
        } else if (_i == pf.lx - 1 && pow(_j, 2.0) + pow(_k, 2.0) < pow(0.15*pf.lx, 2.0)) {
            return 2;
        } else {
            return 0;
        }
    });
    double *uxbc = new double[pf.nbc], *uybc = new double[pf.nbc], *uzbc = new double[pf.nbc], *rhobc = new double[pf.nbc], *usbc = new double[pf.nbc], *utbc = new double[pf.nbc];
    for (int idxbc = 0; idxbc < pf.nbc; ++idxbc) {
        uxbc[idxbc] = 0.0;  uybc[idxbc] = 0.0;  uzbc[idxbc] = 0.0;  rhobc[idxbc] = 1.0; usbc[idxbc] = 0.0;  utbc[idxbc] = 0.0;
    }
    for (int j = 0; j < pf.ny; ++j) {
        for (int k = 0; k < pf.nz; ++k) {
            double r = sqrt(pow(j + pf.offsety, 2.0) + pow(k + pf.offsetz, 2.0));
            if (r < 0.15*pf.lx) {
                uxbc[pf.IndexBCx(j, k) + pf.offsetxmin] = -u0*(r - 0.33*pf.ly)*(r + 0.33*pf.ly)/(0.15*pf.lx*0.15*pf.lx);
            }
        }
    }

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy, uz);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, uz, nu, alpha, true);
        pf.Swap();
        pf.Synchronize();
        pf.BoundaryCondition();
        NS::BoundaryConditionSetU(pf, uxbc, uybc, uzbc, boundaryup);
        NS::BoundaryConditionSetRho(pf, rhobc, usbc, utbc, boundaryup);
        pf.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    ANS::InitialCondition(pf, ux, uy, uz, irho, iux, iuy, iuz);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        ANS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, uz, irho, iux, iuy, iuz, imx, imy, imz, nu, alpha, true);
        pf.Swap();
        pf.iSynchronize();
        pf.iBoundaryCondition();
        ANS::BoundaryConditionSetiU(pf, uxbc, uybc, uzbc, boundaryup);
        ANS::BoundaryConditionSetiRho3D<double>(pf, boundaryup);
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
    MPI_Allreduce(&dfdsmax, &dfdsmaxall, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    for (int idx = 0; idx < pf.nxyz; ++idx) {  
        dfds[idx] /= dfdsmaxall;
    }

    //--------------------Export result--------------------
    VTKXMLExport file("result/adjoint3D", MyRank, lx, ly, lz, mx, my, mz);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j, _k)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uz[pf.Index(_i, _j, _k)]; }
    );
    //file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j, _k)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j, _k)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return iuz[pf.Index(_i, _j, _k)]; }
    );
    //file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j, _k)]; });
    //file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j, _k)]; });
    
    delete[] rho, ux, uy, uz, irho, iux, iuy, iuz, imx, imy, imz, s, alpha, dfds;
    delete[] boundaryup, uxbc, uybc, uzbc, rhobc, usbc, utbc;

    MPI_Finalize();

    return 0;
}