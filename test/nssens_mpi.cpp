#include <iostream>
#include <cmath>
#include "mpi.h"

#define _USE_MPI_DEFINES
#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int mx = atoi(argv[1]), my = atoi(argv[2]);
    assert(mx*my == PeTot);

    //--------------------Set parameters--------------------
    int lx = 101, ly = 51, nt = 10000, dt = 100;
    double nu = 0.1, u0 = 0.0109, rho0 = 1.0, epsdu = 1.0e-4, q = 1e-2, amax = 50.0, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(lx, ly);
    double *rho = new double[pf.nxy], *ux = new double[pf.nxy], *uy = new double[pf.nxy];
    double *irho = new double[pf.nxy], *iux = new double[pf.nxy], *iuy = new double[pf.nxy], *imx = new double[pf.nxy], *imy = new double[pf.nxy];
    double *s = new double[pf.nxy], *alpha = new double[pf.nxy], *sensitivity = new double[pf.nxy];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow(i - 0.5*pf.lx, 2.0) + pow(j, 2.0) < pow(0.15*pf.lx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxy; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0;
        alpha[idx] = amax/(double)(pf.nx - 1)*q*(1.0 - s[idx])/(s[idx] + q);
    }

    pf.SetBoundary([&](int _i, int _j) {    return _j == 0 ? 2 : (_j >= 0.33*pf.ly ? 1 : 0); });
    int *boundaryup = new int[pf.nxy];
    pf.SetBoundary(boundaryup, [&](int _i, int _j) {    return (_i == 0 && _j < 0.33*pf.ly) ? 1 : ((_i == pf.lx - 1 && _j < 0.33*pf.ly) ? 2 : 0); });
    double *uxbc = new double[pf.nxy], *uybc = new double[pf.nxy], *rhobc = new double[pf.nxy], *usbc = new double[pf.nxy];
    for (int j = 0; j < pf.ny; ++j) {
        uxbc[j + pf.offsetxmin] = (j + pf.offsetx) < 0.33*pf.ly ? -u0*((j + pf.offsetx) - 0.33*pf.ly)*((j + pf.offsetx)  + 0.33*pf.ly)/(0.33*pf.ly*0.33*pf.ly) : 0.0;
        uybc[j + pf.offsetxmin] = 0.0;  rhobc[j + pf.offsetxmin] = 1.0; usbc[j + pf.offsetxmin] = 0.0;
        uxbc[j + pf.offsetxmax] = 0.0;  uybc[j + pf.offsetxmax] = 0.0;  rhobc[j + pf.offsetxmax] = 1.0; usbc[j + pf.offsetxmax] = 0.0;
    }
    for (int i = 0; i < pf.nx; ++i) {
        uxbc[i + pf.offsetymin] = 0.0;  uybc[i + pf.offsetymin] = 0.0;  rhobc[i + pf.offsetymin] = 1.0; usbc[i + pf.offsetymin] = 0.0;
        uxbc[i + pf.offsetymax] = 0.0;  uybc[i + pf.offsetymax] = 0.0;  rhobc[i + pf.offsetymax] = 1.0; usbc[i + pf.offsetymax] = 0.0;
    }

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, nu, alpha, true);
        pf.Swap();
        pf.Synchronize();
        pf.BoundaryCondition();
        NS::BoundaryConditionSetU(pf, uxbc, uybc, boundaryup);
        NS::BoundaryConditionSetRho(pf, rhobc, usbc, boundaryup);
        pf.SmoothCorner();
    }

    //--------------------Invert analyze--------------------
    ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
    for (int t = 1; t <= nt; ++t) {
        if (MyRank == 0 && t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        ANS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, irho, iux, iuy, imx, imy, nu, alpha, true);
        pf.Swap();
        pf.iSynchronize();
        pf.iBoundaryCondition();
        ANS::BoundaryConditionSetiU(pf, uxbc, uybc, boundaryup);
        ANS::BoundaryConditionSetiRho<double>(pf, boundaryup);
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int idx = 0; idx < pf.nxy; ++idx) {
        sensitivity[idx] = 3.0*(imx[idx]*ux[idx] + imy[idx]*uy[idx])*(-amax/(double)(pf.nx - 1)*q*(q + 1.0)/pow(q + s[idx], 2.0));
        if (sensitivitymax < fabs(sensitivity[idx])) {
            sensitivitymax = fabs(sensitivity[idx]);
        }
    }
    for (int idx = 0; idx < pf.nxy; ++idx) {  
        sensitivity[idx] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointns_at" + std::to_string(MyRank) + ".vtk", pf.nx, pf.ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[pf.Index(_i, _j)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });
    
    delete[] rho, ux, uy, irho, iux, iuy, imx, imy, s, alpha, sensitivity;
    delete[] boundaryup, uxbc, uybc, rhobc, usbc;

    MPI_Finalize();

    return 0;
}