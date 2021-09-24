#define _USE_MPI_DEFINES
#include <iostream>
#include <cmath>
#include <chrono>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/equation/adjointelastic.h"
#include "../src/utility/vtkexport.h"

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

    //--------------------Setting parameters--------------------
    int lx = 81, ly = 61, nt = 10000, dt = 1000;
    double rho0 = 1e6, stress0 = -1.0, p = 0.25, smin = 0.1, smax = 0.9;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    double rho[pf.nxyz], ux[pf.nxyz], uy[pf.nxyz], rx[pf.nxyz], ry[pf.nxyz], sxx[pf.nxyz], sxy[pf.nxyz], syx[pf.nxyz], syy[pf.nxyz];
    double irho[pf.nxyz], imx[pf.nxyz], imy[pf.nxyz], isxx[pf.nxyz], isxy[pf.nxyz], isyx[pf.nxyz], isyy[pf.nxyz];
    double s[pf.nxyz], gamma[pf.nxyz], sensitivity[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            s[idx] = pow((i + _p.offsetx) - 0.5*pf.lx, 2.0) + pow((j + _p.offsety), 2.0) < pow(0.15*pf.lx, 2.0) ? smin : smax;
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = rho0;    
        ux[idx] = 0.0;  uy[idx] = 0.0;
        rx[idx] = 0.0;  ry[idx] = 0.0;
        sxx[idx] = 0.0; sxy[idx] = 0.0; syx[idx] = 0.0; syy[idx] = 0.0;
        irho[idx] = 0.0;
        imx[idx] = 0.0; imy[idx] = 0.0;
        isxx[idx] = 0.0;    isxy[idx] = 0.0;    isyx[idx] = 0.0;    isyy[idx] = 0.0;
        gamma[idx] = pow(s[idx], p);
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    
    //--------------------Direct analyze--------------------
    EL::InitialCondition(pf, rho, ux, uy, sxx, sxy, syx, syy);
    for (int t = 1; t <= nt; ++t) {
        EL::MacroExtendedCollideStream(pf, rho, ux, uy, sxx, sxy, syx, syy, 0.8, gamma, true);
        pf.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return _i == 0 ? 1 : 0; });
        EL::BoundaryConditionSetStress(pf, 
            [=](int _i, int _j) { return (_i == lx - 1 && fabs(_j - 0.5*ly) < 0.08*ly) ? stress0 : 0.0; }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _i != 0; }
        );
        pf.SmoothCorner();

        for (int idx = 0; idx < pf.nxyz; ++idx) {
            rx[idx] += ux[idx]; ry[idx] += uy[idx];
        }
    }

    //--------------------Invert analyze--------------------
    AEL::InitialCondition(pf, irho, imx, imy, isxx, isxy, isyx, isyy, gamma);
    for (int t = 1; t <= nt; ++t) {
        AEL::MacroCollideStream(pf, irho, imx, imy, isxx, isxy, isyx, isyy, 0.8, gamma, true);
        pf.iStream();
        pf.iBoundaryCondition([=](int _i, int _j) { return _i == 0 ? 1 : 0; });
        AEL::iBoundaryConditionSetStress(pf, 
            [=](int _i, int _j) { return (_i == lx - 1 && fabs(_j - 0.5*ly) < 0.08*ly) ? stress0 : 0.0; }, 
            [=](int _i, int _j) { return 0.0; }, 
            rho, 
            [=](int _i, int _j) { return _i != 0; }
        );
        pf.SmoothCorner();
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        sensitivity[idx] = (4.5*(sxx[idx]*isxx[idx] + sxy[idx]*isxy[idx] + syx[idx]*isyx[idx] + syy[idx]*isyy[idx]) - 1.5*irho[idx]*(sxx[idx] + syy[idx]));
        if (sensitivitymax < fabs(sensitivity[idx])) {
            sensitivitymax = fabs(sensitivity[idx]);
        }
    }
    for (int idx = 0; idx < pf.nxyz; ++idx) {  
        sensitivity[idx] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKXMLExport file(pf, "result/adjointelastic");
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "r", 
        [&](int _i, int _j, int _k) { return rx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return ry[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "s", 
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
    file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "is", 
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
    file.AddPointData(pf, "sensitivity", [&](int _i, int _j, int _k) { return sensitivity[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "gamma", [&](int _i, int _j, int _k) { return gamma[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "ss", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}