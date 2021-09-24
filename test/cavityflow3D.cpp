#define _USE_MATH_DEFINES
//#define _USE_MPI_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d3q15.h"
#include "../src/equation/navierstokes.h"
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
    int MyRank = 0, mx = 1, my = 1, mz = 1;
#endif

    //--------------------Set parameters--------------------
    int lx = 31, ly = 31, lz = 31, nt = 10000, dt = 100;
    double nu = 0.1, u0 = 0.1, theta = 90.0;
    D3Q15<double> pf(lx, ly, lz, MyRank, mx, my, mz);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz];
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;  uy[idx] = 0.0;  uz[idx] = 0.0;
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy, uz);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::MacroCollideStream(pf, rho, ux, uy, uz, nu, true);
        pf.Stream();
        pf.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _i == lx - 1 || _j == 0 || _j == ly - 1 || _k == 0) ? 1 : 0; });
        NS::BoundaryConditionSetU(pf, 
            [=](int _i, int _j, int _k) { return u0*cos(theta*M_PI/180.0); }, 
            [=](int _i, int _j, int _k) { return u0*sin(theta*M_PI/180.0); }, 
            [=](int _i, int _j, int _k) { return 0.0; }, 
            [=](int _i, int _j, int _k) { return _k == lz - 1; }
        );
        pf.SmoothCorner();
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }

    //--------------------Export result--------------------
    VTKXMLExport file(pf, "result/cavity3D");
    file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j, _k)]; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uz[pf.Index(_i, _j, _k)]; }
    );
    
    delete[] rho, ux, uy, uz;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}