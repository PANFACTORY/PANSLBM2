#define _USE_MATH_DEFINES
#define _USE_MPI_DEFINES
#include <cmath>
#include <iostream>
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
    int lx = 30, ly = 30, lz = 30, nt = 10000, dt = 100;
    double nu = 0.1, u0 = 0.1, theta = 90.0;
    D3Q15<double> pf(lx, ly, lz, MyRank, mx, my, mz);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz];
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;  uy[idx] = 0.0;  uz[idx] = 0.0;
    }

    pf.SetBoundary([&](int _i, int _j, int _k) { return (_i == 0 || _i == pf.lx - 1 || _j == 0 || _j == pf.ly - 1 || _k == 0) ? 1 : 0; });
    int *boundaryu = new int[pf.nbc];
    pf.SetBoundary(boundaryu, [&](int _i, int _j, int _k) { return _k == pf.lz - 1 ? 1 : 0; });
    double *uxbc = new double[pf.nbc], *uybc = new double[pf.nbc], *uzbc = new double[pf.nbc];
    for (int idxbc = 0; idxbc < pf.nbc; ++idxbc) {
        uxbc[idxbc] = u0*cos(theta*M_PI/180.0);
        uybc[idxbc] = u0*sin(theta*M_PI/180.0);
        uzbc[idxbc] = 0.0;
    }

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy, uz);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt == 0 && MyRank == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        NS::Macro_Collide_Stream(pf, rho, ux, uy, uz, nu, true);
        pf.Swap();
        pf.Synchronize();
        pf.BoundaryCondition();
        NS::BoundaryConditionSetU(pf, uxbc, uybc, uzbc, boundaryu);
        pf.SmoothCorner();
    }

    //--------------------Export result--------------------
    VTKXMLExport file("result/cavity3D", MyRank, lx, ly, lz, mx, my, mz);
    file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j, _k)]; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j, _k)]; },
        [&](int _i, int _j, int _k) { return uz[pf.Index(_i, _j, _k)]; }
    );
    file.AddPointScaler("bctype", [&](int _i, int _j, int _k) { return pf.GetBoundary(_i, _j, _k); });

    delete[] rho, ux, uy, uz, uxbc, uybc, uzbc;
    delete[] boundaryu;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}