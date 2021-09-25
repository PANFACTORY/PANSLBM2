//#define _USE_MPI_DEFINES
#include <iostream>
#include <chrono>
#include <cassert>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
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
    int lx = 101, ly = 101, nt = 100000, dt = 1000;
    double nu = 0.1, u0 = 0.1, Re = u0*(lx - 1)/nu;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    double rho[pf.np], ux[pf.np], uy[pf.np];
    for (int idx = 0; idx < pf.np; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;
        uy[idx] = 0.0;
    }
     
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt != 0) {
            NS_AVX::MacroCollideStream(pf, rho, ux, uy, nu);
        } else {
            NS_AVX::MacroCollideStream(pf, rho, ux, uy, nu, true);

            if (MyRank == 0) {
                std::cout << "t = " << t/dt << std::endl;
            }
            /*VTKXMLExport file(pf, "result/cavity_" + std::to_string(t/dt));
            file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointData(pf, "u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );*/
        }

        pf.Stream();
        pf.BoundaryCondition([=](int _i, int _j) { return (_i == 0 || _i == lx - 1 || _j == 0) ? 1 : 0; });
        NS::BoundaryConditionSetU(pf, 
            [=](int _i, int _j) { return u0; }, 
            [=](int _i, int _j) { return 0.0; }, 
            [=](int _i, int _j) { return _j == ly - 1; }
        );
        pf.SmoothCorner();
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}