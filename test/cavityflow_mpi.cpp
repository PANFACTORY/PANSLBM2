#include <iostream>
#include <chrono>
#include <cassert>
#include "mpi.h"

#include "../src/particle/d2q9mpi.h"
#include "../src/equation/navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
    assert(PeTot == 2);

    //--------------------Set parameters--------------------
    int nx = 51, ny = 101, nt = 100000, dt = 1000;
    double nu = 0.1, u0 = 0.1, Re = u0*(nx - 1)/nu;
    D2Q9<double> pf(nx, ny);
    double rho[pf.nxy], ux[pf.nxy], uy[pf.nxy];
    for (int idx = 0; idx < pf.nxy; ++idx) {
        rho[idx] = 1.0;
        ux[idx] = 0.0;
        uy[idx] = 0.0;
    }

    int boundaryu[pf.nxy] = { 0 };
    double uxbc[pf.nbc] = { 0 }, uybc[pf.nbc] = { 0 };
    for (int idxbc = 0; idxbc < pf.nbc; ++idxbc) {
        uxbc[idxbc] = u0;
        uybc[idxbc] = 0.0;
    }  
    if (MyRank == 0) {
        pf.SetNeighborId(-1, 1, -1, -1, -1, -1, -1, -1);
        pf.SetBoundaryAlongXmin([&](int _j) {   return 1;   });
        pf.SetBoundaryAlongYmin([&](int _i) {   return 1;   });
        pf.SetBoundaryAlongYmax(boundaryu, [&](int _i) {    return 1;   });        
    } else if (MyRank == 1) {
        pf.SetNeighborId(0, -1, -1, -1, -1, -1, -1, -1);
        pf.SetBoundaryAlongXmax([&](int _j) {   return 1;   });
        pf.SetBoundaryAlongYmin([&](int _i) {   return 1;   });
        pf.SetBoundaryAlongYmax(boundaryu, [&](int _i) {    return 1;   });        
    }
     
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt != 0) {
            NS::Macro_Collide_Stream(pf, rho, ux, uy, nu);
        } else {
            NS::Macro_Collide_Stream(pf, rho, ux, uy, nu, true);

            if (MyRank == 0) {
                std::cout << "t = " << t/dt << std::endl;
            }
            VTKExport file("result/ns_at" + std::to_string(MyRank) + "_" + std::to_string(t/dt) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
        }

        pf.Swap();
        pf.Synchronize();
        pf.BoundaryCondition();
        NS::BoundaryConditionSetU(pf, uxbc, uybc, boundaryu);
        pf.SmoothCorner();
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }
    
    MPI_Finalize();
}