#define _USE_MPI_DEFINES
#include <iostream>
#include <chrono>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
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
    int lx = 101, ly = 51, nt = 100000, dt = 1000;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    D2Q9<double> pg(lx, ly, MyRank, mx, my);
    double rho[pf.nxy], ux[pf.nxy], uy[pf.nxy], tem[pg.nxy], qx[pg.nxy], qy[pg.nxy];
    for (int idx = 0; idx < pf.nxy; ++ idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;
        tem[idx] = 0.5*(Th + Tl);   qx[idx] = 0.0;  qy[idx] = 0.0;
    }

    pf.SetBoundary([&](int _i, int _j) {    return (_j == 0 || _j == pf.ly - 1) ? 2 : 0;    });
    pg.SetBoundary([&](int _i, int _j) {    return 0;    });
    int boundaryt[pg.nbc];
    pg.SetBoundary(boundaryt, [&](int _i, int _j) {    return (_j == 0 || _j == pg.ly - 1) ? 1 : 0;    });
    double tembc[pg.nbc] = { 0 };
    for (int i = 0; i < pg.nx; ++i) {
        tembc[i + pg.offsetymin] = Th;
        tembc[i + pg.offsetymax] = Tl;
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::InitialCondition(pf, rho, ux, uy);
    AD::InitialCondition(pg, tem, ux, uy);
    for (int t = 1; t <= nt; ++t) {
        if (t%dt == 0) {
            std::cout << "t = " << t/dt << std::endl;
        }
        AD::Macro_Collide_Stream_NaturalConvection(
            pf, rho, ux, uy, nu,
            pg, tem, qx, qy, alpha, 
            0.0, 1.6e-5, 0.5*(Th + Tl), true
        );
        pf.Swap();
        pg.Swap();
        pf.Synchronize();
        pg.Synchronize();
        pf.BoundaryCondition();
        pg.BoundaryCondition();
        AD::BoundaryConditionSetT(pg, tembc, ux, uy, boundaryt);
        pf.SmoothCorner();
        pg.SmoothCorner();
    }

    //--------------------Export result--------------------
    VTKXMLExport file("result/naturalconvection", MyRank, lx, ly, lz, mx, my, mz);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    //file.AddPointScaler("bctype", [&](int _i, int _j, int _k) { return pf.GetBoundary(_i, _j, _k); });
    file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
    file.AddPointVector("q", 
        [&](int _i, int _j, int _k) { return qx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}