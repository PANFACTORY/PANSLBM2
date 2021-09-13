#define _USE_MPI_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/utility/densityfilter.h"
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

    int lx = 101, ly = 101;
    double R = 1.5;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    std::vector<double> v(pf.nxyz);
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            v[idx] = (i + pf.offsetx)%2 == (j + pf.offsety)%2 ? 1.0 : 0.0;
        }
    }
    std::vector<double> fv = DensityFilter::GetFilteredValue(pf, 1.5, v);
    std::cout << (int)R << std::endl;
    VTKXMLExport file("result/densityfilter", MyRank, lx, ly, 1, mx, my, 1);
    file.AddPointScaler("v", [&](int _i, int _j, int _k) { return v[pf.Index(_i, _j)]; });
    file.AddPointScaler("fv", [&](int _i, int _j, int _k) { return fv[pf.Index(_i, _j)]; });
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}