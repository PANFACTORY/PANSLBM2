#define _USE_MPI_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d3q15.h"
#include "../src/utility/densityfilter.h"
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

    int lx = 51, ly = 51, lz = 51;
    double R = 1.8;
    D3Q15<double> pf(lx, ly, lz, MyRank, mx, my, mz);
    std::vector<double> v(pf.nxyz);
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            for (int k = 0; k < pf.nz; ++k) {
                int idx = pf.Index(i, j, k);
                v[idx] = ((i + pf.offsetx) + (j + pf.offsety) + (k + pf.offsetz))%2 ? 1.0 : 0.0;
            }
            
        }
    }
    std::vector<double> fv = DensityFilter::GetFilteredValue(pf, 1.5, v);
    VTKXMLExport file("result/densityfilter", MyRank, lx, ly, lz, mx, my, mz);
    file.AddPointScaler("v", [&](int _i, int _j, int _k) { return v[pf.Index(_i, _j, _k)]; });
    file.AddPointScaler("fv", [&](int _i, int _j, int _k) { return fv[pf.Index(_i, _j, _k)]; });
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}