#include <cmath>
#include <cassert>
#include "../src/particle/d2q9.h"
#include "../src/utility/vtkimport.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    assert(argc == 2);

    VTKImport model(argv[1]);
    int nx = model.GetNx(), ny = model.GetNy(), nxyz = nx*ny;
    D2Q9<double> pf(nx, ny);

    double *ss = new double[nxyz], *tem = new double[nxyz], *qx = new double[nxyz], *qy = new double[nxyz];
    model.GetPointScalar("ss", ss);
    model.GetPointScalar("T", tem);
    model.GetPointVector("q", qx, qy);
    
    double S = 0.0, V = 0.0;
    double *dssdx = new double[nxyz], *dssdy = new double[nxyz], *htc = new double[nxyz];
    for (int idx = 0; idx < nxyz; ++idx) {
        dssdx[idx] = 0.0;
        dssdy[idx] = 0.0;
        htc[idx] = 0.0;

        if (ss[idx] < 0.5) {
            V += 1.0;
        }
    } 
    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            int idx = pf.Index(i, j);

            for (int c = 0; c < pf.nc; ++c) {
                dssdx[idx] += pf.cx[c]*ss[pf.Index(i + pf.cx[c], j + pf.cy[c])];
                dssdy[idx] += pf.cy[c]*ss[pf.Index(i + pf.cx[c], j + pf.cy[c])];
            }
            dssdx[idx] /= 6.0;
            dssdy[idx] /= 6.0;
            double dssnorm = sqrt(pow(dssdx[idx], 2.0) + pow(dssdy[idx], 2.0));
            if (ss[idx] > 0.5 || dssnorm < 0.1) {
                dssdx[idx] = 0.0;
                dssdy[idx] = 0.0;
            } else {
                dssdx[idx] /= dssnorm;
                dssdy[idx] /= dssnorm;
                htc[idx] = (dssdx[idx]*qx[idx] + dssdy[idx]*qy[idx])/(tem[idx] - 0.0);
                S += 1.0; 
            }
        }
    }

    VTKExport file("result/heatsink_biot.vtk", nx, ny);
    file.AddPointScaler("htc", [&](int _i, int _j, int _k) { return htc[pf.Index(_i, _j)]; });
    file.AddPointVector("dss", 
        [&](int _i, int _j, int _k) { return dssdx[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return dssdy[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return 0.0; }
    );

    std::cout << V/S << std::endl;

    delete[] ss; delete[] tem; delete[] qx; delete[] qy; delete[] dssdx; delete[] dssdy; delete[] htc;
}