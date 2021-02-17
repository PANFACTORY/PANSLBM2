#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    int tmax = 100000, nx = 100, ny = 50, nb = 64;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0, q0, Td = 1.7, v0 = 0.0, dv = 62.8/(double)nb;
    double bpx[nb], bpy[nb], bux[nb], buy[nb], bgx[nb], bgy[nb], bvx[nb], bvy[nb], btem[nb], bgtem[nb], btemd[nb];

    for (int k = 0; k < nb; k++) {
        bpx[k] = 10.0*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0;
        bpy[k] = 10.0*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0;
        bvx[k] = -v0*sin(2.0*M_PI*k/(double)nb);
        bvy[k] = v0*cos(2.0*M_PI*k/(double)nb);
        btemd[k] = Td;
    }

    VTKExportIB fout("result/model.vtk");
    fout.AddPoint(nb, bpx, bpy);

    return 0;
}