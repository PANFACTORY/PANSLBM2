#include "../src/utility/vtkexport.h"
#include <iostream>

using namespace PANSLBM2;

int main() {
    const int lx = 141, ly = 161, h = 101, w = 81, a = 21, n = 5;
    double weightlimit = 0.5;
    int b = (weightlimit*(2*(w - 1) + 1)*h - a*(2*(w - 1) + 1))/((h - a)*n) + 1;
    VTKExport file("result/heatsink_pinfin.vtk", lx, ly);
    file.AddPointScaler("ss", [&](int _i, int _j, int _k) { 
        return (_i < w && _j < h && (_j < a || (w - _i - 1)%(((2*(w - 1) + 1) - b)/(n - 1)) < b)) ? 0.0 : 1.0; 
    });
    std::cout << (a*(2*(w - 1) + 1) + (h - a)*n*b)/(double)(h*(2*(w - 1) + 1)) << std::endl;
    return 0;
}