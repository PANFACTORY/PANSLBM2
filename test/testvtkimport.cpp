#include "vtkimport.h"

using namespace PANSLBM2;

int main() {
    VTKImport model("result/elastic/elastic99.vtk");
    
    double sxx[10000], sxy[10000], syy[10000];
    model.GetPointTensor("s", sxx, sxy, sxy, syy);

    for (int i = 0; i < 10000; i++) {
        std::cout << sxx[i] << "\t" << sxy[i] << "\t" << syy[i] << std::endl;
    }
}