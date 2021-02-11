#include "vtkimport.h"

using namespace PANSLBM2;

int main() {
    VTKImport model("result/ns/ns99.vtk");
    
    double ux[10000], uy[10000];
    model.GetPointVector("u", ux, uy);

    for (int i = 0; i < 10000; i++) {
        std::cout << ux[i] << "\t" << uy[i] << std::endl;
    }
}