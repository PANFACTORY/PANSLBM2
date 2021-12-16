#include "../src/particle/d3q15.h"
#include "../src/utility/vtkxmlimport.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    VTKXMLImport modelin("result/heavisidefilter.pvts");
    D3Q15<double> p(modelin.GetNx(), modelin.GetNy(), modelin.GetNz());
    std::vector<double> v(p.nxyz);
    modelin.GetPointScalar(p, "v", v.data());
    VTKExport modelout("result/heavisidefilter.vtk", p.nx, p.ny, p.nz);
    modelout.AddPointScaler("v", [&](int _i, int _j, int _k){ return v[p.Index(_i, _j, _k)]; });
}