#include <cassert>

#include "../src/particle/d3q15.h"
#include "../src/utility/vtkxmlimport.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    assert(argc == 4);
    VTKXMLImport modelin(argv[1]);
    D3Q15<double> p(modelin.GetNx(), modelin.GetNy(), modelin.GetNz());
    std::vector<double> v(p.nxyz);
    modelin.GetPointScalar(p, argv[2], v.data());
    VTKExport modelout(argv[3], p.nx, p.ny, p.nz);
    modelout.AddPointScaler(argv[2], [&](int _i, int _j, int _k){ return v[p.Index(_i, _j, _k)]; });
}