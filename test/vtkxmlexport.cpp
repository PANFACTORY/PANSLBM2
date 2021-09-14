#include "../src/particle/d2q9.h"
#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main() {
    D2Q9<double> p(5, 5);
    VTKXMLExport model1(p, "result/testxml");
    
    model1.AddPointData(p, "u", 
        [](int _i, int _j, int _k) { return _i;   }, 
        [](int _i, int _j, int _k) { return _j;   }, 
        [](int _i, int _j, int _k) { return _k;   }
    );

    model1.AddPointData(p, "v", [](int _i, int _j, int _k) { return _j;   });
}