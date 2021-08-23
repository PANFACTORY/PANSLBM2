#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main() {
    VTKXMLExport model1("result/testxml", 0, 5, 5, 1, 2, 1, 1);
    VTKXMLExport model2("result/testxml", 1, 5, 5, 1, 2, 1, 1);
    
    model1.AddPointScaler("u", [](int _i, int _j, int _k) { return _i;   });
    model2.AddPointScaler("u", [](int _i, int _j, int _k) { return _i;   });

    model1.AddPointScaler("v", [](int _i, int _j, int _k) { return _j;   });
    model2.AddPointScaler("v", [](int _i, int _j, int _k) { return _j;   });
}