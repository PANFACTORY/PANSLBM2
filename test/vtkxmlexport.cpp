#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main() {
    VTKXMLExport model1("result/testxml", 0, 5, 5, 1, 2, 1, 1);
    VTKXMLExport model2("result/testxml", 1, 5, 5, 1, 2, 1, 1);
    
    model1.AddPointData("u", 
        [](int _i, int _j, int _k) { return _i;   }, 
        [](int _i, int _j, int _k) { return _j;   }, 
        [](int _i, int _j, int _k) { return _k;   }
    );
    model2.AddPointData("u", 
        [](int _i, int _j, int _k) { return _i;   }, 
        [](int _i, int _j, int _k) { return _j;   }, 
        [](int _i, int _j, int _k) { return _k;   }
    );

    model1.AddPointData("v", [](int _i, int _j, int _k) { return _j;   });
    model2.AddPointData("v", [](int _i, int _j, int _k) { return _j;   });
}