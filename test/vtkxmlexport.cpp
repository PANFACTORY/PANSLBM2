#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main() {
    VTKXMLExport model1("result/testxml", 0, 10, 10, 1, 2, 1, 1);
    VTKXMLExport model2("result/testxml", 1, 10, 10, 1, 2, 1, 1);
}