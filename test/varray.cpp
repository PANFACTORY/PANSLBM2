#include <iostream>
#include "../src/utility/varray.h"

using namespace PANSLBM2;

int main() {
    const int n = 10;
    VArray<double> a(n);
    std::cout << "array generated" << std::endl;
    
    //  Write values
    for (int i = 0; i < n; ++i) {
        a[i] = i + 1;
    }

    //  Read values
    for (int i = 0; i < n; ++i) {
        std::cout << a[i] << std::endl;
    }
}