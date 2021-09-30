#include <iostream>
#include "../src/utility/varray.h"

using namespace PANSLBM2;

int main() {
    const long n = 100000UL*141*161*8;
    VArray<double> a(n);
    
    //  Write values
    for (long i = 0; i < n; ++i) {
        a[i] = i + 1;
    }

    //  Read values
    for (long i = 0; i < n; ++i) {
        std::cout << a[i] << std::endl;
    }
}