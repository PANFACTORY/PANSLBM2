//  This program is unit test of particle.

#include <iostream>
#include "../src/particle.h"

using namespace PANSLBM2;

int main() {
    //  Test of Constructor()
    D2Q9<double> particle0 = D2Q9<double>(100, 100);

    //  Test of CopyConstructor()
    D2Q9<double> particle1 = particle0;

    //  Test of SetBarrier()
    particle0.SetBarrier(0, 0, true);

    //  Test of SetBoundary()
    particle0.SetBoundary(0, 0, OTHER);

    //  Test of Stream()
    particle0.Stream();

    //  Test of GetBarrier()
    
}