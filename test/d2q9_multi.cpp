#include "../src/particle/d2q9_multi.h"

using namespace PANSLBM2;

int main() {
    D2Q9_multi<double> p(2, 3);
    p.GenerateChild(1, 2);
}