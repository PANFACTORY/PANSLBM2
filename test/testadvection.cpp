#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    
    D2Q9<double> particlef(nx, ny);
    D2Q9<double> particleg(nx, ny);
    for (int j = 0; j < ny; j++) {
        particlef.SetBoundary(0, j, PERIODIC);
        particlef.SetBoundary(nx - 1, j, PERIODIC);
        particleg.SetBoundary(0, j, PERIODIC);
        particleg.SetBoundary(nx - 1, j, PERIODIC);
    }
    for (int i = 0; i < nx; i++) {
        particlef.SetBoundary(i, 0, MIRROR);
        particlef.SetBoundary(i, ny - 1, MIRROR);
        particleg.SetBoundary(i, 0, OTHER);
        particleg.SetBoundary(i, ny - 1, OTHER);
    }
        
    AD<double, D2Q9, D2Q9> dsolver(&particlef, &particleg, nu, alpha);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        dsolver.UpdateMacro();      //  Update macroscopic values
        dsolver.Collision();        //  Collision
        particlef.Stream();         //  Stream
        particleg.Stream();
        for (int i = 0; i < nx; i++) {
            particleg.SetTemperature(i, 0, Th);
            particleg.SetTemperature(i, ny - 1, Tl);
        }                           //  Boundary condition (Fix temperature)
        dsolver.ExternalForce();    //  External force by thermal
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/thermal" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return dsolver.GetRho(particlef.GetIndex(_i, _j)); });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return dsolver.GetU(0, particlef.GetIndex(_i, _j)); },
                [&](int _i, int _j, int _k) { return dsolver.GetU(1, particlef.GetIndex(_i, _j)); },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return dsolver.GetTemperature(particleg.GetIndex(_i, _j)); });
            file.AddPointScaler("boundaryf", [&](int _i, int _j, int _k) { return particlef.GetBoundary(_i, _j); });
            file.AddPointScaler("boundaryg", [&](int _i, int _j, int _k) { return particleg.GetBoundary(_i, _j); });
        } 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}