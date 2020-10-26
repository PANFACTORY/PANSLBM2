#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 100;
    double nu = 0.1, u0 = 0.1;
    
    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, BARRIER);
        particle.SetBoundary(nx - 1, j, BARRIER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, OTHER);
    }
    
    NS<double, D2Q9> dsolver(&particle, nu);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        dsolver.UpdateMacro();      //  Update macroscopic values
        dsolver.Collision();        //  Collision
        particle.Stream();          //  Stream
        for (int i = 0; i < nx; i++) {
            particle.SetU(i, ny - 1, u0, 0.0);
        }                           //  Boundary condition (inlet)
        dsolver.ExternalForce();    //  External force by thermal
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/ns" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return dsolver.GetRho(particle.GetIndex(_i, _j)); });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return dsolver.GetU(0, particle.GetIndex(_i, _j)); },
                [&](int _i, int _j, int _k) { return dsolver.GetU(1, particle.GetIndex(_i, _j)); },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j); });
        } 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}