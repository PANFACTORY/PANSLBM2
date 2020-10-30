#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/advection.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    double rho[nx*ny], u[nx*ny], v[nx*ny], T[nx*ny];
    
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
    }                                               //  Set boundary condition
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 1.5, 0.0, 0.0);
    }                                               //  Set initial condition
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        NS::UpdateMacro(particlef, rho, u, v);
        AD::UpdateMacro(particleg, T);              //  Update macroscopic values
        NS::Collision(nu, particlef, rho, u, v);
        AD::Collision(alpha, particleg, T, u, v);   //  Collision
        particlef.Stream();
        particleg.Stream();                         //  Stream
        for (int i = 0; i < nx; i++) {
            particleg.SetTemperature(i, 0, Th);
            particleg.SetTemperature(i, ny - 1, Tl);
        }                                           //  Boundary condition (Fix temperature)
        NS::ExternalForceNaturalConvection(0.0, 1.6e-5, 1.5, particlef, particleg);    //  External force with natural convection
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/thermal" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particlef.GetIndex(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return u[particlef.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return v[particlef.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return T[particleg.GetIndex(_i, _j)]; });
            file.AddPointScaler("boundaryf", [&](int _i, int _j, int _k) { return particlef.GetBoundary(_i, _j); });
            file.AddPointScaler("boundaryg", [&](int _i, int _j, int _k) { return particleg.GetBoundary(_i, _j); });
        } 
    }                                               //  Export result per 1000 steps 

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}