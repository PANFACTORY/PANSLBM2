#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 1000, ny = 1000;
    double nu = 0.005, u0 = 0.05;
    double *rho = new double[nx*ny], *u = new double[nx*ny], *v = new double[nx*ny];
    
    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, BARRIER);
        particle.SetBoundary(nx - 1, j, BARRIER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, OTHER);
    }                                               //  Set boundary condition
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particle, 1.0, 0.0, 0.0);
    }                                               //  Set initial condition
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        NS::UpdateMacro(particle, rho, u, v);       //  Update macroscopic values
        //NS::CollisionMRT(1.19, 1.4, 1.2, 1.2, 1.0/(3.0*nu + 0.5), 1.0/(3.0*nu + 0.5), particle, rho, u, v);     //  Collision
        //NS::CollisionMRT(1.0, 1.0, 1.0, 1.0, 1.0/(3.0*nu + 0.5), 1.0/(3.0*nu + 0.5), particle, rho, u, v);     //  Collision
        //NS::Collision(nu, particle, rho, u, v);     //  Collision
        NS::CollisionTRT(nu, 0.25, particle, rho, u, v);     //  Collision
        particle.Stream();                          //  Stream
        for (int i = 0; i < nx; i++) {
            particle.SetU(i, ny - 1, u0, 0.0);
        }                                           //  Boundary condition (inlet)
        particle.SmoothCorner();
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/TRT_" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particle.GetIndex(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return u[particle.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return v[particle.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j); });
        }                                           //  Export result per 1000 steps 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    delete[] rho, u, v;
}