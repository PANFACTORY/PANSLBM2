#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 300, ny = 30;
    double elasticy = 0.1, rho0 = 0.1, stress0 = -1.0;
    double rho[nx*ny], ux[nx*ny], uy[nx*ny], rx[nx*ny], ry[nx*ny], sxx[nx*ny], sxy[nx*ny], syx[nx*ny], syy[nx*ny];
    for (int i = 0; i < nx*ny; i++) {
        rho[i] = rho0;
        rx[i] = 0.0;
        ry[i] = 0.0;
    }
    
    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, OTHER);
        particle.SetBoundary(nx - 1, j, OTHER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, OTHER);
        particle.SetBoundary(i, ny - 1, OTHER);
    }                                                                       //  Set boundary condition
    for (int i = 0; i < nx*ny; i++) {
        EL::InitialCondition(i, particle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }                                                                       //  Set initial condition
    EL::UpdateMacro(particle, rho, ux, uy, sxx, sxy, syx, syy);             //  Update macroscopic values
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        EL::Collision(elasticy, particle, rho, ux, uy, sxx, sxy, syx, syy); //  Collision
        particle.Stream();                                                  //  Stream
        for (int i = 0; i < nx; i++) {
            particle.SetStress(i, 0, 0.0, 0.0);
            particle.SetStress(i, ny - 1, 0.0, stress0);
        }
        for (int j = 0; j < ny; j++) {
            particle.SetRS(0, j);
            particle.SetStress(nx - 1, j, 0.0, 0.0);
        }                                                                   //  Boundary condition
        EL::UpdateMacro(particle, rho, ux, uy, sxx, sxy, syx, syy);         //  Update macroscopic values
        for (int i = 0; i < nx*ny; i++) {
            rx[i] += ux[i];
            ry[i] += uy[i];
        }
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/elastic" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[particle.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[particle.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointVector("r", 
                [&](int _i, int _j, int _k) { return rx[particle.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return ry[particle.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j); });
        }                                                                   //  Export result per 1000 steps 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}