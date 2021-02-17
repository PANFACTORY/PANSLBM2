#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/ib_navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 10000, nx = 100, ny = 100, nb = 100;
    double nu = 0.1, u0 = 0.1, v0 = 0.1, dv = 95.0/(double)nb;
    double rho[nx*ny], ux[nx*ny], uy[nx*ny], uxl[nx*ny], uyl[nx*ny], gx[nx*ny], gy[nx*ny];
    double bpx[nb], bpy[nb], bux[nb], buy[nb], bgx[nb], bgy[nb], bvx[nb], bvy[nb];
    
    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, BARRIER);
        particle.SetBoundary(nx - 1, j, BARRIER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, BARRIER);//OTHER);
    }                                           //  Set boundary condition
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particle, 1.0, 0.0, 0.0);
    }                                           //  Set initial condition

    for (int k = 0; k < nb; k++) {
        bpx[k] = 15.0*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0;
        bpy[k] = 15.0*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0;
        bvx[k] = -v0*sin(2.0*M_PI*k/(double)nb);
        bvy[k] = v0*cos(2.0*M_PI*k/(double)nb);
    }
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::UpdateMacro(particle, rho, ux, uy);       //  Update macroscopic values
    for (int t = 0; t < tmax; t++) {
        NS::Collision(nu, particle, rho, ux, uy); //  Collision
        particle.Stream();                      //  Stream
        /*for (int i = 0; i < nx; i++) {
            particle.SetU(i, ny - 1, u0, 0.0);
        }                                       //  Boundary condition (inlet)
        particle.SetU(0, ny - 1, u0, 0.0);
        particle.SetU(nx - 1, ny - 1, u0, 0.0); //  Boundary condition (corner node)*/
        particle.SmoothCorner();
        NS::UpdateMacro(particle, rho, ux, uy);   //  Update macroscopic values
        NS::ExternalForceIB(particle, ux, uy, uxl, uyl, gx, gy, nb, bpx, bpy, bux, buy, bgx, bgy, bvx, bvy, dv);
        NS::UpdateMacro(particle, rho, ux, uy);   //  Update macroscopic values

        if (t%100 == 0) {
            std::cout << t/100 << std::endl;
            VTKExport file("result/ibns" + std::to_string(t/100) + ".vtk", nx, ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particle.GetIndex(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[particle.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[particle.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j); });
        }                                       //  Export result per 1000 steps 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}