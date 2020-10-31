#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d3q15.h"
#include "../src/equation/navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 30, ny = 30, nz = 30;
    double nu = 0.1, u0 = 0.1, theta = 30.0;
    double *rho = new double[nx*ny*nz], *u = new double[nx*ny*nz], *v = new double[nx*ny*nz], *w = new double[nx*ny*nz];
    
    D3Q15<double> particle(nx, ny, nz);
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            particle.SetBoundary(0, j, k, BARRIER);
            particle.SetBoundary(nx - 1, j, k, BARRIER);
        } 
    }                                               //  Set boundary condition of x
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            particle.SetBoundary(i, 0, k, BARRIER);
            particle.SetBoundary(i, ny - 1, k, BARRIER);
        } 
    }                                               //  Set boundary condition ok y
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            particle.SetBoundary(i, j, 0, BARRIER);
            particle.SetBoundary(i, j, nz - 1, OTHER);
        } 
    }                                               //  Set boundary condition ok z
    for (int i = 0; i < nx*ny*nz; i++) {
        NS::InitialCondition(i, particle, 1.0, 0.0, 0.0, 0.0);
    }                                               //  Set initial condition
    
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        NS::UpdateMacro(particle, rho, u, v, w);    //  Update macroscopic values
        NS::Collision(nu, particle, rho, u, v, w);  //  Collision
        particle.Stream();                          //  Stream
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                particle.SetU(i, j, nz - 1, u0*cos(theta*M_PI/180.0), u0*sin(theta*M_PI/180.0), 0.0);
            } 
        }                                           //  Boundary condition (inlet)
      
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/ns3D" + std::to_string(t/1000) + ".vtk", nx, ny, nz);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particle.GetIndex(_i, _j, _k)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return u[particle.GetIndex(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return v[particle.GetIndex(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return w[particle.GetIndex(_i, _j, _k)];; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j, _k); });
        }                                           //  Export result per 1000 steps 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    delete[] rho;
    delete[] u;
    delete[] v;
    delete[] w;
}