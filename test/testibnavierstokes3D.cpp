#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <chrono>

#include "../src/particle/d3q15.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/ib_navierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 1000, nx = 300, ny = 30, nz = 30, ntheta = 32, nphi = 30, nb = ntheta*nphi;
    double nu = 0.1, u0 = 0.1, theta = 0.0, dv = 1.0;
    double *rho = new double[nx*ny*nz], *ux = new double[nx*ny*nz], *uy = new double[nx*ny*nz], *uz = new double[nx*ny*nz];
    double *uxl = new double[nx*ny*nz], *uyl = new double[nx*ny*nz], *uzl = new double[nx*ny*nz], *gx = new double[nx*ny*nz], *gy = new double[nx*ny*nz], *gz = new double[nx*ny*nz];
    double *bpx = new double[nb], *bpy = new double[nb], *bpz = new double[nb], *bvx = new double[nb], *bvy = new double[nb], *bvz = new double[nb];
    double *bux = new double[nb], *buy = new double[nb], *buz = new double[nb], *bgx = new double[nb], *bgy = new double[nb], *bgz = new double[nb];

    D3Q15<double> particle(nx, ny, nz);
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            particle.SetBoundary(0, j, k, OTHER);
            particle.SetBoundary(nx - 1, j, k, OUTLET);
        } 
    }                                                   //  Set boundary condition of x
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            particle.SetBoundary(i, 0, k, PERIODIC);
            particle.SetBoundary(i, ny - 1, k, PERIODIC);
        } 
    }                                                   //  Set boundary condition ok y
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            particle.SetBoundary(i, j, 0, PERIODIC);
            particle.SetBoundary(i, j, nz - 1, PERIODIC);
        } 
    }                                                   //  Set boundary condition ok z
    for (int i = 0; i < nx*ny*nz; i++) {
        NS::InitialCondition(i, particle, 1.0, 0.0, 0.0, 0.0);
    }                                                   //  Set initial condition
    
    for (int i = 0; i < ntheta; i++) {
        for (int j = 0; j < nphi; j++) {
            int ij = i + ntheta*j;
            bpx[ij] = 5.0*cos(2.0*M_PI*i/(double)ntheta) + (double)nx/4.0;
            bpy[ij] = 5.0*sin(2.0*M_PI*i/(double)ntheta) + (double)ny/2.0;
            bpz[ij] = (double)j;
            bvx[ij] = 0.0;
            bvy[ij] = 0.0;
            bvz[ij] = 0.0;
        }
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::UpdateMacro(particle, rho, ux, uy, uz);         //  Update macroscopic values
    for (int t = 0; t < tmax; t++) {
        NS::Collision(nu, particle, rho, ux, uy, uz);   //  Collision
        particle.Stream();                              //  Stream
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                particle.SetU(0, j, k, u0*cos(theta*M_PI/180.0), u0*sin(theta*M_PI/180.0), 0.0);
            } 
        }                                               //  Boundary condition (inlet)
        particle.SmoothCorner();                        //  Smooth corner nodes
        NS::UpdateMacro(particle, rho, ux, uy, uz);     //  Update macroscopic values
        NS::ExternalForceIB(particle, ux, uy, uz, uxl, uyl, uzl, gx, gy, gz, nb, bpx, bpy, bpz, bux, buy, buz, bgx, bgy, bgz, bvx, bvy, bvz, dv);
        NS::UpdateMacro(particle, rho, ux, uy, uz);     //  Update macroscopic values

        if (t%10 == 0) {
            std::cout << t/10 << std::endl;
            VTKExport file("result/ibns3D_" + std::to_string(t/10) + ".vtk", nx, ny, nz);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particle.GetIndex(_i, _j, _k)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[particle.GetIndex(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return uy[particle.GetIndex(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return uz[particle.GetIndex(_i, _j, _k)];; }
            );
            file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j, _k); });
        }                                               //  Export result per 1000 steps 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

    delete[] rho;   delete[] ux;    delete[] uy;    delete[] uz;
    delete[] uxl;   delete[] uyl;   delete[] uzl;   delete[] gx;    delete[] gy;    delete[] gz;
    delete[] bpx;   delete[] bpy;   delete[] bpz;   delete[] bvx;   delete[] bvy;   delete[] bvz;
    delete[] bux;   delete[] buy;   delete[] buz;   delete[] bgx;   delete[] bgy;   delete[] bgz;
}