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
    int tmax = 30000, nx = 100, ny = 100, nk = 20;
    double nu = 0.1, v0 = 0.0, K = 0.0, epsu = 1.0e-6, radius = 30.0, Remax = 100.0, Remin = 0.1;
    double L = 2.0*radius, dpmax = pow(nu*Remax/L, 2.0), dpmin = pow(nu*Remin/L, 2.0), surface = 2.0*M_PI*radius, alpha = 2.0*log(Remin/Remax)/(double)(nk - 1);
    double rho[nx*ny], ux[nx*ny], uy[nx*ny], tmpux[nx*ny] = { 0 }, tmpuy[nx*ny] = { 0 }, uxl[nx*ny], uyl[nx*ny], gx[nx*ny], gy[nx*ny];
    const int nb = (int)surface;
    double dv = surface/(double)nb;
    double bpx[nb], bpy[nb], bux[nb], buy[nb], bgx[nb], bgy[nb], bvx[nb], bvy[nb];

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, OTHER);
        particle.SetBoundary(nx - 1, j, OTHER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, PERIODIC);
        particle.SetBoundary(i, ny - 1, PERIODIC);
    }                                               //  Set boundary condition

    for (int k = 0; k < nb; k++) {
        bpx[k] = radius*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0;
        bpy[k] = radius*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0;
        bvx[k] = -v0*sin(2.0*M_PI*k/(double)nb);
        bvy[k] = v0*cos(2.0*M_PI*k/(double)nb);
    }
    
    for (int k = 0; k < nk; k++) {
        double dp = dpmax*exp(alpha*k);
        std::cout << "Re = " << sqrt(dp)*L/nu;

        std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

        //--------------------Direct analyze loop--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS::InitialCondition(i, particle, 1.0, 0.0, 0.0);
        }                                               //  Set initial condition
        NS::UpdateMacro(particle, rho, ux, uy);         //  Update macroscopic values
        for (int t = 0; t < tmax; t++) {
            //NS::CollisionMRT(1.19, 1.4, 1.2, 1.2, 1.0/(3.0*nu + 0.5), 1.0/(3.0*nu + 0.5), particle, rho, ux, uy);     //  Collision
            NS::Collision(nu, particle, rho, ux, uy);   //  Collision
            particle.Stream();                          //  Stream
            for (int j = 0; j < ny; j++) {
                particle.SetDP(0, j, nx - 1, j, dp);
            }                                           //  Boundary condition (inlet)
            particle.SmoothCorner();
            NS::UpdateMacro(particle, rho, ux, uy);     //  Update macroscopic values
            NS::ExternalForceIB(particle, ux, uy, uxl, uyl, gx, gy, nb, bpx, bpy, bux, buy, bgx, bgy, bvx, bvy, dv);
            NS::UpdateMacro(particle, rho, ux, uy);     //  Update macroscopic values

            double sumu = 0.0, sumdu = 0.0;
            for (int i = 0; i < particle.np; i++) {
                sumu += ux[i]*ux[i] + uy[i]*uy[i];
                sumdu += (ux[i] - tmpux[i])*(ux[i] - tmpux[i]) + (uy[i] - tmpuy[i])*(uy[i] - tmpuy[i]);
                tmpux[i] = ux[i];
                tmpuy[i] = uy[i];
            }
            if (t > 0 && sqrt(sumdu/sumu) < epsu) {
                std::cout << "\tConvergence at t = " << t;
                break;
            }
        }

        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        std::cout << "\tTime cost : " << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

        //--------------------Get permeation--------------------
        double sumv = 0.0;
        for (int j = 0; j < ny; j++) {
            sumv += sqrt(pow(ux[particle.GetIndex(0, j)], 2.0) + pow(uy[particle.GetIndex(0, j)], 2.0));
        }
        sumv /= (double)ny;

        std::cout << "\tv mean = " << sumv << "\tdp : " << dp << std::endl;

        //--------------------Export result--------------------
        VTKExport file("result/ibns" + std::to_string(k) + ".vtk", nx, ny);
        file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[particle.GetIndex(_i, _j)]; });
        file.AddPointVector("u", 
            [&](int _i, int _j, int _k) { return ux[particle.GetIndex(_i, _j)]; },
            [&](int _i, int _j, int _k) { return uy[particle.GetIndex(_i, _j)]; },
            [](int _i, int _j, int _k) { return 0.0; }
        );
        file.AddPointScaler("boundary", [&](int _i, int _j, int _k) { return particle.GetBoundary(_i, _j); });
    }

    VTKExportIB model("result/model.vtk");
    model.AddPoint(nb, bpx, bpy);
}