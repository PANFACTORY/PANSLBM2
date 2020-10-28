#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 20000, nx = 100, ny = 100, tmax = nt;
    double nu = 0.1, u0 = 0.1, rho0 = 1.0;    
    double q = 0.1, alpha0 = 1.0;
    double *alpha = new double[nx*ny];  //  Inverse permeation
    double *rho = new double[nt*nx*ny];
    double *ux = new double[nt*nx*ny];
    double *uy = new double[nt*nx*ny];  //  State variable
    double *qho = new double[nx*ny];
    double *vx = new double[nx*ny];
    double *vy = new double[nx*ny];     //  Adjoint variable

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (0.3*ny < j && j < 0.7*ny) {
            particle.SetBoundary(0, j, OTHER);
            particle.SetBoundary(nx - 1, j, OTHER);
        } else {
            particle.SetBoundary(0, j, BARRIER);
            particle.SetBoundary(nx - 1, j, BARRIER);
        }
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, BARRIER);
    }                                                   //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 0.9;
            alpha[particle.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
        }
    }

    //--------------------Direct analyze--------------------
    for (int t = 0; t < nt; t++) {
        NS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);      //  Update macroscopic values
        NS2::Collision(nu, particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);    //  Collision
        particle.Stream();                              //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                double uj = -u0*(j - 0.3*ny)*(j - 0.7*ny)/(0.04*ny*ny);
                particle.SetU(0, j, uj, 0.0);
                particle.SetRho(nx - 1, j, rho0, 0.0);
            }
        }                                               //  Boundary condition (inlet)
        NS2::ExternalForceBrinkman(particle, alpha);    //  External force by Brinkman model

        if (t > 0 && NS2::CheckConvergence(particle, 1.0e-4, &ux[nx*ny*t], &uy[nx*ny*t], &ux[nx*ny*(t - 1)], &uy[nx*ny*(t - 1)])) {
            std::cout << "\tt = " << t << "\t";
            tmax = t;
            break;
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS2::InitialCondition(i, particle, 0.0, 0.0, 0.0);
    }                                                   //  Set initial condition
    for (int t = tmax; t >= 0; t--) {
        ANS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);    //  Update macroscopic values
        ANS2::Collision(nu, particle, &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);                 //  Collision
        particle.iStream();                             //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                double uj = -u0*(j - 0.3*ny)*(j - 0.7*ny)/(0.04*ny*ny);
                particle.SetiU(0, j, uj, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                                               //  Boundary condition (inlet)
        ANS2::ExternalForce(particle, alpha, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);        //  External force by Brinkman model
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjoint.vtk", nx, ny);
    file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[nx*ny*tmax + particle.GetIndex(_i, _j)]; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) { 
        return 3.0*(vx[particle.GetIndex(_i, _j)]*ux[nx*ny*tmax + particle.GetIndex(_i, _j)] + vy[particle.GetIndex(_i, _j)]*uy[nx*ny*tmax + particle.GetIndex(_i, _j)]); 
    });
    file.AddPointScaler("q", [&](int _i, int _j, int _k) { return qho[particle.GetIndex(_i, _j)]; });
    file.AddPointVector("v", 
        [&](int _i, int _j, int _k) { return vx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return vy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    
    delete[] alpha;
    delete[] rho;
    delete[] ux;
    delete[] uy;
    delete[] qho;
    delete[] vx;
    delete[] vy;

    return 0;
}