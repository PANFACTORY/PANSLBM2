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
    double nu = 0.1, u0 = 0.0218, rho0 = 1.0;    
    double q = 0.1, alpha0 = 4e2, *alpha = new double[nx*ny];                                   //  Inverse permeation
    double *rho = new double[nt*nx*ny], *ux = new double[nt*nx*ny], *uy = new double[nt*nx*ny]; //  State variable
    double *qho = new double[nx*ny], *vx = new double[nx*ny], *vy = new double[nx*ny];          //  Adjoint variable
    double *sensitivity = new double[nx*ny];
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 0.0;
    }                                                                                           //  Sensitivity

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (0.35*ny < j && j < 0.65*ny) {
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
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 1.0;
            alpha[particle.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
        }
    }

    //--------------------Direct analyze--------------------
    for (int t = 0; t < nt; t++) {
        NS::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);       //  Update macroscopic values
        NS::Collision(nu, particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);     //  Collision
        particle.Stream();                              //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.35*ny < j && j < 0.65*ny) {
                double uj = -u0*(j - 0.35*ny)*(j - 0.65*ny)/(0.0225*ny*ny);
                particle.SetU(0, j, uj, 0.0);
                particle.SetRho(nx - 1, j, rho0, 0.0);
            }
        }                                               //  Boundary condition (inlet)
        NS::ExternalForceBrinkman(particle, alpha, alpha);     //  External force by Brinkman model

        if (t > 0 && NS::CheckConvergence(particle, 1.0e-4, &ux[nx*ny*t], &uy[nx*ny*t], &ux[nx*ny*(t - 1)], &uy[nx*ny*(t - 1)])) {
            std::cout << "\tt = " << t << "\t";
            tmax = t;
            break;
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particle, 0.0, 0.0, 0.0);
    }                                                   //  Set initial condition
    for (int t = tmax; t >= 0; t--) {
        ANS::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);     //  Update macroscopic values
        ANS::Collision(nu, particle, &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);                  //  Collision
        particle.iStream();                             //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.35*ny < j && j < 0.65*ny) {
                double uj = -u0*(j - 0.35*ny)*(j - 0.65*ny)/(0.0225*ny*ny);
                particle.SetiU(0, j, uj, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                                               //  Boundary condition (inlet)
        ANS::ExternalForceBrinkman(particle, alpha, alpha, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]); //  External force by Brinkman model
        ANS::UpdateSensitivity(particle, &ux[nx*ny*t], &uy[nx*ny*t], sensitivity);              //  Update sensitivity
    }
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 1.0;
            sensitivity[particle.GetIndex(i, j)] *= 3.0*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[particle.GetIndex(i, j)])) {
                sensitivitymax = fabs(sensitivity[particle.GetIndex(i, j)]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjoint.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[nx*ny*tmax + particle.GetIndex(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particle.GetIndex(_i, _j)];  });
    file.AddPointScaler("q", [&](int _i, int _j, int _k) { return qho[particle.GetIndex(_i, _j)]; });
    file.AddPointVector("v", 
        [&](int _i, int _j, int _k) { return vx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return vy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    
    delete[] alpha;
    delete[] rho;   delete[] ux;    delete[] uy;
    delete[] qho;   delete[] vx;    delete[] vy;
    delete[] sensitivity;

    return 0;
}