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
    int nt = 9500, nx = 200, ny = 100, tmax = nt - 1;
    double nu = 0.1, u0 = 0.0109, rho0 = 1.0, epsdu = 1.0e-5;    
    double q = 0.1, alpha0 = 1e0, *alpha = new double[nx*ny];                           //  Inverse permeation
    double **rho = new double*[nt], **ux = new double*[nt], **uy = new double*[nt];     //  State variable
    for (int t = 0; t < nt; t++) {
        rho[t] = new double[nx*ny];
        ux[t] = new double[nx*ny];
        uy[t] = new double[nx*ny];
    }
    double *qho = new double[nx*ny], *vx = new double[nx*ny], *vy = new double[nx*ny];  //  Adjoint variable
    double *sensitivity = new double[nx*ny];
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 0.0;
    }                                                                                   //  Sensitivity

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (j < 0.33*ny) {
            particle.SetBoundary(0, j, OTHER);
            particle.SetBoundary(nx - 1, j, OTHER);
        } else {
            particle.SetBoundary(0, j, BARRIER);
            particle.SetBoundary(nx - 1, j, BARRIER);
        }
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, MIRROR);
        particle.SetBoundary(i, ny - 1, BARRIER);
    }                                                       //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
            alpha[particle.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
        }
    }

    //--------------------Direct analyze--------------------
    for (int t = 0; t < nt; t++) {
        NS::UpdateMacro(particle, rho[t], ux[t], uy[t]);    //  Update macroscopic values
        NS::Collision(nu, particle, rho[t], ux[t], uy[t]);  //  Collision
        particle.Stream();                                  //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particle.SetU(0, j, uj, 0.0);
                particle.SetRho(nx - 1, j, rho0, 0.0);
            }
        }                                                   //  Boundary condition (inlet)
        NS::ExternalForceBrinkman(particle, alpha, alpha);  //  External force by Brinkman model

        if (t > 0 && NS::CheckConvergence(particle, epsdu, ux[t], uy[t], ux[t - 1], uy[t - 1])) {
            std::cout << "\tt = " << t << "\t";
            tmax = t;
            break;
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particle, 0.0, 0.0, 0.0);
    }                                                       //  Set initial condition
    for (int t = tmax; t >= 0; t--) {
        ANS::UpdateMacro(particle, rho[t], ux[t], uy[t], qho, vx, vy);  //  Update macroscopic values
        ANS::Collision(nu, particle, ux[t], uy[t], qho, vx, vy);        //  Collision
        particle.iStream();                                 //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particle.SetiUPressureDrop(0, j, uj, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                                                   //  Boundary condition (inlet)
        ANS::ExternalForceBrinkman(particle, alpha, alpha, rho[t], ux[t], uy[t]);   //  External force by Brinkman model
        ANS::UpdateSensitivity(particle, ux[t], uy[t], sensitivity);    //  Update sensitivity
    }
    ANS::UpdateMacro(particle, rho[0], ux[0], uy[0], qho, vx, vy);      //  Update macroscopic values
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
            int ij = particle.GetIndex(i, j);
            sensitivity[ij] = 3.0*(vx[ij]*ux[tmax][ij] + vy[ij]*uy[tmax][ij])*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[ij])) {
                sensitivitymax = fabs(sensitivity[ij]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjoint.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[tmax][particle.GetIndex(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[tmax][particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[tmax][particle.GetIndex(_i, _j)]; },
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
    for (int t = 0; t < nt; t++) {
        delete[] rho[t];
        delete[] ux[t];
        delete[] uy[t];
    }
    delete[] rho;   delete[] ux;    delete[] uy;
    delete[] qho;   delete[] vx;    delete[] vy;
    delete[] sensitivity;

    return 0;
}