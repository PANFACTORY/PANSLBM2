#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/equation/adjointadvection.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 15000, nx = 100, ny = 100, tmax = nt - 1;
    double viscosity = 0.1, diffusivity = 0.1; 
    double u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-5, epsdq = 1.0e-5;    
    double q = 0.1, alpha0 = 1e0, beta0 = 0.1, *alpha = new double[nx*ny], *beta = new double[nx*ny];    //  Inverse permeation
    double **rho = new double*[nt], **ux = new double*[nt], **uy = new double*[nt];
    double **tem = new double*[nt], **qx = new double*[nt], **qy = new double*[nt];         //  State variable
    for (int t = 0; t < nt; t++) {
        rho[t] = new double[nx*ny]; ux[t] = new double[nx*ny];  uy[t] = new double[nx*ny];
        tem[t] = new double[nx*ny]; qx[t] = new double[nx*ny];  qy[t] = new double[nx*ny];
    }
    double *qho = new double[nx*ny], *vx = new double[nx*ny], *vy = new double[nx*ny];
    double *sem = new double[nx*ny], *rx = new double[nx*ny], *ry = new double[nx*ny];      //  Adjoint variable
    double *sensitivity = new double[nx*ny];                                                //  Sensitivity

    D2Q9<double> particlef(nx, ny);
    D2Q9<double> particleg(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (0.35*ny < j && j < 0.65*ny) {
            particlef.SetBoundary(0, j, OTHER);
            particlef.SetBoundary(nx - 1, j, OTHER);
        } else {
            particlef.SetBoundary(0, j, BARRIER);
            particlef.SetBoundary(nx - 1, j, BARRIER);
        }
        particleg.SetBoundary(0, j, OTHER);
        particleg.SetBoundary(nx - 1, j, OTHER);
    }
    for (int i = 0; i < nx; i++) {
        particlef.SetBoundary(i, 0, BARRIER);
        particlef.SetBoundary(i, ny - 1, BARRIER);
        particleg.SetBoundary(i, 0, OTHER);
        particleg.SetBoundary(i, ny - 1, OTHER);
    }                                                   //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 1.0;
            alpha[particlef.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
            beta[particleg.GetIndex(i, j)] = beta0*q*(1.0 - gamma)/(gamma + q);;
        }
    }

    //--------------------Direct analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 1.0, 0.0, 0.0);
    }                                                                       //  Set initial condition
    for (int t = 0; t < nt; t++) {
        NS::UpdateMacro(particlef, rho[t], ux[t], uy[t]);
        AD::UpdateMacro(particleg, tem[t], qx[t], qy[t], ux[t], uy[t]);     //  Update macroscopic values
        NS::Collision(viscosity, particlef, rho[t], ux[t], uy[t]);
        AD::Collision(diffusivity, particleg, tem[t], ux[t], uy[t]);        //  Collision
        particlef.Stream();
        particleg.Stream();                                                 //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.35*ny < j && j < 0.65*ny) {
                double uj = -u0*(j - 0.35*ny)*(j - 0.65*ny)/(0.0225*ny*ny);
                particlef.SetU(0, j, uj, 0.0);
                particlef.SetRho(nx - 1, j, rho0, 0.0);
                particleg.SetTemperature(0, j, tem0);
                int ij = particlef.GetIndex(nx - 1, j);
                particleg.SetFlux(nx - 1, j, ux[t][ij], uy[t][ij], q0);
            } else {
                particleg.SetFlux(0, j, 0.0, 0.0, q0);
                particleg.SetFlux(nx - 1, j, 0.0, 0.0, q0);
            }
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetFlux(i, 0, 0.0, 0.0, q0);
            particleg.SetFlux(i, ny - 1, 0.0, 0.0, q0);
        }                                                                   //  Boundary condition (inlet)
        NS::ExternalForceBrinkman(particlef, alpha, alpha);
        AD::ExternalForceHeatgeneration(particleg, beta);                   //  External force term
        if (t > 0 && NS::CheckConvergence(particlef, epsdu, ux[t], uy[t], ux[t - 1], uy[t - 1])  && NS::CheckConvergence(particleg, epsdq, qx[t], qy[t], qx[t - 1], qy[t - 1])) {
            std::cout << "\tt = " << t << "\t";
            tmax = t;
            break;
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 0.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 0.0, 0.0, 0.0);
    }                                                                       //  Set initial condition
    for (int t = tmax; t >= 0; t--) {
        ANS::UpdateMacro(particlef, rho[t], ux[t], uy[t], qho, vx, vy);
        AAD::UpdateMacro(particleg, sem, rx, ry);                           //  Update macroscopic values
        ANS::Collision(viscosity, particlef, ux[t], uy[t], qho, vx, vy);
        AAD::Collision(diffusivity, particleg, ux[t], uy[t], sem, rx, ry);  //  Collision
        particlef.iStream();
        particleg.iStream();                                                //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.35*ny < j && j < 0.65*ny) {
                double uj = -u0*(j - 0.35*ny)*(j - 0.65*ny)/(0.0225*ny*ny);
                particlef.SetiU(0, j, uj, 0.0);
                particleg.SetiTemperature(0, j);
                int ij = particlef.GetIndex(nx - 1, j);
                particlef.SetiRhoFlux(particleg, nx - 1, j, rho[t][ij], ux[t][ij], uy[t][ij], tem[t][ij]);
            } else {
                particleg.SetiFlux(0, j, 0.0, 0.0);
            }
            particleg.SetiFlux(nx - 1, j, 0.0, 0.0);
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetiFlux(i, 0, 0.0, 0.0);
            particleg.SetiFlux(i, ny - 1, 0.0, 0.0);
        }                                                                   //  Boundary condition (inlet)
        ANS::ExternalForceHeatexchange(diffusivity, particlef, particleg, rho[t], ux[t], uy[t], tem[t]);
        ANS::ExternalForceBrinkman(particlef, alpha, alpha, rho[t], ux[t], uy[t]);
        AAD::ExternalForceHeatexchange(particleg, sem, beta);               //  External force by Brinkman model
    }

    //--------------------Analyse sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 1.0;
            int ij = particlef.GetIndex(i, j);
            sensitivity[ij] = 3.0*(vx[ij]*ux[tmax][ij] + vy[ij]*uy[tmax][ij])*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0)) - (1.0 - tem[tmax][ij])*(1.0 + sem[ij])*(-beta0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[ij])) {
                sensitivitymax = fabs(sensitivity[ij]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointadvection.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[tmax][particlef.GetIndex(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[tmax][particlef.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[tmax][particlef.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[tmax][particleg.GetIndex(_i, _j)]; });
    file.AddPointVector("q", 
        [&](int _i, int _j, int _k) { return qx[tmax][particleg.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qy[tmax][particleg.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particlef.GetIndex(_i, _j)];  });
    file.AddPointScaler("qho", [&](int _i, int _j, int _k) { return qho[particlef.GetIndex(_i, _j)]; });
    file.AddPointVector("v", 
        [&](int _i, int _j, int _k) { return vx[particlef.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return vy[particlef.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("S", [&](int _i, int _j, int _k) { return sem[particleg.GetIndex(_i, _j)]; });
    file.AddPointVector("r", 
        [&](int _i, int _j, int _k) { return rx[particleg.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return ry[particleg.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    
    delete[] alpha; delete[] beta;
    for (int t = 0; t < nt; t++) {
        delete[] rho[t];    delete[] ux[t]; delete[] uy[t]; delete[] tem[t];    delete[] qx[t]; delete[] qy[t];
    }
    delete[] rho;   delete[] ux;    delete[] uy;    delete[] tem;   delete[] qx;    delete[] qy;
    delete[] qho;   delete[] vx;    delete[] vy;    delete[] sem;   delete[] rx;    delete[] ry;
    delete[] sensitivity;

    return 0;
}