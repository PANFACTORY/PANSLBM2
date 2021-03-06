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
    int nt = 30000, nx = 100, ny = 50, tmax = nt - 1;
    double viscosity = 0.1, diffusivity = 0.1/6.0; 
    double u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-8, epsdq = 1.0e-8;    
    double q = 0.1, alpha0 = 100.0/(double)nx, beta0 = 0.1/(double)nx;
    double *alpha = new double[nx*ny], *beta = new double[nx*ny];                       //  Inverse permeation
    double *rho = new double[nx*ny], *ux = new double[nx*ny], *uy = new double[nx*ny];
    double *tem = new double[nx*ny], *qx = new double[nx*ny], *qy = new double[nx*ny];  //  State variable
    double *sensitivity = new double[nx*ny], *gamma = new double[nx*ny];                //  Sensitivity
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 0.0;
    }
    double dgamma = -1.0e-5;  

    D2Q9<double> particlef(nx, ny);
    D2Q9<double> particleg(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (j < 0.33*ny) {
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
        particlef.SetBoundary(i, 0, MIRROR);
        particlef.SetBoundary(i, ny - 1, BARRIER);
        particleg.SetBoundary(i, 0, MIRROR);
        particleg.SetBoundary(i, ny - 1, OTHER);
    }                                                               //  Set boundary condition

    //--------------------Function to get objective--------------------
    auto getObjective = [&](double *_gamma) {
        //--------------------Set design variable--------------------
        for (int i = 0; i < nx*ny; i++) {
            alpha[i] = alpha0*q*(1.0 - _gamma[i])/(_gamma[i] + q);
            beta[i] = beta0*q*(1.0 - _gamma[i])/(_gamma[i] + q);
        }

        //--------------------Direct analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
            AD::InitialCondition(i, particleg, 1.0, 0.0, 0.0);
        }                                                           //  Set initial condition
        for (int t = 0; t < nt; t++) {
            NS::UpdateMacro(particlef, rho, ux, uy);
            AD::UpdateMacro(particleg, tem, qx, qy, ux, uy);        //  Update macroscopic values
            NS::Collision(viscosity, particlef, rho, ux, uy);
            AD::Collision(diffusivity, particleg, tem, ux, uy);     //  Collision
            particlef.Stream();
            particleg.Stream();                                     //  Stream
            for (int j = 0; j < ny; j++) {
                if (j < 0.33*ny) {
                    double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                    particlef.SetU(0, j, uj, 0.0);
                    particlef.SetRho(nx - 1, j, rho0, 0.0);
                    particleg.SetTemperature(0, j, tem0);
                    int ij = particlef.GetIndex(nx - 1, j);
                    particleg.SetFlux(nx - 1, j, ux[ij], uy[ij], q0);
                } else {
                    particleg.SetFlux(0, j, 0.0, 0.0, q0);
                    particleg.SetFlux(nx - 1, j, 0.0, 0.0, q0);
                }
            }
            for (int i = 0; i < nx; i++) {
                particleg.SetFlux(i, ny - 1, 0.0, 0.0, q0);
            }                                                       //  Boundary condition (inlet)
            NS::ExternalForceBrinkman(particlef, alpha, alpha);
            AD::ExternalForceHeatgeneration(particleg, beta);       //  External force term
        }

        //--------------------Get objective--------------------
        double objective = 0.0;
        for (int i = 0; i < nx*ny; i++) {
            objective += -beta[i]*(1.0 - tem[i]);
        }
        return objective; 
    };

    //--------------------Get objective of neutral--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            gamma[particlef.GetIndex(i, j)] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
        }
    }
    double objective_base = getObjective(gamma);
    std::cout << objective_base << std::endl;

    //--------------------FDM loop--------------------
    for (int l = 0; l < nx; l++) {
        std::cout << l << "\t";
        int k = particlef.GetIndex(l, 5);
        
        //--------------------Set design variable--------------------
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                gamma[particlef.GetIndex(i, j)] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
            }
        }
        gamma[k] += dgamma;

        //--------------------Get objective--------------------
        double objective = getObjective(gamma);
        sensitivity[k] = (objective - objective_base)/dgamma;
        std::cout << objective << "\t" << sensitivity[k] << std::endl;
    }

    //--------------------Normalize sensitivity--------------------
    double sensitivity_max = 0.0;
    for (int i = 0; i < nx*ny; i++) {
        if (sensitivity_max < fabs(sensitivity[i])) {
            sensitivity_max = fabs(sensitivity[i]);
        }
    }
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] /= sensitivity_max;
    }

    //--------------------Export result--------------------
    VTKExport file("result/FDMadjointadvectionY5.vtk", nx, ny);
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particlef.GetIndex(_i, _j)];  });
    
    delete[] alpha; delete[] beta;
    delete[] rho;   delete[] ux;    delete[] uy;    delete[] tem;   delete[] qx;    delete[] qy;
    delete[] sensitivity;           delete[] gamma;

    return 0;
}