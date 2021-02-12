#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation_avx/navierstokes_avx.h"
#include "../src/equation_avx/adjointnavierstokes_avx.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 10000, nx = 200, ny = 100, tmax = nt - 1;
    double nu = 0.1, u0 = 0.0109, rho0 = 1.0, epsdu = 1.0e-15;    
    double q = 0.01, alpha0 = 50.0/(double)nx, *alpha = new double[nx*ny];              //  Inverse permeation
    double *rhot = new double[nx*ny], *uxt = new double[nx*ny], *uyt = new double[nx*ny], *uxtm1 = new double[nx*ny], *uytm1 = new double[nx*ny];   //  State variable
    double *irho = new double[nx*ny], *iux = new double[nx*ny], *iuy = new double[nx*ny], *imx = new double[nx*ny], *imy = new double[nx*ny];       //  Adjoint variable
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
    }                                                                                   //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
            alpha[particle.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
        }
    }

    //--------------------Direct analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particle, 1.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    NS::UpdateMacro(particle, rhot, uxt, uyt);                                          //  Update macroscopic values
    for (int t = 0; t < nt - 1; t++) {
        NS::Collision(nu, particle, rhot, uxt, uyt);                                    //  Collision
        particle.Stream();                                                              //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particle.SetU(0, j, uj, 0.0);
                particle.SetRho(nx - 1, j, rho0, 0.0);
            }
        }                                                                               //  Boundary condition (inlet)
        NS::UpdateMacro(particle, rhot, uxt, uyt);                                      //  Update macroscopic values
        NS::ExternalForceBrinkman(particle, rhot, uxt, uyt, alpha, alpha);              //  External force by Brinkman model
        if (t > 0 && NS::CheckConvergence(particle, epsdu, uxt, uyt, uxtm1, uytm1)) {
            std::cout << "\tt = " << t << "\t";
            tmax = t;
            break;
        }
        for (int i = 0; i < nx*ny; i++) {
            uxtm1[i] = uxt[i];
            uytm1[i] = uyt[i];
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        ANS::InitialCondition(i, particle, 0.0, 0.0, 0.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    ANS::UpdateMacro(particle, rhot, uxt, uyt, irho, iux, iuy, imx, imy);               //  Update macroscopic values
    for (int t = tmax; t >= 0; t--) {
        ANS::Collision(nu, particle, uxt, uyt, irho, iux, iuy);                         //  Collision
        particle.iStream();                                                             //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particle.SetiUPressureDrop(0, j, uj, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                                                                               //  Boundary condition (inlet)
        ANS::UpdateMacro(particle, rhot, uxt, uyt, irho, iux, iuy, imx, imy);           //  Update macroscopic values
        ANS::ExternalForceBrinkman(particle, rhot, iux, iuy, imx, imy, alpha, alpha);   //  External force by Brinkman model
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
            int ij = particle.GetIndex(i, j);
            sensitivity[ij] = 3.0*(imx[ij]*uxt[ij] + imy[ij]*uyt[ij])*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[ij])) {
                sensitivitymax = fabs(sensitivity[ij]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjoint_AL50_q1e-2_avx.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rhot[particle.GetIndex(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return uxt[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uyt[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particle.GetIndex(_i, _j)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[particle.GetIndex(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[particle.GetIndex(_i, _j)]; });
    
    delete[] alpha;
    delete[] rhot;  delete[] uxt;   delete[] uyt;   delete[] uxtm1; delete[] uytm1;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] imx;   delete[] imy;
    delete[] sensitivity;

    return 0;
}