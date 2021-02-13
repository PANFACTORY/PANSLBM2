#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation_avx/navierstokes_avx.h"
#include "../src/equation_avx/advection_avx.h"
#include "../src/equation_avx/adjointnavierstokes_avx.h"
#include "../src/equation_avx/adjointadvection_avx.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 30000, nx = 100, ny = 50, tmax = nt - 1;
    double viscosity = 0.1, diffusivity = 0.1/6.0; 
    double u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-8, epsdq = 1.0e-8;    
    double q = 0.01, alpha0 = 50.0/(double)nx, beta0 = 0.1/(double)nx, *alpha = new double[nx*ny], *beta = new double[nx*ny];                       //  Inverse permeation
    double *rhot = new double[nx*ny], *uxt = new double[nx*ny], *uyt = new double[nx*ny], *uxtm1 = new double[nx*ny], *uytm1 = new double[nx*ny];
    double *temt = new double[nx*ny], *qxt = new double[nx*ny], *qyt = new double[nx*ny], *qxtm1 = new double[nx*ny], *qytm1 = new double[nx*ny];   //  State variable
    double *irho = new double[nx*ny], *iux = new double[nx*ny], *iuy = new double[nx*ny], *imx = new double[nx*ny], *imy = new double[nx*ny];
    double *item = new double[nx*ny], *iqx = new double[nx*ny], *iqy = new double[nx*ny];                                                           //  Adjoint variable
    double *sensitivity = new double[nx*ny];                                                                                                        //  Sensitivity

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
    }                                                                                   //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
            alpha[particlef.GetIndex(i, j)] = alpha0*q*(1.0 - gamma)/(gamma + q);
            beta[particleg.GetIndex(i, j)] = beta0*q*(1.0 - gamma)/(gamma + q);;
        }
    }

    //--------------------Direct analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 1.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    NS::UpdateMacro(particlef, rhot, uxt, uyt);
    AD::UpdateMacro(particleg, temt, qxt, qyt, uxt, uyt);                               //  Update macroscopic values
    for (int t = 0; t < nt - 1; t++) {
        //  Navier-Stokes
        NS::Collision(viscosity, particlef, rhot, uxt, uyt);                            //  Collision
        particlef.Stream();                                                             //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particlef.SetU(0, j, uj, 0.0);
                particlef.SetRho(nx - 1, j, rho0, 0.0);   
            }
        }                                                                               //  Boundary condition (inlet)
        NS::UpdateMacro(particlef, rhot, uxt, uyt);                                     //  Update macroscopic values
        NS::ExternalForceBrinkman(particlef, rhot, uxt, uyt, alpha, alpha);             //  External force term
        
        //  Advection
        AD::Collision(diffusivity, particleg, temt, uxt, uyt);                          //  Collision
        particleg.Stream();                                                             //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particleg.SetTemperature(0, j, uj, 0.0, tem0);
                int ij = particlef.GetIndex(nx - 1, j);
                particleg.SetFlux(nx - 1, j, uxt[ij], uyt[ij], q0);
            } else {
                particleg.SetFlux(0, j, 0.0, 0.0, q0);
                particleg.SetFlux(nx - 1, j, 0.0, 0.0, q0);
            }
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetFlux(i, ny - 1, 0.0, 0.0, q0);
        }                                                                               //  Boundary condition (inlet)
        AD::UpdateMacro(particleg, temt, qxt, qyt, uxt, uyt);                           //  Update macroscopic values
        AD::ExternalForceHeatgeneration(particleg, temt, beta);                         //  External force term
        
        //  Check convergence
        if (t > 0 && NS::CheckConvergence(particlef, epsdu, uxt, uyt, uxtm1, uytm1) && NS::CheckConvergence(particleg, epsdq, qxt, qyt, qxtm1, qytm1)) {
            std::cout << "\tt = " << t << "\t";
            tmax = t + 1;
            break;
        }
        for (int i = 0; i < nx*ny; i++) {
            uxtm1[i] = uxt[i];
            uytm1[i] = uyt[i];
            qxtm1[i] = qxt[i];
            qytm1[i] = qyt[i];
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        ANS::InitialCondition(i, particlef, 0.0, 0.0, 0.0, 0.0, 0.0);
        AAD::InitialCondition(i, particleg, 0.0, 0.0, 0.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    AAD::UpdateMacro(particleg, item, iqx, iqy);
    ANS::UpdateMacro(particlef, rhot, uxt, uyt, irho, iux, iuy, imx, imy);              //  Update macroscopic values
    for (int t = tmax; t >= 0; t--) {
        //  Adjoint advection
        AAD::Collision(diffusivity, particleg, uxt, uyt, item, iqx, iqy);               //  Collision
        particleg.iStream();                                                            //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particleg.SetiTemperature(0, j, uj, 0.0);
                int ij = particleg.GetIndex(nx - 1, j);
                particleg.SetiFlux(nx - 1, j, uxt[ij], uyt[ij]);
            } else {
                particleg.SetiFlux(0, j, 0.0, 0.0);
                particleg.SetiFlux(nx - 1, j, 0.0, 0.0);
            }
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetiFlux(i, ny - 1, 0.0, 0.0);
        }                                                                               //  Boundary condition (inlet)
        AAD::UpdateMacro(particleg, item, iqx, iqy);                                    //  Update macroscopic values
        AAD::ExternalForceHeatexchange(particleg, item, beta);                          //  External force by Brinkman model
        
        //  Adjoint Navier-Stokes
        ANS::Collision(viscosity, particlef, uxt, uyt, irho, iux, iuy);                 //  Collision
        particlef.iStream();                                                            //  Stream
        ANS::ExternalForceHeatexchange(diffusivity, particlef, rhot, uxt, uyt, temt, iqx, iqy);
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particlef.SetiU(0, j, uj, 0.0);
                int ij = particlef.GetIndex(nx - 1, j);
                particlef.SetiRhoFlux(particleg, nx - 1, j, rho0, uxt[ij], uyt[ij], temt[ij]);
            }
        }                                                                               //  Boundary condition (inlet)
        ANS::UpdateMacro(particlef, rhot, uxt, uyt, irho, iux, iuy, imx, imy);          //  Update macroscopic values 
        ANS::ExternalForceBrinkman(particlef, rhot, iux, iuy, imx, imy, alpha, alpha);  //  External force by Brinkman model
    }

    //--------------------Analyse sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
            int ij = particlef.GetIndex(i, j);
            sensitivity[ij] = 3.0*(imx[ij]*uxt[ij] + imy[ij]*uyt[ij])*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0)) - (1.0 - temt[ij])*(1.0 + item[ij])*(-beta0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[ij])) {
                sensitivitymax = fabs(sensitivity[ij]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointadvection_AL50_BE1e-1_q1e-2_avx.vtk", nx, ny);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rhot[particlef.GetIndex(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return uxt[particlef.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uyt[particlef.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("T", [&](int _i, int _j, int _k) { return temt[particleg.GetIndex(_i, _j)]; });
    file.AddPointVector("q", 
        [&](int _i, int _j, int _k) { return qxt[particleg.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return qyt[particleg.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particlef.GetIndex(_i, _j)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[particlef.GetIndex(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[particlef.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[particlef.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointVector("im", 
                [&](int _i, int _j, int _k) { return imx[particlef.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return imy[particlef.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
    file.AddPointScaler("iT", [&](int _i, int _j, int _k) { return item[particleg.GetIndex(_i, _j)]; });
    file.AddPointVector("iq", 
        [&](int _i, int _j, int _k) { return iqx[particleg.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iqy[particleg.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("boundaryf", [&](int _i, int _j, int _k) { return particlef.GetBoundary(_i, _j); });
    file.AddPointScaler("boundaryg", [&](int _i, int _j, int _k) { return particleg.GetBoundary(_i, _j); });
    file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[particlef.GetIndex(_i, _j)]; });
    file.AddPointScaler("beta", [&](int _i, int _j, int _k) { return beta[particleg.GetIndex(_i, _j)]; });
    
    delete[] alpha; delete[] beta;
    delete[] rhot;  delete[] uxt;   delete[] uyt;   delete[] uxtm1; delete[] uytm1; delete[] temt;  delete[] qxt;   delete[] qyt;   delete[] qxtm1; delete[] qytm1;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] imx;   delete[] imy;   delete[] item;  delete[] iqx;   delete[] iqy;
    delete[] sensitivity;

    return 0;
}