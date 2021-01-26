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
    double q = 0.1, alpha0 = 50.0/(double)nx, beta0 = 0.1/(double)nx, *alpha = new double[nx*ny], *beta = new double[nx*ny];   //  Inverse permeation
    double **rho = new double*[nt], **ux = new double*[nt], **uy = new double*[nt];
    double **tem = new double*[nt], **qx = new double*[nt], **qy = new double*[nt];                     //  State variable
    for (int t = 0; t < nt; t++) {
        rho[t] = new double[nx*ny]; ux[t] = new double[nx*ny];  uy[t] = new double[nx*ny];
        tem[t] = new double[nx*ny]; qx[t] = new double[nx*ny];  qy[t] = new double[nx*ny];
    }
    double *irho = new double[nx*ny], *iux = new double[nx*ny], *iuy = new double[nx*ny], *imx = new double[nx*ny], *imy = new double[nx*ny];
    double *item = new double[nx*ny], *iqx = new double[nx*ny], *iqy = new double[nx*ny];               //  Adjoint variable
    double *sensitivity = new double[nx*ny];                                                            //  Sensitivity

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
    }                                                                           //  Set boundary condition

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
    }                                                                                           //  Set initial condition
    NS::UpdateMacro(particlef, rho[0], ux[0], uy[0]);
    AD::UpdateMacro(particleg, tem[0], qx[0], qy[0], ux[0], uy[0]);                             //  Update macroscopic values
    for (int t = 0; t < nt - 1; t++) {
        //  Navier-Stokes
        NS::Collision(viscosity, particlef, rho[t], ux[t], uy[t]);                              //  Collision
        particlef.Stream();                                                                     //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particlef.SetU(0, j, uj, 0.0);
                particlef.SetRho(nx - 1, j, rho0, 0.0);   
            }
        }                                                                                       //  Boundary condition (inlet)
        NS::UpdateMacro(particlef, rho[t + 1], ux[t + 1], uy[t + 1]);                           //  Update macroscopic values
        NS::ExternalForceBrinkman(particlef, rho[t + 1], ux[t + 1], uy[t + 1], alpha, alpha);   //  External force term
        
        //  Advection
        AD::Collision(diffusivity, particleg, tem[t], ux[t], uy[t]);                            //  Collision
        particleg.Stream();                                                                     //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                particleg.SetTemperature(0, j, tem0);
                int ij = particlef.GetIndex(nx - 1, j);
                particleg.SetFlux(nx - 1, j, ux[t][ij], uy[t][ij], q0);
            } else {
                particleg.SetFlux(0, j, 0.0, 0.0, q0);
                particleg.SetFlux(nx - 1, j, 0.0, 0.0, q0);
            }
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetFlux(i, ny - 1, 0.0, 0.0, q0);
        }                                                                                       //  Boundary condition (inlet)
        AD::UpdateMacro(particleg, tem[t + 1], qx[t + 1], qy[t + 1], ux[t + 1], uy[t + 1]);     //  Update macroscopic values
        AD::ExternalForceHeatgeneration(particleg, tem[t + 1], beta);                           //  External force term
        
        //  Check convergence
        if (t > 0 && NS::CheckConvergence(particlef, epsdu, ux[t + 1], uy[t + 1], ux[t], uy[t]) && NS::CheckConvergence(particleg, epsdq, qx[t + 1], qy[t + 1], qx[t], qy[t])) {
            std::cout << "\tt = " << t << "\t";
            tmax = t + 1;
            break;
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 0.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 0.0, 0.0, 0.0);
    }                                                                                           //  Set initial condition
    AAD::UpdateMacro(particleg, item, iqx, iqy);
    ANS::UpdateMacro(particlef, rho[tmax], ux[tmax], uy[tmax], irho, iux, iuy, imx, imy);       //  Update macroscopic values
    for (int t = tmax; t >= 0; t--) {
        //  Adjoint advection
        AAD::Collision(diffusivity, particleg, ux[tmax], uy[tmax], item, iqx, iqy);                   //  Collision
        particleg.iStream();                                                                    //  Stream
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                particleg.SetiTemperature(0, j);
                int ij = particleg.GetIndex(nx - 1, j);
                particleg.SetiFlux(nx - 1, j, ux[tmax][ij], uy[tmax][ij]);
            } else {
                particleg.SetiFlux(0, j, 0.0, 0.0);
                particleg.SetiFlux(nx - 1, j, 0.0, 0.0);
            }
        }
        for (int i = 0; i < nx; i++) {
            particleg.SetiFlux(i, ny - 1, 0.0, 0.0);
        }                                                                                       //  Boundary condition (inlet)
        AAD::UpdateMacro(particleg, item, iqx, iqy);                                            //  Update macroscopic values
        AAD::ExternalForceHeatexchange(particleg, item, beta);                                  //  External force by Brinkman model
        
        //  Adjoint Navier-Stokes
        ANS::Collision(viscosity, particlef, ux[tmax], uy[tmax], irho, iux, iuy);               //  Collision
        particlef.iStream();                                                                    //  Stream
        ANS::ExternalForceHeatexchange(diffusivity, particlef, rho[tmax], ux[tmax], uy[tmax], tem[tmax], iqx, iqy);
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                particlef.SetiU(0, j, uj, 0.0);
                int ij = particlef.GetIndex(nx - 1, j);
                particlef.SetiRhoFlux(particleg, nx - 1, j, rho0, ux[tmax][ij], uy[tmax][ij], tem[tmax][ij]);
                //particlef.SetiUFlux(particleg, nx - 1, j, rho0, ux[t][ij], uy[t][ij], tem[t][ij]);
                //particlef.SetiUFlux(particleg, nx - 1, j, rho0, uj, 0.0, tem[t][ij]);
                //particlef.SetiU(nx - 1, j, ux[t][ij], uy[t][ij]);
                //particlef.SetiRho(nx - 1, j);
            }
        }                                                                                       //  Boundary condition (inlet)
        ANS::UpdateMacro(particlef, rho[tmax], ux[tmax], uy[tmax], irho, iux, iuy, imx, imy);   //  Update macroscopic values 
        ANS::ExternalForceBrinkman(particlef, rho[tmax], iux, iuy, imx, imy, alpha, alpha);     //  External force by Brinkman model

        if (t%300 == 0) {
            VTKExport file("result/adjointadvection_AL50_BE1e-1_q1e-1_outrhoflux_static_log" + std::to_string(t/300) + ".vtk", nx, ny);
            file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[t][particlef.GetIndex(_i, _j)]/3.0; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[t][particlef.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[t][particlef.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[t][particleg.GetIndex(_i, _j)]; });
            file.AddPointVector("q", 
                [&](int _i, int _j, int _k) { return qx[t][particleg.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[t][particleg.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            //file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[particlef.GetIndex(_i, _j)];  });
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
        }
    }

    //--------------------Analyse sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
            int ij = particlef.GetIndex(i, j);
            sensitivity[ij] = 3.0*(imx[ij]*ux[tmax][ij] + imy[ij]*uy[tmax][ij])*(-alpha0*q*(q + 1.0)/pow(q + gamma, 2.0)) - (1.0 - tem[tmax][ij])*(1.0 + item[ij])*(-beta0*q*(q + 1.0)/pow(q + gamma, 2.0));
            if (sensitivitymax < fabs(sensitivity[ij])) {
                sensitivitymax = fabs(sensitivity[ij]);
            }
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointadvection_AL50_BE1e-1_q1e-1_outrhoflux_static.vtk", nx, ny);
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
    for (int t = 0; t < nt; t++) {
        delete[] rho[t];    delete[] ux[t]; delete[] uy[t]; delete[] tem[t];    delete[] qx[t]; delete[] qy[t];
    }
    delete[] rho;   delete[] ux;    delete[] uy;    delete[] tem;   delete[] qx;    delete[] qy;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] imx;   delete[] imy;   delete[] item;  delete[] iqx;   delete[] iqy;
    delete[] sensitivity;

    return 0;
}