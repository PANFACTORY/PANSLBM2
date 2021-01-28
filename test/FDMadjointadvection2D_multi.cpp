#include <iostream>
#include <cmath>
#include <thread>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/equation/adjointadvection.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int con = std::thread::hardware_concurrency();
    int nt = 30000, nx = 100, ny = 50, tmax = nt - 1;
    double viscosity = 0.1, diffusivity = 0.1/6.0; 
    double u0 = 0.0218, rho0 = 1.0, q0 = 0.0, tem0 = 0.0, epsdu = 1.0e-8, epsdq = 1.0e-8;    
    double q = 0.1, alpha0 = 100.0/(double)nx, beta0 = 0.1/(double)nx, dgamma = -1.0e-5;
    double **alpha = new double*[con], **beta = new double*[con], **gamma = new double*[con];   //  Inverse permeation
    double **rho = new double*[con], **ux = new double*[con], **uy = new double*[con];
    double **tem = new double*[con], **qx = new double*[con], **qy = new double*[con];          //  State variable
    for (int i = 0; i < con; i++) {
        alpha[i] = new double[nx*ny];   beta[i] = new double[nx*ny];    gamma[i] = new double[nx*ny];
        rho[i] = new double[nx*ny];     ux[i] = new double[nx*ny];      uy[i] = new double[nx*ny];
        tem[i] = new double[nx*ny];     qx[i] = new double[nx*ny];      qy[i] = new double[nx*ny];
    }
    double *sensitivity = new double[nx*ny];                                                    //  Sensitivity
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 0.0;
    }

    //--------------------Function to get objective--------------------
    auto getObjective = [=](double *_alpha, double *_beta, double *_gamma, double *_rho, double *_ux, double *_uy, double *_tem, double *_qx, double *_qy) {
        //--------------------Set particles--------------------
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
        }                                                           //  Set boundary condition

        //--------------------Set design variable--------------------
        for (int i = 0; i < nx*ny; i++) {
            _alpha[i] = alpha0*q*(1.0 - _gamma[i])/(_gamma[i] + q);
            _beta[i] = beta0*q*(1.0 - _gamma[i])/(_gamma[i] + q);
        }

        //--------------------Direct analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
            AD::InitialCondition(i, particleg, 1.0, 0.0, 0.0);
        }                                                           //  Set initial condition
        for (int t = 0; t < nt; t++) {
            NS::UpdateMacro(particlef, _rho, _ux, _uy);
            AD::UpdateMacro(particleg, _tem, _qx, _qy, _ux, _uy);        //  Update macroscopic values
            NS::Collision(viscosity, particlef, _rho, _ux, _uy);
            AD::Collision(diffusivity, particleg, _tem, _ux, _uy);     //  Collision
            particlef.Stream();
            particleg.Stream();                                     //  Stream
            for (int j = 0; j < ny; j++) {
                if (j < 0.33*ny) {
                    double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                    particlef.SetU(0, j, uj, 0.0);
                    particlef.SetRho(nx - 1, j, rho0, 0.0);
                    particleg.SetTemperature(0, j, tem0);
                    int ij = particlef.GetIndex(nx - 1, j);
                    particleg.SetFlux(nx - 1, j, _ux[ij], _uy[ij], q0);
                } else {
                    particleg.SetFlux(0, j, 0.0, 0.0, q0);
                    particleg.SetFlux(nx - 1, j, 0.0, 0.0, q0);
                }
            }
            for (int i = 0; i < nx; i++) {
                particleg.SetFlux(i, ny - 1, 0.0, 0.0, q0);
            }                                                       //  Boundary condition (inlet)
            NS::ExternalForceBrinkman(particlef, _alpha, _alpha);
            AD::ExternalForceHeatgeneration(particleg, _beta);       //  External force term
        }

        //--------------------Get objective--------------------
        double objective = 0.0;
        for (int i = 0; i < nx*ny; i++) {
            objective += -_beta[i]*(1.0 - _tem[i]);
        }
        return objective; 
    };

    //--------------------Get objective of neutral--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            gamma[0][ny*i + j] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
        }
    }
    double objective_base = getObjective(alpha[0], beta[0], gamma[0], rho[0], ux[0], uy[0], tem[0], qx[0], qy[0]);
    std::cout << objective_base << std::endl;

    //--------------------FDM loop--------------------
    std::thread *threads = new std::thread[con];
    for (int th = 0; th < con; th++) {
        threads[th] = std::thread([=, &sensitivity](double *_alpha, double *_beta, double *_gamma, double *_rho, double *_ux, double *_uy, double *_tem, double *_qx, double *_qy){
            int lbegin = nx/con*th + std::min(nx%con, th);
            int lend = nx/con*(th + 1) + std::min(nx%con, th + 1);
            for (int l = lbegin; l < lend; l++) {
                std::cout << l << "\t";
                int k = ny*l + 25;
                
                //--------------------Set design variable--------------------
                for (int i = 0; i < nx; i++) {
                    for (int j = 0; j < ny; j++) {
                        _gamma[ny*i + j] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 1.0;
                    }
                }
                _gamma[k] += dgamma;
                double objective = getObjective(_alpha, _beta, _gamma, _rho, _ux, _uy, _tem, _qx, _qy);
                sensitivity[k] = (objective - objective_base)/dgamma;
                std::cout << objective << "\t" << sensitivity[k] << std::endl;
            }
        }, alpha[th], beta[th], gamma[th], rho[th], ux[th], uy[th], tem[th], qx[th], qy[th]);

        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(th, &cpuset);
        int rc = pthread_setaffinity_np(threads[th].native_handle(), sizeof(cpu_set_t), &cpuset);
        if (rc != 0) {
            std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";
        }
    }
    for (int th = 0; th < con; th++) {
        threads[th].join();
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
    VTKExport file("result/FDMadjointadvectionY15.vtk", nx, ny);
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[ny*_i + _j];  });
    
    delete[] threads;
    for (int i = 0; i < con; i++) {
        delete[] alpha[i];  delete[] beta[i];   delete[] gamma[i];
        delete[] rho[i];    delete[] ux[i];     delete[] uy[i];
        delete[] tem[i];    delete[] qx[i];     delete[] qy[i];
    }
    delete[] alpha; delete[] beta;  delete[] gamma;
    delete[] rho;   delete[] ux;    delete[] uy;
    delete[] tem;   delete[] qx;    delete[] qy;
    delete[] sensitivity;

    return 0;
}