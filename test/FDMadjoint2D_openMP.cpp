#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 10000, nx = 200, ny = 100, tmax = nt - 1;
    double viscosity = 0.1; 
    double u0 = 0.0109, rho0 = 1.0, q0 = 0.0, epsdu = 1.0e-8;    
    double q = 0.01, alpha0 = 50.0/(double)nx, dgamma = 1.0e-5;
    double *sensitivity = new double[nx*ny];
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 0.0;
    }

    //--------------------Function to get objective--------------------
    auto getObjective = [=](double *_alpha, double *_gamma, double *_rho, double *_ux, double *_uy) {
        //--------------------Set particles--------------------
        D2Q9<double> particlef(nx, ny);
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                particlef.SetBoundary(0, j, OTHER);
                particlef.SetBoundary(nx - 1, j, OTHER);
            } else {
                particlef.SetBoundary(0, j, BARRIER);
                particlef.SetBoundary(nx - 1, j, BARRIER);
            }
        }
        for (int i = 0; i < nx; i++) {
            particlef.SetBoundary(i, 0, MIRROR);
            particlef.SetBoundary(i, ny - 1, BARRIER);
        }                                                           //  Set boundary condition

        //--------------------Set design variable--------------------
        for (int i = 0; i < nx*ny; i++) {
            _alpha[i] = alpha0*q*(1.0 - _gamma[i])/(_gamma[i] + q);
        }

        //--------------------Direct analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
        }                                                           //  Set initial condition
        NS::UpdateMacro(particlef, _rho, _ux, _uy);    
        for (int t = 0; t < nt - 1; t++) {
            NS::Collision(viscosity, particlef, _rho, _ux, _uy);    //  Collision
            particlef.Stream();                                     //  Stream
            for (int j = 0; j < ny; j++) {
                if (j < 0.33*ny) {
                    double uj = -u0*(j - 0.33*ny)*(j + 0.33*ny)/(0.33*ny*0.33*ny);
                    particlef.SetU(0, j, uj, 0.0);
                    particlef.SetRho(nx - 1, j, rho0, 0.0);
                }
            }
            NS::UpdateMacro(particlef, _rho, _ux, _uy);             //  Update macroscopic values
            NS::ExternalForceBrinkman(particlef, _rho, _ux, _uy, _alpha, _alpha);   //  External force term
        }

        //--------------------Get objective--------------------
        double objective = 0.0;
        for (int j = 0; j < ny; j++) {
            if (j < 0.33*ny) {
                objective += (_rho[particlef.GetIndex(0, j)] - _rho[particlef.GetIndex(nx - 1, j)]);
            }
        }
        return objective/3.0; 
    };

    //--------------------FDM loop--------------------
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    #pragma omp parallel
    {
        double *alpha = new double[nx*ny], *gamma = new double[nx*ny];                      //  Inverse permeation
        double *rho = new double[nx*ny], *ux = new double[nx*ny], *uy = new double[nx*ny];  //  State variable

        int id = omp_get_thread_num();
        int con = omp_get_num_threads();
        
        int begin = (nx*ny)/con*id + std::min((nx*ny)%con, id);
        int end = (nx*ny)/con*(id + 1) + std::min((nx*ny)%con, id + 1);
        for (int k = begin; k < end; k++) {
            std::cout << k << "\t";

            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    gamma[ny*i + j] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
                }
            }
            
            gamma[k] += dgamma;
            double objectivep = getObjective(alpha, gamma, rho, ux, uy);

            gamma[k] -= 2.0*dgamma;
            double objectivem = getObjective(alpha, gamma, rho, ux, uy);

            sensitivity[k] = (objectivep - objectivem)/(2.0*dgamma);
        }

        delete[] alpha; delete[] gamma;
        delete[] rho;   delete[] ux;    delete[] uy;
    }
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::endl << "time:" << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;

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
    VTKExport file("result/FDMadjoint_AL50_q1e-2.vtk", nx, ny);
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[ny*_i + _j];  });
    
    delete[] sensitivity;

    return 0;
}