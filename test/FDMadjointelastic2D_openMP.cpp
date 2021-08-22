#include <iostream>
#include <cmath>
#include <omp.h>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/equation/adjointelastic.h"
#include "../src/utility/vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 10000, nx = 80, ny = 60, tmax = nt - 1;
    double elasticy = 0.1, rho0 = 10000.0, stress0 = 1.0;
    double sensitivity[nx*ny] = { 0.0 }, dgamma = 1.0e-5;                               //  Desugn variable and Sensitivity
    
    //--------------------Function to get objective--------------------
    auto getObjective = [=](double *_gamma, double *_rho, double *_ux, double *_uy, double *_sxx, double *_sxy, double *_syx, double *_syy) {
        //--------------------Set particles--------------------
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
        }                                                                               //  Set boundary condition
        double rx[nx*ny] = { 0.0 }, ry[nx*ny] = { 0.0 };

        //--------------------Direct analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            EL::InitialCondition(i, particle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        }                                                                               //  Set initial condition
        EL::UpdateMacro(particle, _rho, _ux, _uy, _sxx, _sxy, _syx, _syy);              //  Update macroscopic values
        for (int t = 0; t < nt - 1; t++) {
            EL::ExpandMacro(particle, _sxx, _sxy, _syx, _syy, _gamma);
            EL::Collision(elasticy, particle, _rho, _ux, _uy, _sxx, _sxy, _syx, _syy);  //  Collision
            particle.Stream();                                                          //  Stream
            for (int i = 0; i < nx; i++) {
                particle.SetStress(i, 0, 0.0, 0.0);
                particle.SetStress(i, ny - 1, 0.0, 0.0);
            }
            for (int j = 0; j < ny; j++) {
                particle.SetRS(0, j);
                if (fabs(j - 0.5*ny) < 5) {
                    particle.SetStress(nx - 1, j, 0.0, stress0);
                } else {
                    particle.SetStress(nx - 1, j, 0.0, 0.0);
                }
            }                                                                           //  Boundary condition
            EL::UpdateMacro(particle, _rho, _ux, _uy, _sxx, _sxy, _syx, _syy);          //  Update macroscopic values
            for (int i = 0; i < nx*ny; i++) {
                rx[i] += _ux[i];
                ry[i] += _uy[i];
            }
        }

        //--------------------Get objective--------------------
        double objective = 0.0;
        for (int j = 0; j < ny; j++) {
            if (fabs(j - 0.5*ny) < 5) {
                int ij = particle.GetIndex(nx - 1, j);
                objective += stress0*ry[ij];
            }
        }
        return objective; 
    };

    //--------------------FDM loop--------------------
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();
    #pragma omp parallel
    {
        double rho[nx*ny], ux[nx*ny], uy[nx*ny], sxx[nx*ny], sxy[nx*ny], syx[nx*ny], syy[nx*ny];   //  State variable
        double gamma[nx*ny];                                                              //  Desugn variable and Sensitivity
        for (int i = 0; i < nx*ny; i++) {
            rho[i] = rho0;
        }

        int id = omp_get_thread_num();
        int con = omp_get_num_threads();
        
        int begin = (nx*ny)/con*id + std::min((nx*ny)%con, id);
        int end = (nx*ny)/con*(id + 1) + std::min((nx*ny)%con, id + 1);
        for (int k = begin; k < end; k++) {
            std::cout << k << "\t";

            for (int i = 0; i < nx*ny; i++) {
                gamma[i] = 0.5;
            }
            
            gamma[k] += dgamma;
            double objectivep = getObjective(gamma, rho, ux, uy, sxx, sxy, syx, syy);

            gamma[k] -= 2.0*dgamma;
            double objectivem = getObjective(gamma, rho, ux, uy, sxx, sxy, syx, syy);

            sensitivity[k] = (objectivep - objectivem)/(2.0*dgamma);
        }
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
    VTKExport file("result/FDMadjointelastic.vtk", nx, ny);
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return sensitivity[ny*_i + _j];  });
    
    return 0;
}