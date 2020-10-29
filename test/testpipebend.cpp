#include <iostream>
#include <vector>
#include <cmath>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "src/Optimize/Solver/MMA.h"
#include "vtkexport.h"

using namespace PANSLBM2;
using namespace PANSFEM2;

int main() {
    //********************Setting parameters********************
    int nt = 20000, nx = 100, ny = 100, nk = 100, tmax = nt;
    double nu = 0.1, u0 = 0.0025, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0, scale0 = 1.0e0, weightlimit = 0.25, *alpha = new double[nx*ny];
    double *rho = new double[nt*nx*ny], *ux = new double[nt*nx*ny], *uy = new double[nt*nx*ny]; //  State variable
    double *qho = new double[nx*ny], *vx = new double[nx*ny], *vy = new double[nx*ny];
    double *sensitivity = new double[nx*ny];

    std::vector<double> s(nx*ny, 1.0);
    MMA<double> optimizer(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.5, 0.5, 0.7, 1.2, 1.0e-6);

    for (int k = 0; k < nk; k++) {
        std::cout << "k = " << k << "\t";

        //********************Get weight********************
        double g = -1.0;
        std::vector<double> dgds(s.size(), 0.0);
        for(int i = 0; i < nx*ny; i++){
            g += s[i]/(weightlimit*s.size());
            dgds[i] = 1.0/(weightlimit*s.size()); 
        }

        //********************Get PressureDrop********************
        D2Q9<double> particle(nx, ny);
        for (int j = 0; j < ny; j++) {
            if (0.7*ny < j && j < 0.9*ny) {
                particle.SetBoundary(0, j, OTHER);
            } else {
                particle.SetBoundary(0, j, BARRIER);
            }
            particle.SetBoundary(nx - 1, j, BARRIER);
        }
        for (int i = 0; i < nx; i++) {
            if (0.7*nx < i && i < 0.9*nx) {
                particle.SetBoundary(i, 0, OTHER);
            } else {
                particle.SetBoundary(i, 0, BARRIER);
            }
            particle.SetBoundary(i, ny - 1, BARRIER);
        }

        //--------------------Set inverse permeation--------------------        
        for (int i = 0; i < s.size(); i++) {            
            alpha[i] = alpha0*q*(1.0 - s[i])/(s[i] + q);
            sensitivity[i] = 0.0;
        }

        //--------------------Direct analyze--------------------
        for (int t = 0; t < nt; t++) {
            NS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);      //  Update macroscopic values
            NS2::Collision(nu, particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);    //  Collision
            particle.Stream();              //  Stream
            for (int i = 0; i < nx; i++) {
                if (0.7*nx < i && i < 0.9*nx) {
                    particle.SetRho(i, 0, rho0, 0.0);
                }
            }
            for (int j = 0; j < ny; j++) {
                if (0.7*ny < j && j < 0.9*ny) {
                    double uj = -u0*(j - 0.7*ny)*(j - 0.9*ny)/(0.01*ny*ny);
                    particle.SetU(0, j, uj, 0.0);
                }
            }                               //  Boundary condition (inlet)
            NS2::ExternalForceBrinkman(particle, alpha);        //  External force by Brinkman model

            if (t > 0 && NS2::CheckConvergence(particle, 1.0e-4, &ux[nx*ny*t], &uy[nx*ny*t], &ux[nx*ny*(t - 1)], &uy[nx*ny*(t - 1)])) {
                std::cout << "\tt = " << t << "\t";
                tmax = t;
                break;
            }
        }

        //--------------------Invert analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS2::InitialCondition(i, particle, 0.0, 0.0, 0.0);
        }
        for (int t = tmax; t >= 0; t--) {
            ANS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);         //  Update macroscopic values
            ANS2::Collision(nu, particle, &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);           //  Collision
            particle.iStream();              //  Stream
            for (int i = 0; i < nx; i++) {
                if (0.7*nx < i && i < 0.9*nx) {
                    particle.SetiRho(i, 0);
                }
            }
            for (int j = 0; j < ny; j++) {
                if (0.7*ny < j && j < 0.9*ny) {
                    double uj = -u0*(j - 0.7*ny)*(j - 0.9*ny)/(0.01*ny*ny);
                    particle.SetiU(0, j, uj, 0.0);
                }
            }                               //  Boundary condition (inlet)
            ANS2::ExternalForceBrinkman(particle, alpha, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);        //  External force by Brinkman model
            ANS2::UpdateSensitivity(particle, &ux[nx*ny*t], &uy[nx*ny*t], sensitivity);             //  Update sensitivity
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        for (int j = 0; j < ny; j++) {
            if (0.7*ny < j && j < 0.9*ny) {
                f += rho[nx*ny*tmax + particle.GetIndex(0, j)];
            }
        }
        for (int i = 0; i < nx; i++) {
            if (0.7*nx < i && i < 0.9*nx) {
                f -= rho[nx*ny*tmax + particle.GetIndex(i, 0)];
            }
        }
        std::vector<double> dfds(s.size(), 0.0);
        for (int i = 0; i < s.size(); i++) {  
            dfds[i] = scale0*sensitivity[i]*(-alpha0*q*(q + 1.0)/pow(q + s[i], 2.0));
        }

        //--------------------Export result--------------------
        VTKExport file("result/pipebend" + std::to_string(k) + ".vtk", nx, ny);
        file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[nx*ny*tmax + particle.GetIndex(_i, _j)]/3.0; });
        file.AddPointVector("u", 
            [&](int _i, int _j, int _k) { return ux[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
            [&](int _i, int _j, int _k) { return uy[nx*ny*tmax + particle.GetIndex(_i, _j)]; },
            [](int _i, int _j, int _k) { return 0.0; }
        );
        file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[particle.GetIndex(_i, _j)]; });
        file.AddPointScaler("dfds", [&](int _i, int _j, int _k) { return dfds[particle.GetIndex(_i, _j)]; });        
        file.AddPointScaler("q", [&](int _i, int _j, int _k) { return qho[particle.GetIndex(_i, _j)]; });
        file.AddPointVector("v", 
            [&](int _i, int _j, int _k) { return vx[particle.GetIndex(_i, _j)]; },
            [&](int _i, int _j, int _k) { return vy[particle.GetIndex(_i, _j)]; },
            [](int _i, int _j, int _k) { return 0.0; }
        );

        std::cout << "Objective:\t" << f << "\tWeight:\t" << g << std::endl;

        //********************Check convergence********************
        if(optimizer.IsConvergence(f)){
            std::cout << std::endl << "-----Optimized-----" << std::endl;
            break;
        }

        //********************Update variable********************
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
    }

    delete[] alpha;
    delete[] rho;   delete[] ux;    delete[] uy;
    delete[] qho;   delete[] vx;    delete[] vy;
    delete[] sensitivity;

    return 0;
}