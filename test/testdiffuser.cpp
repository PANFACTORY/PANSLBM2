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
    int nt = 20000, nx = 100, ny = 100, nk = 40, tmax = nt;
    double nu = 0.1, u0 = 0.001, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0, scale0 = 1.0e0, weightlimit = 0.5, *alpha = new double[nx*ny];
    double *rho = new double[nt*nx*ny], *ux = new double[nt*nx*ny], *uy = new double[nt*nx*ny]; //  State variable
    double *qho = new double[nx*ny], *vx = new double[nx*ny], *vy = new double[nx*ny];          //  Adjoint variable

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
            particle.SetBoundary(0, j, OTHER);
            if (0.35*ny <= j && j <= 0.66*ny) {
                particle.SetBoundary(nx - 1, j, OTHER);
            } else {
                particle.SetBoundary(nx - 1, j, BARRIER);
            }
        }
        for (int i = 0; i < nx; i++) {
            particle.SetBoundary(i, 0, BARRIER);
            particle.SetBoundary(i, ny - 1, BARRIER);
        }

        //--------------------Set inverse permeation--------------------
        for (int i = 0; i < s.size(); i++) {            
            alpha[i] = alpha0*q*(1.0 - s[i])/(s[i] + q);
        }

        //--------------------Direct analyze--------------------
        for (int t = 0; t < nt; t++) {
            NS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);      //  Update macroscopic values
            NS2::Collision(nu, particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);    //  Collision
            particle.Stream();                                                          //  Stream
            for (int j = 0; j < ny; j++) {
                double uj = -u0*j*(j - ny)/2500.0;
                particle.SetU(0, j, uj, 0.0);
                if (0.35*ny <= j && j <= 0.66*ny) {
                    particle.SetRho(nx - 1, j, rho0, 0.0);
                }
            }                                                                           //  Boundary condition (inlet)
            NS2::ExternalForceBrinkman(particle, alpha);                                        //  External force by Brinkman model

            if (t > 0 && NS2::CheckConvergence(particle, 1.0e-4, &ux[nx*ny*t], &uy[nx*ny*t], &ux[nx*ny*(t - 1)], &uy[nx*ny*(t - 1)])) {
                std::cout << "\tt = " << t << "\t";
                tmax = t;
                break;
            }
        }

        //--------------------Invert analyze--------------------
        for (int i = 0; i < nx*ny; i++) {
            NS2::InitialCondition(i, particle, 0.0, 0.0, 0.0);
        }                                                                                           //  Set initial condition
        for (int t = tmax; t >= 0; t--) {
            ANS2::UpdateMacro(particle, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);    //  Update macroscopic values
            ANS2::Collision(nu, particle, &ux[nx*ny*t], &uy[nx*ny*t], qho, vx, vy);                 //  Collision
            particle.iStream();                                                                     //  Stream
            for (int j = 0; j < ny; j++) {
                double uj = -u0*j*(j - ny)/2500.0;
                particle.SetiU(0, j, uj, 0.0);
                if (0.35*ny <= j && j <= 0.66*ny) {
                    particle.SetiRho(nx - 1, j);
                }
            }                                                                                       //  Boundary condition (inlet)
            ANS2::ExternalForceBrinkman(particle, alpha, &rho[nx*ny*t], &ux[nx*ny*t], &uy[nx*ny*t]);        //  External force by Brinkman model
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        for (int j = 0; j < ny; j++) {
            f += rho[nx*ny*tmax + particle.GetIndex(0, j)];
            if (0.35*ny <= j && j <= 0.66*ny) {
                f -= rho[nx*ny*tmax + particle.GetIndex(nx - 1, j)];
            }
        }
        std::vector<double> dfds(s.size(), 0.0);
        for (int i = 0; i < s.size(); i++) {
            dfds[i] = scale0*3.0*(vx[i]*ux[nx*ny*tmax + i] + vy[i]*uy[nx*ny*tmax + i])*(-alpha0*q*(q + 1.0)/pow(q + s[i], 2.0));
        }

        //--------------------Export result--------------------
        VTKExport file("result/diffuser" + std::to_string(k) + ".vtk", nx, ny);
        file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[nx*ny*tmax + particle.GetIndex(_i, _j)]; });
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

    return 0;
}