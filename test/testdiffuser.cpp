#include <iostream>
#include <vector>
#include <cmath>

#include "../src/particle/d2q9.h"
#include "../src/equation/nsadjoint.h"
#include "src/Optimize/Solver/MMA.h"
#include "vtkexport.h"

using namespace PANSLBM2;
using namespace PANSFEM2;

int main() {
    //********************Setting parameters********************
    int nt = 20000, nx = 100, ny = 100, nk = 40;
    double nu = 0.1, u0 = 0.001, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0, scale0 = 1.0e0, weightlimit = 0.5;

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

        NSAdjoint<double, D2Q9> dsolver(&particle, nu, nt);
        for (int i = 0; i < s.size(); i++) {            
            dsolver.SetAlpha(i, alpha0*q*(1.0 - s[i])/(s[i] + q));
        }

        //--------------------Direct analyze--------------------
        for (dsolver.t = 0; dsolver.t < nt; dsolver.t++) {
            dsolver.UpdateMacro();          //  Update macroscopic values
            dsolver.Collision();            //  Collision
            particle.Stream();              //  Stream
            for (int j = 0; j < ny; j++) {
                double uj = -u0*j*(j - ny)/2500.0;
                particle.SetU(0, j, uj, 0.0);
                if (0.35*ny <= j && j <= 0.66*ny) {
                    particle.SetRho(nx - 1, j, rho0, 0.0);
                }
            }                               //  Boundary condition (inlet)
            dsolver.ExternalForce();        //  External force by Brinkman model

            if (dsolver.t > 0 && dsolver.CheckConvergence(1.0e-4)) {
                std::cout << "\tt = " << dsolver.t << "\t";
                break;
            }
        }

        //--------------------Invert analyze--------------------
        for (int i = 0; i < s.size(); i++) {            
            dsolver.SetFt(i, 0.0, 0.0, 0.0);
        }
        for (dsolver.t = dsolver.t >= nt ? nt - 1 : dsolver.t; dsolver.t >= 0; dsolver.t--) {
            dsolver.iUpdateMacro();         //  Update macroscopic values
            dsolver.iCollision();           //  Collision
            particle.iStream();              //  Stream
            for (int j = 0; j < ny; j++) {
                double uj = -u0*j*(j - ny)/2500.0;
                particle.SetiU(0, j, uj, 0.0);
                if (0.35*ny <= j && j <= 0.66*ny) {
                    particle.SetiRho(nx - 1, j);
                }
            }                               //  Boundary condition (inlet)
            dsolver.iExternalForce();       //  External force by Brinkman model   
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        for (int j = 0; j < ny; j++) {
            f += dsolver.GetRho(particle.GetIndex(0, j), dsolver.tmax);
            if (0.35*ny <= j && j <= 0.66*ny) {
                f -= dsolver.GetRho(particle.GetIndex(nx - 1, j), dsolver.tmax);
            }
        }
        std::vector<double> dfds(s.size(), 0.0);
        for (int i = 0; i < s.size(); i++) {  
            dfds[i] = scale0*dsolver.GetSensitivity(i)*(-alpha0*q*(q + 1.0)/pow(q + s[i], 2.0));
        }

        //--------------------Export result--------------------
        VTKExport file("result/optimizep" + std::to_string(k) + ".vtk", nx, ny);
        file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return dsolver.GetRho(particle.GetIndex(_i, _j), dsolver.tmax); });
        file.AddPointVector("u", 
            [&](int _i, int _j, int _k) { return dsolver.GetU(0, particle.GetIndex(_i, _j), dsolver.tmax); },
            [&](int _i, int _j, int _k) { return dsolver.GetU(1, particle.GetIndex(_i, _j), dsolver.tmax); },
            [](int _i, int _j, int _k) { return 0.0; }
        );
        file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[particle.GetIndex(_i, _j)]; });
        file.AddPointScaler("dfds", [&](int _i, int _j, int _k) { return dfds[particle.GetIndex(_i, _j)]; });
        file.AddPointScaler("q", [&](int _i, int _j, int _k) { return dsolver.GetQ(particle.GetIndex(_i, _j)); });
        file.AddPointVector("v", 
            [&](int _i, int _j, int _k) { return dsolver.GetV(0, particle.GetIndex(_i, _j)); },
            [&](int _i, int _j, int _k) { return dsolver.GetV(1, particle.GetIndex(_i, _j)); },
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

    return 0;
}