#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/nsadjoint.h"
#include "src/Optimize/Solver/MMA.h"

using namespace PANSLBM2;
using namespace PANSFEM2;

int main() {
    //********************Setting parameters********************
    int nt = 30000, nx = 100, ny = 100;
    double nu = 0.1, u0 = 0.025, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0, scale0 = 1.0e0, weightlimit = 0.25;

    std::vector<double> s = std::vector<double>(nx*ny, 1.0);
    MMA<double> optimizer = MMA<double>(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), std::vector<double>(s.size(), 1.0));
	optimizer.SetParameters(1.0e-5, 0.1, 0.2, 0.5, 0.7, 1.2, 1.0e-6);

    for (int k = 0; k < 50; k++) {
        std::cout << "k = " << k << "\t";

        //********************Get weight********************
        double g = -1.0;
        std::vector<double> dgds = std::vector<double>(s.size(), 0.0);
        for(int i = 0; i < nx*ny; i++){
            g += s[i]/(weightlimit*s.size());
            dgds[i] = 1.0/(weightlimit*s.size()); 
        }

        //********************Get PressureDrop********************
        D2Q9<double> particle = D2Q9<double>(nx, ny);
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

        NSAdjoint<double, D2Q9> dsolver = NSAdjoint<double, D2Q9>(&particle, nu, nt);
        
        for (int i = 0; i < s.size(); i++) {            
            dsolver.SetAlpha(i, alpha0*q*(1.0 - s[i])/(s[i] + q));
        }

        //--------------------Direct analyze--------------------
        for (dsolver.t = 0; dsolver.t < nt; dsolver.t++) {
            dsolver.UpdateMacro();          //  Update macroscopic values
            dsolver.Collision();            //  Collision
            particle.Stream();              //  Stream
            for (int i = 0; i < nx; i++) {
                if (0.7*nx < i && i < 0.9*nx) {
                    particle.SetRho(i, 0, rho0, 0.0);
                    particle.SetU(0, i, u0, 0.0);
                }
            }                               //  Boundary condition (inlet)
            dsolver.ExternalForce();        //  External force by Brinkman model

            if (dsolver.t > 0 && dsolver.CheckConvergence(1.0e-4)) {
                std::cout << "\tt = " << dsolver.t << "\t";
                break;
            }
        }

        //--------------------Invert analyze--------------------
        dsolver.SwitchDirection();
        for (dsolver.t = dsolver.t >= nt ? nt - 1 : dsolver.t; dsolver.t >= 0; dsolver.t--) {
            dsolver.iUpdateMacro();         //  Update macroscopic values
            dsolver.iCollision();           //  Collision
            particle.Stream();              //  Stream
            for (int i = 0; i < nx; i++) {
                if (0.7*nx < i && i < 0.9*nx) {
                    particle.SetiRho(i, 0);
                    particle.SetiU(0, i, u0, 0.0);
                }
            }                               //  Boundary condition (inlet)
            dsolver.iExternalForce();       //  External force by Brinkman model   
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        std::vector<double> dfds = std::vector<double>(s.size(), 0.0);
        for (int i = 0; i < s.size(); i++) {  
            //f += scale0*dsolver.GetU(0, i, j)*alpha0*q*(1.0 - s[ny*i + j])/(s[ny*i + j] + q);
            dfds[i] = scale0*dsolver.GetSensitivity(i)*(-alpha0*q*(q + 1.0)/pow(q + s[i], 2.0));
        }

        //--------------------Export result--------------------
        std::ofstream fout("result/optimizep" + std::to_string(k) + ".vtk");
        fout << "# vtk DataFile Version 3.0" << std::endl;
        fout << "2D flow" << std::endl;
        fout << "ASCII" << std::endl;
        fout << "DATASET\tSTRUCTURED_GRID" << std::endl;
        fout << "DIMENSIONS\t" << nx << "\t" << ny << "\t" << 1 << std::endl;
        
        fout << "POINTS\t" << nx*ny << "\t" << "float" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << i << "\t" << j << "\t" << 0.0 << std::endl;
            }
        }

        fout << "POINT_DATA\t" << nx*ny << std::endl;

        fout << "SCALARS\ts\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << s[particle.GetIndex(i, j)] << std::endl;
            }
        }

        fout << "SCALARS\tsensitivity\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dfds[particle.GetIndex(i, j)] << std::endl;
            }
        }

        /*fout << "SCALARS\trho\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetRho(i, j) << std::endl;
            }
        }

        fout << "VECTORS\tu\tfloat" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetU(0, i, j) << "\t" << dsolver.GetU(1, i, j) << "\t" << 0.0 << std::endl;
            }
        }*/

        fout << "SCALARS\tq\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetQ(particle.GetIndex(i, j)) << std::endl;
            }
        }

        fout << "VECTORS\tv\tfloat" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetV(0, particle.GetIndex(i, j)) << "\t" << dsolver.GetV(1, particle.GetIndex(i, j)) << "\t" << 0.0 << std::endl;
            }
        }

        std::cout << "Objective:\t" << f << "\tWeight:\t" << g << std::endl;

        //********************Check convergence********************
        if(optimizer.IsConvergence(f)){
            //std::cout << std::endl << "-----Optimized-----" << std::endl;
            //break;
        }

        //********************Update variable********************
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
    }

    return 0;
}