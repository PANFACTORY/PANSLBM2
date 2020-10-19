#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "nsadjointd2q9_2.h"
#include "src/Optimize/Solver/MMA.h"

using namespace PANSLBM2;
using namespace PANSFEM2;

int main() {
    //********************Setting parameters********************
    int nt = 20000;
    int nx = 100;
    int ny = 100;
    double nu = 0.1;
    double q = 0.1;
    double u0 = 0.025;
    double rho0 = 1.0;
    double alpha0 = 1.0;
    double scale0 = 1.0e0;
    double weightlimit = 0.25;

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
        D2Q9<double> partial = D2Q9<double>(nx, ny);
        NSAdjointd2q9<double> dsolver = NSAdjointd2q9<double>(&partial, nt, nu);
        dsolver.SetAlpha([=](int _i, int _j) {
            return alpha0*q*(1.0 - s[ny*_i + _j])/(s[ny*_i + _j] + q);
        });
        partial.SetBoundary([=](int _i, int _j) {
            if ((_i == 0 && ((int)(0.7*ny) < _j && _j < (int)(0.9*ny))) || (((int)(0.7*nx) < _i && _i < (int)(0.9*nx)) && _j == 0)) {
                return INLET;
            } else {
                return BARRIER;
            }
        });

        //--------------------Direct analyze--------------------
        for (dsolver.t = 0; dsolver.t < nt; dsolver.t++) {
            dsolver.UpdateMacro();          //  Update macroscopic values
            dsolver.Collision();            //  Collision
            partial.Stream();               //  Stream
            dsolver.Inlet(u0, 0.0, rho0);   //  Boundary condition (inlet)
            dsolver.ExternalForce();        //  External force by Brinkman model

            if (dsolver.t > 0 && dsolver.CheckConvergence(1.0e-4)) {
                std::cout << "\tt = " << dsolver.t << "\t";
                break;
            }
        }

        //--------------------Invert analyze--------------------
        dsolver.SwitchDirection();
        for (dsolver.t = dsolver.t >= nt ? nt - 1 : dsolver.t; dsolver.t >= 0; dsolver.t--) {
            dsolver.UpdateMacro();          //  Update macroscopic values
            dsolver.Collision();            //  Collision
            partial.Stream();               //  Stream
            dsolver.Inlet(u0, 0.0, rho0);   //  Boundary condition (inlet)
            dsolver.ExternalForce();        //  External force by Brinkman model   
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        std::vector<double> dfds = std::vector<double>(s.size(), 0.0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                f += scale0*dsolver.GetUx(i, j)*alpha0*q*(1.0 - s[ny*i + j])/(s[ny*i + j] + q);
                dfds[ny*i + j] = scale0*dsolver.GetSensitivity(i, j)*(-alpha0*q*(q + 1.0)/pow(q + s[ny*i + j], 2.0));
            }
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
                fout << s[ny*i + j] << std::endl;
            }
        }

        fout << "SCALARS\tsensitivity\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dfds[ny*i + j] << std::endl;
            }
        }

        fout << "SCALARS\trho\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetRho(i, j) << std::endl;
            }
        }

        fout << "VECTORS\tu\tfloat" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetUx(i, j) << "\t" << dsolver.GetUy(i, j) << "\t" << 0.0 << std::endl;
            }
        }

        fout << "SCALARS\tq\tfloat" << std::endl;
        fout << "LOOKUP_TABLE\tdefault" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetQ(i, j) << std::endl;
            }
        }

        fout << "VECTORS\tv\tfloat" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetVx(i, j) << "\t" << dsolver.GetVy(i, j) << "\t" << 0.0 << std::endl;
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