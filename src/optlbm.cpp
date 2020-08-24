#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "lbm.h"
#include "adjoint.h"
#include "src/Optimize/Solver/MMA.h"

using namespace PANSLBM2;
using namespace PANSFEM2;

int main() {
    //********************Setting parameters********************
    int tmax = 1000;
    int nx = 200;
    int ny = 100;
    double nu = 0.01;
    double r = 10.0;
    double q = 0.1;
    double u0 = 0.1;
    double alpha0 = 1.0;
    double weightlimit = 0.95;

    std::vector<double> s = std::vector<double>(nx*ny, 1.0);
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            if (pow(i - nx/3, 2.0) + pow(j - ny/2, 2.0) < pow(r, 2.0)) {
                s[ny*i + j] = 0.0;
            }
        }
    }
    MMA<double> optimizer = MMA<double>(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.01), std::vector<double>(s.size(), 1.0));
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

        //********************Get Drag********************
        LBM<double> dsolver = LBM<double>(nx, ny, nu);
        dsolver.SetPermeation([=](int _i, int _j) {
            return alpha0*q*(1.0 - s[ny*_i + _j])/(s[ny*_i + _j] + q);
        });
        dsolver.SetBoundary([=](int _i, int _j) {
            if (_i == 0) {
                return INLET;
            } else if (_i == nx - 1) {
                return OUTLET;
            } else if (_j == 0 || _j == ny - 1) {
                return MIRROR;
            } else {
                return PERIODIC;
            }
        });
        AdjointLBM<double> isolver = AdjointLBM<double>(dsolver, tmax);

        //--------------------Direct analyze--------------------
        for (isolver.t = 0; isolver.t < tmax; isolver.t++) {
            dsolver.UpdateMacro();          //  Update macroscopic values
            isolver.CaptureMacro(dsolver);  //  Capture macroscopic values
            dsolver.Collision();            //  Collision
            dsolver.Stream();               //  Stream
            dsolver.Inlet(u0, 0.0);         //  Boundary condition (inlet)
            dsolver.ExternalForce();        //  External force by Brinkman model   
        }

        //--------------------Invert analyze--------------------
        for (isolver.t = tmax - 1; isolver.t >= 0; isolver.t--) {
            isolver.UpdateMacro();          //  Update macroscopic values
            isolver.Collision();            //  Collision
            isolver.Stream();               //  Stream
            isolver.Inlet(u0, 0.0);         //  Boundary condition (inlet)
            isolver.ExternalForce();        //  External force by Brinkman model   
        }

        //--------------------Get sensitivity--------------------
        double f = 0.0;
        std::vector<double> dfds = std::vector<double>(s.size(), 0.0);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                f += dsolver.GetU(i, j)*alpha0*q*(1.0 - s[ny*i + j])/(s[ny*i + j] + q);
                dfds[ny*i + j] = isolver.GetSensitivity(i, j)*(-alpha0*q*(q + 1.0)/pow(q + s[ny*i + j], 2.0));
            }
        }

        //--------------------Export result--------------------
        std::ofstream fout("result/result" + std::to_string(k) + ".vtk");
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
        fout << "VECTORS\tu\tfloat" << std::endl;
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                fout << dsolver.GetU(i, j) << "\t" << dsolver.GetV(i, j) << "\t" << 0.0 << std::endl;
            }
        }

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