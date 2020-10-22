#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/nsadjoint.h"

using namespace PANSLBM2;

int main() {
    //********************Setting parameters********************
    int nt = 20000, nx = 100, ny = 100;
    double nu = 0.1, u0 = 0.25, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0;

    D2Q9<double> particle = D2Q9<double>(nx, ny);
    for (int j = 0; j < ny - 1; j++) {
        if (0.3*ny < j && j < 0.7*ny) {
            particle.SetBoundary(0, j, OTHER);
            particle.SetBoundary(nx - 1, j, OTHER);
        } else {
            particle.SetBoundary(0, j, BARRIER);
            particle.SetBoundary(nx - 1, j, BARRIER);
        }
    }
    for (int i = 0; i < nx - 1; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, BARRIER);
    }

    NSAdjoint<double, D2Q9> dsolver = NSAdjoint<double, D2Q9>(&particle, nu, nt);

    dsolver.SetAlpha([=](int _i, int _j) {
        double gamma = 0.9;
        if (pow(_i - 50, 2.0) + pow(_j - 50, 2.0) < pow(15.0, 2.0)) {
            gamma = 0.1;
        }
        return alpha0*q*(1.0 - gamma)/(gamma + q);
    });

    //--------------------Direct analyze--------------------
    for (dsolver.t = 0; dsolver.t < nt; dsolver.t++) {
        dsolver.UpdateMacro();          //  Update macroscopic values
        dsolver.Collision();            //  Collision
        particle.Stream();              //  Stream
        for (int j = 0; j < ny - 1; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                particle.SetU(0, j, u0, 0.0);
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
    dsolver.SwitchDirection();
    for (dsolver.t = dsolver.t >= nt ? nt - 1 : dsolver.t; dsolver.t >= 0; dsolver.t--) {
        dsolver.iUpdateMacro();         //  Update macroscopic values
        dsolver.iCollision();           //  Collision
        particle.Stream();              //  Stream
        for (int j = 0; j < ny - 1; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                particle.SetiU(0, j, u0, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                               //  Boundary condition (inlet)
        dsolver.iExternalForce();       //  External force by Brinkman model

        //--------------------Export result--------------------
        if (dsolver.t%100 == 0) {
            std::ofstream fout("result/adjoint" + std::to_string(dsolver.t/100) + ".vtk");
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

            fout << "SCALARS\trho\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetRho(i, j, dsolver.t) << std::endl;
                }
            }

            fout << "VECTORS\tu\tfloat" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetU(0, i, j, dsolver.t) << "\t" << dsolver.GetU(1, i, j, dsolver.t) << "\t" << 0.0 << std::endl;
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
                    fout << dsolver.GetV(0, i, j) << "\t" << dsolver.GetV(1, i, j) << "\t" << 0.0 << std::endl;
                }
            }

            fout << "SCALARS\tsensitivity\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetSensitivity(i, j) << std::endl;
                }
            }
        }
    }

    return 0;
}