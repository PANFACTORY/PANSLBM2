#include <iostream>
#include <fstream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 100;
    double nu = 0.05, u0 = 0.1;
    
    D2Q9<double> particle = D2Q9<double>(nx, ny);
    for (int j = 0; j < ny - 1; j++) {
        particle.SetBoundary(0, j, BARRIER);
        particle.SetBoundary(nx - 1, j, BARRIER);
    }
    for (int i = 0; i < nx - 1; i++) {
        particle.SetBoundary(i, ny - 1, BARRIER);
        particle.SetBoundary(i, 0, OTHER);
    }

    NS<double, D2Q9> dsolver = NS<double, D2Q9>(&particle, nu);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        dsolver.UpdateMacro();      //  Update macroscopic values
        dsolver.Collision();        //  Collision
        particle.Stream();          //  Stream
        for (int i = 0; i < nx - 1; i++) {
            particle.SetU(i, 0, u0, 0.0);
        }                           //  Boundary condition (inlet)
        dsolver.ExternalForce();    //  External force by thermal
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            std::ofstream fout("result/ns" + std::to_string(t/1000) + ".vtk");
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
                    fout << dsolver.GetRho(i, j) << std::endl;
                }
            }

            fout << "VECTORS\tu\tfloat" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetU(0, i, j) << "\t" << dsolver.GetU(1, i, j) << "\t" << 0.0 << std::endl;
                }
            }
        } 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}