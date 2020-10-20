#include <iostream>
#include <fstream>
#include <chrono>

#include "../src/particle.h"
#include "../src/navierstokes.h"
#include "../src/advection.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    
    D2Q9<double> particlef = D2Q9<double>(nx, ny);
    D2Q9<double> particleg = D2Q9<double>(nx, ny);
    for (int j = 0; j < ny - 1; j++) {
        partialf.SetBoundary(0, j, PERIODIC);
        partialf.SetBoundary(nx - 1, j, PERIODIC);
        partialg.SetBoundary(0, j, PERIODIC);
        partialg.SetBoundary(nx - 1, j, PERIODIC);
    }
    for (int i = 0; i < nx - 1; i++) {
        partialf.SetBoundary(i, 0, MIRROR);
        partialf.SetBoundary(i, ny - 1, MIRROR);
        partialg.SetBoundary(i, 0, OTHER);
        partialg.SetBoundary(i, ny - 1, OTHER);
    }

    NS<double> solverf = NS<double>(&particlef, nu);
    AD<double> solverg = AD<double>(&particleg, alpha);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        std::cout << t << std::endl;
        solverf.UpdateMacro();      //  Update macroscopic values
        solverg.UpdateMacro();
        solverf.Collision();        //  Collision
        solverg.Collision();
        solverf.Stream();           //  Stream
        solverg.Stream();
        for (int i = 0; i < nx - 1; i++) {
            solverg.SetTemperature(i, 0, Tl);
            solverg.SetTemperature(i, ny - 1, Th);
        }                           //  Boundary condition (Fix temperature)
        solverf.ExternalForce();    //  External force by thermal
        
        if (t%1000 == 0) {
            std::ofstream fout("result/thermal" + std::to_string(t/1000) + ".vtk");
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
                    fout << solverf.GetRho(i, j) << std::endl;
                }
            }

            fout << "VECTORS\tu\tfloat" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << solverf.GetU(0, i, j) << "\t" << solverf.GetU(1, i, j) << "\t" << 0.0 << std::endl;
                }
            }

            fout << "SCALARS\tT\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << solverg.GetTemperature(i, j) << std::endl;
                }
            }
        } 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}