#include <iostream>
#include <fstream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0;
    
    D2Q9<double> particlef = D2Q9<double>(nx, ny);
    D2Q9<double> particleg = D2Q9<double>(nx, ny);
    for (int j = 0; j < ny - 1; j++) {
        particlef.SetBoundary(0, j, PERIODIC);
        particlef.SetBoundary(nx - 1, j, PERIODIC);
        particleg.SetBoundary(0, j, PERIODIC);
        particleg.SetBoundary(nx - 1, j, PERIODIC);
    }
    for (int i = 0; i < nx - 1; i++) {
        particlef.SetBoundary(i, 0, MIRROR);
        particlef.SetBoundary(i, ny - 1, MIRROR);
        particleg.SetBoundary(i, 0, OTHER);
        particleg.SetBoundary(i, ny - 1, OTHER);
    }

    AD<double, D2Q9, D2Q9> dsolver = AD<double, D2Q9, D2Q9>(&particlef, &particleg, nu, alpha);

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        dsolver.UpdateMacro();      //  Update macroscopic values
        dsolver.Collision();        //  Collision
        particlef.Stream();         //  Stream
        particleg.Stream();
        for (int i = 0; i < nx - 1; i++) {
            particleg.SetTemperature(i, 0, Th);
            particleg.SetTemperature(i, ny - 1, Tl);
        }                           //  Boundary condition (Fix temperature)
        dsolver.ExternalForce();    //  External force by thermal
        
        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
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
                    fout << dsolver.GetRho(i, j) << std::endl;
                }
            }

            fout << "VECTORS\tu\tfloat" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetU(0, i, j) << "\t" << dsolver.GetU(1, i, j) << "\t" << 0.0 << std::endl;
                }
            }

            fout << "SCALARS\tT\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetTemperature(i, j) << std::endl;
                }
            }
        } 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}