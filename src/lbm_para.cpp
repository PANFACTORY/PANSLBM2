//*****************************************************************************
//  Title       :   src/lbm.cpp
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/07
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>


#include "lbm.h"


using namespace PANSLBM2;


int main() {
    //--------------------Parameters--------------------
    int tmax = 10000;
    double ux0 = 0.1;
    LBM<double> solver = LBM<double>(200, 80, 0.02);
    solver.SetBarrier([=](int _i, int _j) {
        return _i == solver.ny/2 && abs(_j - solver.ny/2) <= 8;
    });


    //--------------------Loop for time step--------------------
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    for (int t = 0; t < tmax; t++) {
        solver.Inlet(0.1, 0.0);         //  Boundary condition (inlet)     
        solver.Stream();                //  Stream
        solver.UpdateMacro();           //  Update macroscopic values

        //..........Export result..........
        if (t%100 == 0) {
            std::cout << t << std::endl;

            std::ofstream fout("result/result" + std::to_string(t/100) + ".vtk");
            fout << "# vtk DataFile Version 3.0" << std::endl;
            fout << "2D flow" << std::endl;
            fout << "ASCII" << std::endl;
            fout << "DATASET\tSTRUCTURED_GRID" << std::endl;
            fout << "DIMENSIONS\t" << solver.nx << "\t" << solver.ny << "\t" << 1 << std::endl;
            
            fout << "POINTS\t" << solver.nx*solver.ny << "\t" << "float" << std::endl;
            for (int j = 0; j < solver.ny; j++) {
                for (int i = 0; i < solver.nx; i++) {
                    fout << i << "\t" << j << "\t" << 0.0 << std::endl;
                }
            }

            fout << "POINT_DATA\t" << solver.nx*solver.ny << std::endl;
            fout << "SCALARS\trho\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < solver.ny; j++) {
                for (int i = 0; i < solver.nx; i++) {
                    fout << solver.rho[i][j] << std::endl;
                }
            }

            fout << "SCALARS\tcurl\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < solver.ny; j++) {
                for (int i = 0; i < solver.nx; i++) {
                    if (i == 0 || j == 0 || i == solver.nx - 1 || j == solver.ny - 1) {
                        fout << 0.0 << std::endl;
                    } else {
                        fout << solver.v[i + 1][j] - solver.v[i - 1][j] - solver.u[i][j + 1] + solver.u[i][j - 1] << std::endl;
                    }
                }
            }

            fout << "VECTORS\tvelocity\tfloat" << std::endl;
            for (int j = 0; j < solver.ny; j++) {
                for (int i = 0; i < solver.nx; i++) {
                    fout << solver.u[i][j] << "\t" << solver.v[i][j] << "\t" << 0.0 << std::endl;
                }
            }
        }

        solver.Collision();             //  Collision 
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}