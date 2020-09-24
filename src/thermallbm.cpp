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

#include "nsadd2q9.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50;
    double nu = 0.02;
    
    NSADd2q9<double> dsolver = NSADd2q9<double>(nx, ny, nu);
    dsolver.SetBoundary([=](int _i, int _j) {
        if (_j == 0 || _j == ny - 1) {
            return MIRROR;
        } else {
            return PERIODIC;
        }
    }, [=](int _i, int _j) {
        if (_j == 0 || _j == ny - 1) {
            return INLET;
        } else {
            return PERIODIC;
        }
    });

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        std::cout << t << std::endl;
        dsolver.UpdateMacro();      //  Update macroscopic values
        dsolver.Collision();        //  Collision
        dsolver.Stream();           //  Stream
        dsolver.Inlet(2.0, 1.0);    //  Boundary condition (inlet)
        dsolver.ExternalForce();    //  External force by thermal
        
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
            fout << "VECTORS\tu\tfloat" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dsolver.GetU(i, j) << "\t" << dsolver.GetV(i, j) << "\t" << 0.0 << std::endl;
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