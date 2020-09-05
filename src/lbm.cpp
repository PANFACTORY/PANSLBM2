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
#include "thermal.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 10000, nx = 200, ny = 80;
    
    LBM<double> dsolver = LBM<double>(nx, ny, 0.02);
    dsolver.SetBarrier([=](int _i, int _j) {
        return _i == ny/2 && abs(_j - ny/2) <= 8;
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

    ThermalLBM<double> tsolver = ThermalLBM<double>(nx, ny, 0.02);
    tsolver.SetBarrier([=](int _i, int _j) {
        return _i == ny/2 && abs(_j - ny/2) <= 8;
    });
    tsolver.SetBoundary([=](int _i, int _j) {
        if (_j == 0) {
            return INLET;
        } else if (_j == ny - 1) {
            return OUTLET;
        } else if (_i == 0 || _i == nx - 1) {
            return OUTLET;
        } else {
            return PERIODIC;
        }
    });

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (int t = 0; t < tmax; t++) {
        std::cout << t << std::endl;
        dsolver.UpdateMacro();          //  Update macroscopic values
        tsolver.CaptureMacro(dsolver);  //  Capture macroscopic values
        tsolver.UpdateMacro();          //  Update macroscopic values

        dsolver.Collision();            //  Collision
        tsolver.Collision();            //  Collision
        
        dsolver.Stream();               //  Stream
        tsolver.Stream();               //  Stream
        
        dsolver.Inlet(0.1, 0.0);        //  Boundary condition (inlet)
        tsolver.Inlet(100.0);           //  Boundary condition (inlet)
        
        if (t%100 == 0) {
            std::ofstream fout("result/result" + std::to_string(t) + ".vtk");
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
                    fout << tsolver.GetTemperature(i, j) << std::endl;
                }
            }
        } 
    }

    //--------------------Export result--------------------
    

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}