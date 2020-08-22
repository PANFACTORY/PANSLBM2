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
#include "adjoint.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 1000, nx = 200, ny = 80;
    LBM<double> dsolver = LBM<double>(nx, ny, 0.02);
    dsolver.SetPermeation([=](int _i, int _j) {
        double ganma = 1.0;
        if (abs(_i - ny/2) <= 0 && abs(_j - ny/2) <= 8) {
            ganma = 0.0;
        }
        return 0.1*(1.0 - ganma)/(ganma + 0.1);
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

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    for (isolver.t = 0; isolver.t < tmax; isolver.t++) {
        std::cout << isolver.t << std::endl;
        dsolver.UpdateMacro();          //  Update macroscopic values
        isolver.CaptureMacro(dsolver);  //  Capture macroscopic values
        dsolver.Collision();            //  Collision
        dsolver.Stream();               //  Stream
        dsolver.Inlet(0.1, 0.0);        //  Boundary condition (inlet)
        dsolver.ExternalForce();        //  External force by Brinkman model   
    }

    //--------------------Invert analyze--------------------
    for (isolver.t = tmax - 1; isolver.t >= 0; isolver.t--) {
        std::cout << isolver.t << std::endl;
        isolver.Collision();            //  Collision
        isolver.Stream();               //  Stream
        isolver.Inlet(0.1, 0.0);        //  Boundary condition (inlet)
        isolver.ExternalForce();        //  External force by Brinkman model   
    }

    //--------------------Export result--------------------
    std::ofstream fout("result/result.vtk");
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

    fout << "SCALARS\tsensitivity\tfloat" << std::endl;
    fout << "LOOKUP_TABLE\tdefault" << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            fout << isolver.GetSensitivity(i, j) << std::endl;
        }
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}