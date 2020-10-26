#include <iostream>
#include <cmath>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/nsadjoint.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //********************Setting parameters********************
    int nt = 20000, nx = 100, ny = 100;
    double nu = 0.1, u0 = 0.1, rho0 = 1.0;
    double q = 0.1, alpha0 = 1.0;

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        if (0.3*ny < j && j < 0.7*ny) {
            particle.SetBoundary(0, j, OTHER);
            particle.SetBoundary(nx - 1, j, OTHER);
        } else {
            particle.SetBoundary(0, j, BARRIER);
            particle.SetBoundary(nx - 1, j, BARRIER);
        }
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, BARRIER);
        particle.SetBoundary(i, ny - 1, BARRIER);
    }

    NSAdjoint<double, D2Q9> dsolver(&particle, nu, nt);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double gamma = pow(i - 50, 2.0) + pow(j - 50, 2.0) < pow(15.0, 2.0) ? 0.1 : 0.9;
            dsolver.SetAlpha(particle.GetIndex(i, j), alpha0*q*(1.0 - gamma)/(gamma + q));
        }
    }

    //--------------------Direct analyze--------------------
    for (dsolver.t = 0; dsolver.t < nt; dsolver.t++) {
        dsolver.UpdateMacro();          //  Update macroscopic values
        dsolver.Collision();            //  Collision
        particle.Stream();              //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                double uj = -u0*(j - 0.3*ny)*(j - 0.7*ny)/(0.04*ny*ny);
                particle.SetU(0, j, uj, 0.0);
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
    dsolver.SetFt(0.0, 0.0, 0.0);
    for (dsolver.t = dsolver.t >= nt ? nt - 1 : dsolver.t; dsolver.t >= 0; dsolver.t--) {
        dsolver.iUpdateMacro();         //  Update macroscopic values
        dsolver.iCollision();           //  Collision
        particle.iStream();             //  Stream
        for (int j = 0; j < ny; j++) {
            if (0.3*ny < j && j < 0.7*ny) {
                double uj = -u0*(j - 0.3*ny)*(j - 0.7*ny)/(0.04*ny*ny);
                particle.SetiU(0, j, uj, 0.0);
                particle.SetiRho(nx - 1, j);
            }
        }                               //  Boundary condition (inlet)
        dsolver.iExternalForce();       //  External force by Brinkman model
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjoint.vtk", nx, ny);
    file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return dsolver.GetRho(particle.GetIndex(_i, _j), dsolver.tmax); });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return dsolver.GetU(0, particle.GetIndex(_i, _j), dsolver.tmax); },
        [&](int _i, int _j, int _k) { return dsolver.GetU(1, particle.GetIndex(_i, _j), dsolver.tmax); },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("dfds", [&](int _i, int _j, int _k) { return dsolver.GetSensitivity(particle.GetIndex(_i, _j)); });
    file.AddPointScaler("q", [&](int _i, int _j, int _k) { return dsolver.GetQ(particle.GetIndex(_i, _j)); });
    file.AddPointVector("v", 
        [&](int _i, int _j, int _k) { return dsolver.GetV(0, particle.GetIndex(_i, _j)); },
        [&](int _i, int _j, int _k) { return dsolver.GetV(1, particle.GetIndex(_i, _j)); },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointVector("r", 
        [&](int _i, int _j, int _k) { return dsolver.GetR(0, particle.GetIndex(_i, _j)); },
        [&](int _i, int _j, int _k) { return dsolver.GetR(1, particle.GetIndex(_i, _j)); },
        [](int _i, int _j, int _k) { return 0.0; }
    );

    return 0;
}