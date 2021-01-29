#include <iostream>
#include <cmath>

#include "../src/particle/d2q9.h"
#include "../src/equation/elastic.h"
#include "../src/equation/adjointelastic.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Setting parameters--------------------
    int nt = 10000, nx = 80, ny = 60, tmax = nt - 1;
    double elasticy = 0.1, rho0 = 10000.0, stress0 = 1.0;
    double rho[nx*ny], ux[nx*ny], uy[nx*ny], rx[nx*ny] = { 0.0 }, ry[nx*ny] = { 0.0 }, sxx[nx*ny], sxy[nx*ny], syx[nx*ny], syy[nx*ny];   //  State variable
    double irho[nx*ny], imx[nx*ny], imy[nx*ny], isxx[nx*ny], isxy[nx*ny], isyx[nx*ny], isyy[nx*ny];                 //  Adjoint variable
    double gamma[nx*ny], sensitivity[nx*ny] = { 0.0 };                                                              //  Desugn variable and Sensitivity
    for (int i = 0; i < nx*ny; i++) {
        rho[i] = rho0;
    }

    D2Q9<double> particle(nx, ny);
    for (int j = 0; j < ny; j++) {
        particle.SetBoundary(0, j, OTHER);
        particle.SetBoundary(nx - 1, j, OTHER);
    }
    for (int i = 0; i < nx; i++) {
        particle.SetBoundary(i, 0, OTHER);
        particle.SetBoundary(i, ny - 1, OTHER);
    }                                                                                   //  Set boundary condition

    //--------------------Set design variable--------------------
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            gamma[particle.GetIndex(i, j)] = pow(i - 0.5*nx, 2.0) + pow(j, 2.0) < pow(0.15*nx, 2.0) ? 0.1 : 0.9;
        }
    }

    //--------------------Direct analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        EL::InitialCondition(i, particle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    EL::UpdateMacro(particle, rho, ux, uy, sxx, sxy, syx, syy);                         //  Update macroscopic values
    for (int t = 0; t < nt - 1; t++) {
        EL::Collision(elasticy, particle, rho, ux, uy, sxx, sxy, syx, syy);             //  Collision
        particle.Stream();                                                              //  Stream
        for (int i = 0; i < nx; i++) {
            particle.SetStress(i, 0, 0.0, 0.0);
            particle.SetStress(i, ny - 1, 0.0, 0.0);
        }
        for (int j = 0; j < ny; j++) {
            particle.SetRS(0, j);
            if (fabs(j - 0.5*ny) < 5) {
                particle.SetStress(nx - 1, j, 0.0, stress0);
            } else {
                particle.SetStress(nx - 1, j, 0.0, 0.0);
            }
        }                                                                               //  Boundary condition
        EL::UpdateMacro(particle, rho, ux, uy, sxx, sxy, syx, syy);                     //  Update macroscopic values
        for (int i = 0; i < nx*ny; i++) {
            rx[i] += ux[i];
            ry[i] += uy[i];
        }
    }

    //--------------------Invert analyze--------------------
    for (int i = 0; i < nx*ny; i++) {
        AEL::InitialCondition(i, particle, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }                                                                                   //  Set initial condition
    AEL::UpdateMacro(particle, irho, imx, imy, isxx, isxy, isyx, isyy);                 //  Update macroscopic values
    for (int t = tmax; t >= 0; t--) {
        AEL::Collision(elasticy, particle, irho, imx, imy, isxx, isxy, isyx, isyy);     //  Collision
        particle.iStream();                                                             //  Stream
        for (int i = 0; i < nx; i++) {
            particle.SetiStress(i, 0, rho0, 0.0, 0.0);
            particle.SetiStress(i, ny - 1, rho0, 0.0, 0.0);
        }
        for (int j = 0; j < ny; j++) {
            particle.SetiRS(0, j);
            if (fabs(j - 0.5*ny) < 5) {
                particle.SetiStress(nx - 1, j, rho0, 0.0, stress0);
            } else {
                particle.SetiStress(nx - 1, j, rho0, 0.0, 0.0);
            }
        }                                                                               //  Boundary condition                                                                          //  Boundary condition (inlet)
        AEL::UpdateMacro(particle, irho, imx, imy, isxx, isxy, isyx, isyy);             //  Update macroscopic values
    }

    //--------------------Get sensitivity--------------------
    double sensitivitymax = 0.0;
    for (int i = 0; i < nx*ny; i++) {
        sensitivity[i] = 13.5*(sxx[i]*isxx[i] + sxy[i]*isxy[i] + syx[i]*isyx[i] + syy[i]*isyy[i]) - 4.5*irho[i]*(sxx[i] + syy[i]);
        if (sensitivitymax < fabs(sensitivity[i])) {
            sensitivitymax = fabs(sensitivity[i]);
        }
    }
    for (int i = 0; i < nx*ny; i++) {  
        sensitivity[i] /= sensitivitymax;
    }

    //--------------------Export result--------------------
    VTKExport file("result/adjointelastic.vtk", nx, ny);
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointVector("r", 
        [&](int _i, int _j, int _k) { return rx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return ry[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointTensor("s", 
        [&](int _i, int _j, int _k) { return sxx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return sxy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [&](int _i, int _j, int _k) { return syx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return syy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("irho", [&](int _i, int _j, int _k) { return irho[particle.GetIndex(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return imx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return imy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointTensor("is", 
        [&](int _i, int _j, int _k) { return isxx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return isxy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [&](int _i, int _j, int _k) { return isyx[particle.GetIndex(_i, _j)]; },
        [&](int _i, int _j, int _k) { return isyy[particle.GetIndex(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("sensitivity", [&](int _i, int _j, int _k) { return sensitivity[particle.GetIndex(_i, _j)]; });

    return 0;
}