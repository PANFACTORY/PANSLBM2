#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/ib_navierstokes.h"
#include "../src/equation/advection.h"
#include "../src/equation/ib_advection.h"
#include "vtkexport.h"

using namespace PANSLBM2;

int main() {
    //--------------------Set parameters--------------------
    int tmax = 100000, nx = 100, ny = 50, nb = 64;
    double nu = 0.02, alpha = 0.02, Th = 2.0, Tl = 1.0, q0, Td = 1.7, v0 = 0.0, dv = 62.8/(double)nb;
    double rho[nx*ny], ux[nx*ny], uy[nx*ny], tem[nx*ny], qx[nx*ny], qy[nx*ny];
    double uxl[nx*ny], uyl[nx*ny], gx[nx*ny], gy[nx*ny], teml[nx*ny], gtem[nx*ny];
    double bpx[nb], bpy[nb], bux[nb], buy[nb], bgx[nb], bgy[nb], bvx[nb], bvy[nb], btem[nb], bgtem[nb], btemd[nb];

    D2Q9<double> particlef(nx, ny);
    D2Q9<double> particleg(nx, ny);
    for (int j = 0; j < ny; j++) {
        particlef.SetBoundary(0, j, PERIODIC);
        particlef.SetBoundary(nx - 1, j, PERIODIC);
        particleg.SetBoundary(0, j, PERIODIC);
        particleg.SetBoundary(nx - 1, j, PERIODIC);
    }
    for (int i = 0; i < nx; i++) {
        particlef.SetBoundary(i, 0, MIRROR);
        particlef.SetBoundary(i, ny - 1, MIRROR);
        particleg.SetBoundary(i, 0, OTHER);
        particleg.SetBoundary(i, ny - 1, OTHER);
    }                                               //  Set boundary condition
    for (int i = 0; i < nx*ny; i++) {
        NS::InitialCondition(i, particlef, 1.0, 0.0, 0.0);
        AD::InitialCondition(i, particleg, 0.5*(Th + Tl), 0.0, 0.0);
    }                                               //  Set initial condition
    
    for (int k = 0; k < nb; k++) {
        bpx[k] = 10.0*cos(2.0*M_PI*k/(double)nb) + (double)nx/2.0;
        bpy[k] = 10.0*sin(2.0*M_PI*k/(double)nb) + (double)ny/2.0;
        bvx[k] = -v0*sin(2.0*M_PI*k/(double)nb);
        bvy[k] = v0*cos(2.0*M_PI*k/(double)nb);
        btemd[k] = Td;
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    //--------------------Direct analyze--------------------
    NS::UpdateMacro(particlef, rho, ux, uy);
    AD::UpdateMacro(particleg, tem, qx, qy, ux, uy);      //  Update macroscopic values
    for (int t = 0; t < tmax; t++) {
        NS::Collision(nu, particlef, rho, ux, uy);
        AD::Collision(alpha, particleg, tem, ux, uy);     //  Collision
        particlef.Stream();
        particleg.Stream();                             //  Stream
        for (int i = 0; i < nx; i++) {
            particleg.SetFlux(i, 0, ux[particlef.GetIndex(i, 0)], uy[particlef.GetIndex(i, 0)], q0);
            particleg.SetFlux(i, ny - 1, ux[particlef.GetIndex(i, ny - 1)], uy[particlef.GetIndex(i, ny - 1)], q0);
        }                                               //  Boundary condition (Fix temperature)
        particlef.SmoothCorner();
        particleg.SmoothCorner();
        NS::ExternalForceNaturalConvection(0.0, 1.6e-5, 0.5*(Th + Tl), particlef, particleg);   //  External force with natural convection
        NS::UpdateMacro(particlef, rho, ux, uy);
        AD::UpdateMacro(particleg, tem, qx, qy, ux, uy);  //  Update macroscopic values
        NS::ExternalForceIB(particlef, ux, uy, uxl, uyl, gx, gy, nb, bpx, bpy, bux, buy, bgx, bgy, bvx, bvy, dv);
        AD::ExternalForceIBThermal(particleg, tem, teml, gtem, nb, bpx, bpy, btem, bgtem, btemd, dv);
        NS::UpdateMacro(particlef, rho, ux, uy);
        AD::UpdateMacro(particleg, tem, qx, qy, ux, uy);  //  Update macroscopic values

        if (t%1000 == 0) {
            std::cout << t/1000 << std::endl;
            VTKExport file("result/thermal" + std::to_string(t/1000) + ".vtk", nx, ny);
            file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[particlef.GetIndex(_i, _j)]/3.0; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[particlef.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[particlef.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[particleg.GetIndex(_i, _j)]; });
            file.AddPointVector("q", 
                [&](int _i, int _j, int _k) { return qx[particlef.GetIndex(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[particlef.GetIndex(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("boundaryf", [&](int _i, int _j, int _k) { return particlef.GetBoundary(_i, _j); });
            file.AddPointScaler("boundaryg", [&](int _i, int _j, int _k) { return particleg.GetBoundary(_i, _j); });
        } 
    }                                                   //  Export result per 1000 steps 

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}