#define _USE_AVX_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <iomanip>
#include <cmath>

#include "../src/particle/d3q15.h"
#include "../src/equation/advection.h"
#include "../src/utility/vtkimport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    assert(argc == 2);

    //********************Parameters********************
    int nt = 60000, dt = 100;
    double Pr = 6.0, Ra = 2e5, nu = 0.1, tem0 = 0.0, qn = 1.0e-2, alphamax = 1.0e4, qf = 1e1, qg = 1e0, weightlimit = 0.05, sth = 0.5;

    VTKImport model(argv[1]);
    int lx = model.GetNx(), ly = model.GetNy(), lz = model.GetNz(), mx = 3*(lx - 1)/4 + 1, my = 3*(ly - 1)/4 + 1, mz = 3*(lz - 1)/4 + 1;
    double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1), gz = 0.0, L = (lx - 1)/10;
    std::cout << "U:" << U << std::endl;
    std::cout << "gy:" << gy << std::endl;

    D3Q15<double> pf(lx, ly, lz), pg(lx, ly, lz);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz], *qz = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0; uz[idx] = 0.0; tem[idx] = 0.0; qx[idx] = 0.0; qy[idx] = 0.0; qz[idx] = 0.0;
    }
    
    std::vector<double> s(pf.nxyz, 1.0);
    model.GetPointScalar("ss", s.data());
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        //s[idx] = s[idx] < sth ? 0.0 : 1.0;
    }
    double *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*s[idx]*(1.0 + qg)/(s[idx] + qg);
        alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - s[idx])/(s[idx] + qf);
    }

    //********************Constraint function********************
    double g = -1.0;
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            for (int k = 0; k < pf.nz; ++k) {
                if ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) {
                    int idx = pf.Index(i, j, k);
                    g += (1.0 - s[idx])/(weightlimit*mx*my*mz);
                }
            } 
        }
    }

    //********************Objective function********************
    double f = 0.0;
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();   
    NS::InitialCondition(pf, rho, ux, uy, uz);
    AD::InitialCondition(pg, tem, ux, uy, uz);
    for (int t = 1; t < nt; ++t) {
        AD::MacroBrinkmanCollideNaturalConvection(
            pf, rho, ux, uy, uz, alpha, nu, 
            pg, tem, qx, qy, qz, diffusivity, 
            gx, gy, gz, tem0, true
        );
        if (t%dt == 0) {
            std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
        }
        pf.Stream();
        pg.Stream(30);
        pf.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 1; });
        AD::BoundaryConditionSetT(pg, 
            [=](int _i, int _j, int _k) { return tem0; }, 
            ux, uy, uz,
            [=](int _i, int _j, int _k) { return _i == lx - 1 || _j == ly - 1 || _k == lz - 1; }
        );
        AD::BoundaryConditionSetQ(pg, 
            [=](int _i, int _j, int _k) { return (_j == 0 && _i < L && _k < L) ? qn : 0.0; }, 
            ux, uy, uz, diffusivity,
            [=](int _i, int _j, int _k) { return _j == 0; }
        );
        pg.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 0; });
        pf.SmoothCorner();
        pg.SmoothCorner();

        for (int k = 0; k < L; ++k) {
            for (int i = 0; i < L; ++i) {
                f += tem[pf.Index(i, 0, k)];
            }
        }
    }
    f /= (double)(L*L*nt);
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();

    //********************Export result********************
    std::cout << std::fixed << std::setprecision(6) << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << " " << f << " " << g << std::endl;

    delete[] rho; delete[] ux; delete[] uy; delete[] uz; delete[] tem; delete[] qx; delete[] qy; delete[] qz; delete[] diffusivity; delete[] alpha;

    return 0;
}