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
    //--------------------Set parameters--------------------
    const int dt = 100, nt0 = 0, nc = 3, nm = 3;
    double Pr = 6.0, nu = 0.1, tem0 = 0.0, qn = 1.5e-2, alphamax = 1.0e4, weightlimit = 0.05, sth = 0.1;
    int ntList[nc] = { 50000, 100000, 100000 };
    double RaList[nc] = { 1e5, 5e4, 1e5 }, qfList[nc] = { 1e-1, 1e-1, 1e-1 }, qgList[nc] = { 1e0, 1e0, 1e0 }, fList[nc][nm] = { 0 }, gList[nm] = { 0 };
    
    if (argc != nm + 1) {
        std::cout << "Error:No vtk file selected." << std::endl;
        exit(1);
    }

    //  Loop of model
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "Model id:" << modelid << " and fname:" << argv[modelid + 1] << std::endl;
        VTKImport model(argv[modelid + 1]);
        int lx = model.GetNx(), ly = model.GetNy(), lz = model.GetNz(), mx = 3*(lx - 1)/4 + 1, my = 3*(ly - 1)/4 + 1, mz = 3*(lz - 1)/4 + 1, L = (lx - 1)/10;
        D3Q15<double> pf(lx, ly, lz), pg(lx, ly, lz);
        double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz];
        double *tem = new double[pg.nxyz], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz], *qz = new double[pg.nxyz];
        double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz];
        model.GetPointScalar("ss", s);
        gList[modelid] = -1.0;
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            //s[idx] = s[idx] < sth ? 0.0 : 1.0;
            gList[modelid] += (1.0 - s[idx])/(weightlimit*mx*my);
        }

        //  Loop of condition
        for (int conditionid = 0; conditionid < nc; ++conditionid) {
            int nt = ntList[conditionid];
            double Ra = RaList[conditionid], qf = qfList[conditionid], qg = qgList[conditionid];
            double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1), gz = 0.0;
            for (int idx = 0; idx < pf.nxyz; idx++) {
                rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0; uz[idx] = 0.0; tem[idx] = 0.0; qx[idx] = 0.0; qy[idx] = 0.0; qz[idx] = 0.0;
                diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*s[idx]*(1.0 + qg)/(s[idx] + qg);
                alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - s[idx])/(s[idx] + qf);
            }
            std::cout << "\rCondition id:" << conditionid << " t = 0";
            NS::InitialCondition(pf, rho, ux, uy, uz);
            AD::InitialCondition(pg, tem, ux, uy, uz);
            for (int t = 1; t < nt; ++t) {
                AD::MacroBrinkmanCollideNaturalConvection(
                    pf, rho, ux, uy, uz, alpha, nu, 
                    pg, tem, qx, qy, qz, diffusivity, 
                    gx, gy, gz, tem0, true
                );
                if (t%dt == 0) {
                    std::cout << "\rCondition id:" << conditionid << " t = " << t << std::string(10, ' ');
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
                        fList[conditionid][modelid] += tem[pf.Index(i, 0, k)];
                    }
                }
            }
            std::cout << "\rCondition id:" << conditionid << std::string(50, ' ') << std::endl;
            fList[conditionid][modelid] /= (double)(L*L*nt);
        }

        delete[] rho; delete[] ux; delete[] uy; delete[] uz; delete[] tem; delete[] qx; delete[] qy; delete[] qz; 
        delete[] s; delete[] diffusivity; delete[] alpha;
    }

    std::cout << std::endl << "**********objective**********" << std::endl << "cnd/mod\t";
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "\t" << modelid;    
    }
    for (int conditionid = 0; conditionid < nc; ++conditionid) {
        std::cout << std::endl << std::scientific << std::setprecision(2) << conditionid << std::setprecision(6);
        for (int modelid = 0; modelid < nm; ++modelid) {
            std::cout << "\t" << fList[conditionid][modelid];
        }    
    }
    std::cout << std::endl << "g";
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "\t" << gList[modelid];    
    }

    std::cout << std::endl;
    return 0;
}