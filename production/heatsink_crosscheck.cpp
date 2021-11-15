#define _USE_MATH_DEFINES
#define _USE_AVX_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/utility/residual.h"
#include "../src/utility/vtkimport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    //--------------------Set parameters--------------------
    const int dt = 100, nt0 = 0, nc = 3, nm = 3;
    double Pr = 6.0, nu = 0.1, tem0 = 0.0, qn = 1.0e-2, alphamax = 1.0e4, eps = 1.0e-5;
    int ntList[nc] = { 1000000, 1000000, 1000000 };
    double RaList[nc] = { 2.5e3, 2.5e4, 2.5e5 }, qfList[nc] = { 1e1, 1e1, 1e1 }, qgList[nc] = { 1e0, 1e0, 1e0 }, fList[nc*nm] = { 0 };

    if (argc != nm + 1) {
        std::cout << "Error:No vtk file selected." << std::endl;
        exit(1);
    }
    
    //  Loop of model
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "Model id:" << modelid << " and fname:" << argv[modelid + 1] << std::endl;
        //  Open model file
        std::string fname(argv[modelid + 1]);
        VTKImport model(fname);
        int lx = model.GetNx(), ly = model.GetNy(), L = (ly - 1)/40;
        D2Q9<double> pf(lx, ly), pg(lx, ly);
        double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
        double *tem = new double[pg.nxyz], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz], *qxp = new double[pg.nxyz], *qyp = new double[pg.nxyz];
        double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz];
        model.GetPointScalar("ss", s);
        
        //  Loop of condition
        for (int conditionid = 0; conditionid < nc; ++conditionid) {
            int nt = ntList[conditionid];
            double Ra = RaList[conditionid], qf = qfList[conditionid], qg = qgList[conditionid];
            double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1);

            for (int idx = 0; idx < pf.nxyz; idx++) {
                rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0; uxp[idx] = 0.0; uyp[idx] = 0.0; tem[idx] = 0.0; qx[idx] = 0.0; qy[idx] = 0.0; qxp[idx] = 0.0; qyp[idx] = 0.0;
                diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*s[idx]*(1.0 + qg)/(s[idx] + qg);
                alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - s[idx])/(s[idx] + qf);
            }
            std::cout << "\rCondition id:" << conditionid << " t = 0";
            NS::InitialCondition(pf, rho, ux, uy);
            AD::InitialCondition(pg, tem, ux, uy);
            int td = nt;
            for (int t = 1; t < nt0 + nt; ++t) {
                AD::MacroBrinkmanCollideNaturalConvection(pf, rho, ux, uy, alpha, nu, pg, tem, qx, qy, diffusivity, gx, gy, tem0, true);
                if (t%dt == 0) {
                    double du = Residual(ux, uy, uxp, uyp, pf.nxyz), dq = Residual(qx, qy, qxp, qyp, pf.nxyz);
                    std::cout << "\rCondition id:" << conditionid << " t = " << t << " du = " << du << " dq = " << dq << std::string(10, ' ');
                    if (du < eps && dq < eps) {
                        td = t;
                        break;
                    }
                }
                pf.Stream();
                pg.Stream();
                pf.BoundaryCondition([=](int _i, int _j) { return _i == 0 ? 2 : 1; });
                AD::BoundaryConditionSetT(pg, 
                    [=](int _i, int _j) { return tem0; }, 
                    ux, uy, 
                    [=](int _i, int _j) { return _i == lx - 1 || _j == ly - 1; }
                );
                AD::BoundaryConditionSetQ(pg, 
                    [=](int _i, int _j) { return (_j == 0 && _i < L) ? qn : 0.0; }, 
                    ux, uy, diffusivity,
                    [=](int _i, int _j) { return _j == 0; } 
                );
                pg.BoundaryCondition([=](int _i, int _j) { return _i == 0 ? 2 : 0; });
                pf.SmoothCorner();
                pg.SmoothCorner();

                std::swap(ux, uxp);
                std::swap(uy, uyp);
                std::swap(qx, qxp);
                std::swap(qy, qyp);
            }
            std::cout << "\rCondition id:" << conditionid << " td = " << td << std::string(50, ' ') << std::endl;
            for (int i = 0; i < pf.nx; ++i) {
                if ((i + pf.offsetx) < L && pf.PEy == 0) {
                    fList[conditionid*nm + modelid] += tem[pf.Index(i, 0)];
                }
            }
            fList[conditionid*nm + modelid] /= (double)L;
        }

        delete[] rho; delete[] ux; delete[] uy; delete[] uxp; delete[] uyp; 
        delete[] tem; delete[] qx; delete[] qy; delete[] qxp; delete[] qyp; 
        delete[] s; delete[] alpha; delete[] diffusivity;
    }

    std::cout << std::endl << "**********objective**********" << std::endl << "cnd/mod\t";
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "\t" << modelid;    
    }
    for (int conditionid = 0; conditionid < nc; ++conditionid) {
        std::cout << std::endl << std::scientific << std::setprecision(2) << conditionid << std::setprecision(6);
        for (int modelid = 0; modelid < nm; ++modelid) {
            std::cout << "\t" << fList[conditionid*nm + modelid];
        }    
    }

    std::cout << std::endl;
}