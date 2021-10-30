#define _USE_MATH_DEFINES
#define _USE_AVX_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/utility/vtkimport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    //--------------------Set parameters--------------------
    const int dt = 1000, nt0 = 1000000, nc = 3, nm = 3;
    int ntList[nc] = { 50000, 50000, 50000 }, periodList[nc] = { 50000, 50000, 50000 }, dutyList[nc] = { 20, 50, 80 };
    double viscosity = 0.1/6.0, diff_fluid = viscosity/1.0, Th = 1.0, Tl = 0.0, gx = 0.0, gy = 1/3.6e6;//1/4.5e5;
    double alphamax = 1e5, diff_solid = diff_fluid*10.0;
    double qfList[nc] = { 1e-5, 1e-5, 1e-5 }, qgList[nc] = { 1e-4, 1e-4, 1e-4 }, fList[nc*nm] = { 0 };

    if (argc != nm + 1) {
        std::cout << "Error:No vtk file selected." << std::endl;
        exit(1);
    }
        
    //auto tembc = [=](int _t, int _period, int _duty) { return Th*(1 - cos(2*M_PI*_t/_period)); };
    auto tembc = [=](int _t, int _period, int _duty) { return _t%_period < _period*_duty/100.0 ? (Th - Tl)*100.0/(double)_duty + Tl : Tl; };

    //  Loop of model
    for (int modelid = 0; modelid < nm; ++modelid) {
        std::cout << "\r" << std::string(80, ' ') << std::endl << "Model id:" << modelid << "\tand fname:" << argv[modelid + 1] << std::endl;
        //  Open model file
        std::string fname(argv[modelid + 1]);
        VTKImport model(fname);
        int lx = model.GetNx(), ly = model.GetNy();
        D2Q9<double> pf(lx, ly), pg(lx, ly);
        double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz];
        double *tem = new double[pg.nxyz], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz];
        double *s = new double[pf.nxyz], *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz];
        double *directionx = new double[pf.nxyz], *directiony = new double[pf.nxyz];
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                int idx = pf.Index(i, j);
                directionx[idx] = ((i + pf.offsetx) == lx/2 && (j + pf.offsety) > 9*ly/10) ? -1.0 : 0.0;
                directiony[idx] = 0.0;
            }
        }
        model.GetPointScalar("ss", s);

        //  Loop of condition
        for (int conditionid = 0; conditionid < nc; ++conditionid) {
            int nt = ntList[conditionid], period = periodList[conditionid], duty = dutyList[conditionid];
            double qf = qfList[conditionid], qg = qgList[conditionid];
            std::cout << "\rCondition id:" << conditionid << std::string(50, ' ') << std::endl;
            for (int idx = 0; idx < pf.nxyz; idx++) {
                rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0; tem[idx] = 0.5*(Tl + Th); qx[idx] = 0.0; qy[idx] = 0.0;
                diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*s[idx]*(1.0 + qg)/(s[idx] + qg);
                alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - s[idx])/(s[idx] + qf);
            }
            std::cout << "Direct analyse t = 0";
            NS::InitialCondition(pf, rho, ux, uy);
            AD::InitialCondition(pg, tem, ux, uy);
            for (int t = 1; t < nt0 + nt; ++t) {
                if (t%dt == 0) {
                    std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                }
                AD::MacroBrinkmanCollideNaturalConvection(pf, rho, ux, uy, alpha, viscosity, pg, tem, qx, qy, diffusivity, gx, gy, 0.5*(Th + Tl), true);

                pf.Stream();
                pg.Stream();
                pf.BoundaryCondition([=](int _i, int _j) { return 1; });
                pg.BoundaryCondition([=](int _i, int _j) { return 0; });
                AD::BoundaryConditionSetT(pg, 
                    [=](int _i, int _j) { return _i == 0 ? tembc(t, period, duty) : Tl; }, 
                    ux, uy, 
                    [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); }
                );
                AD::BoundaryConditionSetQ(pg, 
                    [=](int _i, int _j) { return 0.0; },
                    ux, uy, diffusivity, 
                    [=](int _i, int _j) { return (_i == 0 && ly/2 <= _j) || (_i == lx - 1 && ly/2 <= _j) || _j == 0 || _j == ly - 1; }
                );
                pf.SmoothCorner();
                pg.SmoothCorner();

                pf.BoundaryConditionAlongXEdge(lx/5, 1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                pf.BoundaryConditionAlongYEdge(ly/2, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                pf.BoundaryConditionAlongXEdge(4*lx/5, -1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                pf.BoundaryConditionAlongYEdge(9*ly/10, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                pf.SmoothCornerAt(lx/5, ly/2, -1, -1);
                pf.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
                pf.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
                pf.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
                AD::BoundaryConditionSetQAlongXEdge(pg, lx/5, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                AD::BoundaryConditionSetQAlongYEdge(pg, ly/2, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                AD::BoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, [=](int _i, int _j) { return 0.0; }, ux, uy, diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                AD::BoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, [=](int _i, int _j) { return 0.0; }, ux, uy, diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
                pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
                pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
                pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);

                if (t > nt0) {
                    for (int j = 0; j < pf.ny; ++j) {
                        int i = lx/2 - pf.offsetx;
                        if (0 <= i && i < pf.nx && (j + pf.offsety) > 9*ly/10) {
                            int idx = pf.Index(i, j);
                            fList[conditionid*nm + modelid] += ux[idx]*directionx[idx] + uy[idx]*directiony[idx];
                        }
                    }
                }
            }
            fList[conditionid*nm + modelid] /= (double)(ly/10)*nt;
        }

        delete[] rho; delete[] ux; delete[] uy; delete[] tem; delete[] qx; delete[] qy; delete[] s; delete[] alpha; delete[] diffusivity; delete[] directionx; delete[] directiony;
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
}