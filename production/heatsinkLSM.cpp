#define _USE_AVX_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/residual.h"
#include "../src/utility/reactiondiffusion.h"
#include "../src/utility/normalize.h"
#ifdef _USE_MPI_DEFINES
    #include "../src/utility/vtkxmlexport.h"
#else
    #include "../src/utility/vtkexport.h"
#endif

using namespace PANSLBM2;

int main(int argc, char** argv) {
#ifdef _USE_MPI_DEFINES
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int nPEx = atoi(argv[1]), nPEy = atoi(argv[2]);
    assert(nPEx*nPEy == PeTot);
#else
    int MyRank = 0, nPEx = 1, nPEy = 1;
#endif

    //********************Parameters********************
    const int lx = 141, ly = 161, mx = 81, my = 101, nt = 30000, dt = 100, nk = 50;
    double Pr = 6.0, Ra = 1e3, nu = 0.1, L = 4.0, tem0 = 0.0, qn = 1.0e-2;
    double movelimit = 0.02, weightlimit = 0.5, eps = 1.0e-6, s0 = -1.0, tau = 0.005;

    double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1);
    D2Q9<double> pf(lx, ly, MyRank, nPEx, nPEy), pg(lx, ly, MyRank, nPEx, nPEy);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pf.nxyz], *qy = new double[pf.nxyz], *qxp = new double[pf.nxyz], *qyp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqxp = new double[pg.nxyz], *iqyp = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0; ux[idx] = 0.0; uy[idx] = 0.0; uxp[idx] = 0.0; uyp[idx] = 0.0; tem[idx] = 0.0; qx[idx] = 0.0; qy[idx] = 0.0; qxp[idx] = 0.0; qyp[idx] = 0.0;
        irho[idx] = 0.0; iux[idx] = 0.0; iuy[idx] = 0.0; iuxp[idx] = 0.0; iuyp[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0; iqxp[idx] = 0.0; iqyp[idx] = 0.0;
    }
    double *chi = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dkds = new double[pf.nxyz];
    double *gi = new double[pg.nxyz*pg.nc], *igi = new double[pg.nxyz*pg.nc];
    
    if (MyRank == 0) {
        std::cout << "U:" << U << std::endl;
        std::cout << "gy:" << gy << std::endl;
    }
    
    std::vector<double> s(pf.nxyz, 1.0);
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            if ((i + pf.offsetx) < mx && (j + pf.offsety) < my) {
                int idx = pf.Index(i, j);
                s[idx] = s0;
            }
        }
    }

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    for (int k = 1; k <= nk; k++) {
        //********************Set parameter********************
        for (int idx = 0; idx < s.size(); idx++) {
            chi[idx] = s[idx] >= 0.0 ? 1.0 : 0.0;
            diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*chi[idx];
            dkds[idx] = diff_fluid - diff_solid;
        }
        
        //********************Constraint function********************
        double g_buffer = 0.0, g;
        std::vector<double> dgds(s.size(), 0.0);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                if ((i + pf.offsetx) < mx && (j + pf.offsety) < my) {
                    int idx = pf.Index(i, j);
                    g_buffer += (1.0 - chi[idx])/(weightlimit*mx*my);
                    dgds[idx] = -1.0/(weightlimit*mx*my); 
                }
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&g_buffer, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        g = g_buffer;
#endif
        g -= 1.0;
        
        //********************Objective function********************
        //  Direct analyse
        int td = nt;
        if (MyRank == 0) {
            std::cout << "Direct analyse t = 0";
        }
        NS::InitialCondition(pf, rho, ux, uy);
        AD::InitialCondition(pg, tem, ux, uy);
        for (int t = 1; t <= nt; t++) {
            AD::MacroCollideNaturalConvectionLSM(pf, rho, ux, uy, chi, nu, pg, tem, qx, qy, diffusivity, gx, gy, tem0, true, gi);
            if (t%dt == 0) {
                if (MyRank == 0) {
                    std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                }
                if (Residual(ux, uy, uxp, uyp, pf.nxyz) < eps && Residual(qx, qy, qxp, qyp, pf.nxyz) < eps) {
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
                
        //  Inverse analyse
        int ti = nt;
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
        AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
        for (int t = 1; t <= nt; t++) {
            AAD::MacroCollideNaturalConvectionLSM(
                pf, rho, ux, uy, irho, iux, iuy, imx, imy, chi, nu, 
                pg, tem, item, iqx, iqy, diffusivity,
                gx, gy, true, igi
            );
            if (t%dt == 0) {
                if (MyRank == 0) {
                    std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
                }
                if (Residual(iux, iuy, iuxp, iuyp, pf.nxyz) < eps && Residual(iqx, iqy, iqxp, iqyp, pf.nxyz) < eps) {
                    ti = t;
                    break;
                }
            }
            pf.iStream();
            pg.iStream();
            AAD::iBoundaryConditionSetT(pg, ux, uy, [=](int _i, int _j) { return _i == lx - 1 || _j == ly - 1; });
            AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return _j == 0; });
            AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return _j == 0 && _i < L; }, 1.0);
            pg.iBoundaryCondition([=](int _i, int _j) { return _i == 0 ? 2 : 0; });
            pf.iBoundaryCondition([=](int _i, int _j) { return _i == 0 ? 2 : 1; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            std::swap(iux, iuxp);
            std::swap(iuy, iuyp);
            std::swap(iqx, iqxp);
            std::swap(iqy, iqyp);
        }
               
        //  Get sensitivity
        double f_buffer = 0.0, f;
        for (int i = 0; i < pf.nx; ++i) {
            if ((i + pf.offsetx) < L && pf.PEy == 0) {
                f_buffer += tem[pf.Index(i, 0)];
            }
        }
        #ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        f /= L;
        std::vector<double> dfds(s.size(), 0.0);
        AAD::SensitivityTemperatureAtHeatSourceLSM(
            pg, dfds.data(), rho, ux, uy, iux, iuy, nu, tem, item, iqx, iqy, gi, igi, diffusivity, dkds, chi,
            [=](int _i, int _j) { return (_j == 0 && _i < L) ? qn : 0.0; }, 
            [=](int _i, int _j) { return _j == 0 && _i < L; }
        );
        Normalize(dfds.data(), pf.nxyz);

        //********************Update design variable********************
        if (MyRank == 0) {
            std::cout << "\rUpdate design variable" << std::string(10, ' ');
        }
        ReactionDiffusion::UpdateVariables(
            pf, s, f, dfds, g, dgds, tau, movelimit, 
            [&](int _i, int _j){ 
                return (_i + pf.offsetx) >= mx || (_j + pf.offsety) >= my; 
            }, 
            [](int _i, int _j){ 
                return 1.0;
            }
        );

        //********************Check convergence********************
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << k << " " << f << " " << g << " " << td << " " << ti << std::endl;
        }
        if (k == nk) {
#ifdef _USE_MPI_DEFINES
            VTKXMLExport file(pf, "result/heatsinkLSM");
            file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointData(pf, "u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointData(pf, "T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
            file.AddPointData(pf, "q", 
                [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) {   return irho[pf.Index(_i, _j)];  });
            file.AddPointData(pf, "iu", 
                [&](int _i, int _j, int _k) {   return iux[pf.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return iuy[pf.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointData(pf, "im", 
                [&](int _i, int _j, int _k) {   return imx[pf.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return imy[pf.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointData(pf, "iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j)];  });
            file.AddPointData(pf, "iq", 
                [&](int _i, int _j, int _k) {   return iqx[pg.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return iqy[pg.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)];    });
            file.AddPointData(pf, "chi", [&](int _i, int _j, int _k) { return chi[pf.Index(_i, _j)];    });
            file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) { return dfds[pf.Index(_i, _j)];    });
#else
            VTKExport file("result/heatsinkLSM.vtk", pf.nx, pf.ny);
            file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]; });
            file.AddPointVector("u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j)]; });
            file.AddPointVector("q", 
                [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j)]; },
                [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j)]; },
                [](int _i, int _j, int _k) { return 0.0; }
            );
            file.AddPointScaler("irho", [&](int _i, int _j, int _k) {   return irho[pf.Index(_i, _j)];  });
            file.AddPointVector("iu", 
                [&](int _i, int _j, int _k) {   return iux[pf.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return iuy[pf.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointVector("im", 
                [&](int _i, int _j, int _k) {   return imx[pf.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return imy[pf.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointScaler("iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j)];  });
            file.AddPointVector("iq", 
                [&](int _i, int _j, int _k) {   return iqx[pg.Index(_i, _j)];   },
                [&](int _i, int _j, int _k) {   return iqy[pg.Index(_i, _j)];   },
                [](int _i, int _j, int _k) {   return 0.0;   }
            );
            file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)];    });
            file.AddPointScaler("chi", [&](int _i, int _j, int _k) { return chi[pf.Index(_i, _j)];    });
            file.AddPointScaler("dfds", [&](int _i, int _j, int _k) { return dfds[pf.Index(_i, _j)];    });
#endif
            break;
        }
    }
    
    delete[] rho;   delete[] ux;    delete[] uy;    delete[] uxp;   delete[] uyp;
    delete[] tem;   delete[] qx;    delete[] qy;    delete[] qxp;   delete[] qyp;
    delete[] diffusivity;   delete[] chi; delete[] dkds;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] iuxp;  delete[] iuyp;  delete[] imx;   delete[] imy;
    delete[] item;  delete[] iqx;   delete[] iqy;   delete[] iqxp;  delete[] iqyp;
    delete[] gi;    delete[] igi;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}