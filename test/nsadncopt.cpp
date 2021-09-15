#define _USE_MPI_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.cuh"
#include "../src/equation/advection.cuh"
#include "../src/equation/adjointadvection.cuh"
#include "../src/prepost/vtkxmlexport.h"
#include "../utility/residual.cuh"
#include "../utility/mma.h"
#include "../utility/densityfilter.h"

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
    int lx = 141, ly = 161, mx = 81, my = 101, nt = 100000, dt = 100, nk = 2000, nb = 100;
    double Pr = 6.0, Ra = 2.5e3, nu = 0.1, L = 4.0, tem0 = 0.0, qn = 1.0e-2, alphamax = 1.0e4;
    double qf = 1e-2, qg = 1e0, movelimit = 0.2, weightlimit = 0.5, R = 1.5, eps = 1.0e-5, s0 = 0.5;

    double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1);
    D2Q9<double> pf(lx, ly, MyRank, nPEx, nPEy), pg(lx, ly, MyRank, nPEx, nPEy);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pf.nxyz], *qy = new double[pf.nxyz], *qxp = new double[pf.nxyz], *qyp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqxp = new double[pg.nxyz], *iqyp = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0;   ux[idx] = 0.0;    uy[idx] = 0.0;    uxp[idx] = 0.0;   uyp[idx] = 0.0;   tem[idx] = 0.0;   qx[idx] = 0.0;    qy[idx] = 0.0;    qxp[idx] = 0.0; qyp[idx] = 0.0;
        irho[idx] = 0.0;  iux[idx] = 0.0;   iuy[idx] = 0.0;   iuxp[idx] = 0.0;  iuyp[idx] = 0.0;  imx[idx] = 0.0;   imy[idx] = 0.0;   item[idx] = 0.0;  iqx[idx] = 0.0;   iqy[idx] = 0.0;  iqxp[idx] = 0.0;   iqyp[idx] = 0.0;
    }
    double *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dads = new double[pf.nxyz], *dkds = new double[pf.nxyz];
    double *gi = new double[pg.nxyz*pg.nc], *igi = new double[pg.nxyz*pg.nc];
    
    if (MyRank == 0) {
        std::cout << "U:" << U << std::endl;
        std::cout << "gy:" << gy << std::endl;
    }
    
    std::vector<double> s(pf.nxyz, 1.0), snm1(pf.nxyz, 1.0);
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            if ((i + pf.offsetx) < mx && (j + pf.offsety) < my) {
                int idx = pf.Index(i, j);
                s[idx] = s0;
                snm1[idx] = s0;
            }
        }
    }
    MMA<double> optimizer(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), 
        std::vector<double>(s.size(), 1.0)
    );
    optimizer.move = movelimit;

    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    for (int k = 1; k <= nk; k++) {
        if (k%nb == 0) {
            qf = std::min(1e7, qf*10.0);
        }

        //********************Filter variables********************
        std::vector<double> ss = DensityFilter::GetFilteredValue(pf, R, s);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                int idx = pf.Index(i, j);
                ss[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my) ? ss[idx] : 1.0;
            }
        }
        
        //********************Set alpha and diffusivity********************
        for (int idx = 0; idx < pf.nxyz; idx++) {
            diffusivity[idx] = diff_solid + (diff_fluid - diff_solid)*ss[idx]*(1.0 + qg)/(ss[idx] + qg);
            alpha[idx] = alphamax/(double)(ly - 1)*qf*(1.0 - ss[idx])/(ss[idx] + qf);
            dkds[idx] = (diff_fluid - diff_solid)*qg*(1.0 + qg)/pow(ss[idx] + qg, 2.0);
            dads[idx] = -alphamax/(double)(ly - 1)*qf*(1.0 + qf)/pow(ss[idx] + qf, 2.0);
        }
        
        //********************Constraint function********************
        double g_buffer = 0.0, g;
        std::vector<double> dgdss(s.size(), 0.0);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                if ((i + pf.offsetx) < mx && (j + pf.offsety) < my) {
                    int idx = pf.Index(i, j);
                    g_buffer += (1.0 - ss[idx])/(weightlimit*mx*my);
                    dgdss[idx] = -1.0/(weightlimit*mx*my); 
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
            AD::MacroBrinkmanCollideStreamNaturalConvection(pf, rho, ux, uy, alpha, pg, tem, qx, qy, gx, gy, tem0, nu, diffusivity, true);
            if (t%dt == 0) {
                if (MyRank == 0) {
                    std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                }
                if (Residual(ux, uy, uxp, uyp, pf.nxyz) < eps && Residual(qx, qy, qxp, qyp, pf.nxyz) < eps) {
                    td = t;
                    break;
                }
            }
            pf.Swap();
            pg.Swap();
            pf.BoundaryCondition([=]__device__(int _i, int _j) { return _i == 0 ? 2 : 1; });
            AD::BoundaryConditionSetT(pg, 
                [=]__device__(int _i, int _j) { return tem0; }, 
                ux, uy, 
                [=]__device__(int _i, int _j) { return _i == lx - 1 || _j == ly - 1; }
            );
            AD::BoundaryConditionSetFlux(pg, 
                [=]__device__(int _i, int _j) { return (_j == 0 && _i < L) ? qn : 0.0; }, 
                ux, uy, 
                [=]__device__(int _i, int _j) { return _j == 0; }, 
                diffusivity
            );
            pg.BoundaryCondition([=]__device__(int _i, int _j) { return _i == 0 ? 2 : 0; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            std::swap(ux, uxp);
            std::swap(uy, uyp);
            std::swap(qx, qxp);
            std::swap(qy, qyp);
        }
        for (int idxc = 0; idxc < pg.nxyz*pg.nc; ++idxc) {
            gi[idxc] = pg.f[idxc];
        }
                
        //  Inverse analyse
        int ti = nt;
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, rho, ux, uy, irho, iux, iuy);
        AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
        for (int t = 1; t <= nt; t++) {
            AAD::MacroBrinkmanCollideStreamNaturalConvection(
                pf, rho, ux, uy, irho, iux, iuy, imx, imy, alpha, 
                pg, tem, item, iqx, iqy,
                gx, gy, nu, diffusivity, true
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
            pf.Swap();
            pg.Swap();
            AAD::BoundaryConditionSetiT(pg, ux, uy, [=]__device__(int _i, int _j) { return _i == lx - 1 || _j == ly - 1; });
            AAD::BoundaryConditionSetiFlux(pg, ux, uy, [=]__device__(int _i, int _j) { return _j == 0; });
            AAD::BoundaryConditionSetiFlux(pg, ux, uy, [=]__device__(int _i, int _j) { return _j == 0 && _i < L; }, 1.0);
            pg.iBoundaryCondition([=]__device__(int _i, int _j) { return _i == 0 ? 2 : 0; });
            pf.iBoundaryCondition([=]__device__(int _i, int _j) { return _i == 0 ? 2 : 1; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            std::swap(iux, iuxp);
            std::swap(iuy, iuyp);
            std::swap(iqx, iqxp);
            std::swap(iqy, iqyp);
        }
        for (int idxc = 0; idxc < pg.nxyz*pg.nc; ++idxc) {
            igi[idxc] = pg.f[idxc];
        }
               
        //  Get sensitivity
        double f_buffer = 0.0, f;
        for (int i = 0; i < nx; ++i) {
            if ((i + pf.offsetx) < L) {
                f_buffer += tem[pf.Index(i, 0)];
            }
        }
        #ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        f /= L;
        std::vector<double> dfdss(s.size(), 0.0);
        AAD::SensitivityTemperatureAtHeatSource(pg, 
            [=]__device__(int _i, int _j) { return (_j == 0 && _i < L) ? qn : 0.0; }, 
            ux, uy, imx, imy, tem, item, iqx, iqy, gi, igi, dfds, diffusivity, dads, dkds,
            [=]__device__(int _i, int _j) { return _j == 0 && _i < L; }
        );
        double dfdsmax_buffer = 0.0, dfdsmax;
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            if (dfdssmax_buffer < fabs(dfdss[idx])) {
                dfdssmax_buffer = fabs(dfdss[idx]);
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&dfdsmax_buffer, &dfdsmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
        dfdsmax = dfdsmax_buffer;
#endif
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            dfdss[idx] /= dfdssmax;
        }
        
        //********************Filter sensitivities********************
        std::vector<double> dfds = DensityFilter::GetFilteredValue(pf, R, dfdss);
        std::vector<double> dgds = DensityFilter::GetFilteredValue(pf, R, dgdss);
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int idx = pf.Index(i, j);
                dfds[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my) ? dfds[idx] : 0.0;
                dgds[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my) ? dgds[idx] : 0.0;
            }
        }

        //********************Check grayscale********************
        double mnd_buffer = 0.0, mnd;
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                mnd_buffer += ss[pf.Index(i, j)]*(1.0 - ss[pf.Index(i, j)]);
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&mnd_buffer, &mnd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        mnd = mnd_buffer;
#endif
        mnd *= 4.0/(double)(mx*my)*100.0;

        //********************Update design variable********************
        if (MyRank == 0) {
            std::cout << "\rUpdate design variable" << std::string(10, ' ');
        }
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
        double dsmax = 0.0;
        int imax = 0, jmax = 0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                int idx = pf.Index(i, j);
                s[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my) ? s[idx] : 1.0;
                double tmpds = fabs(s[idx] - snm1[idx]);
                if (dsmax < tmpds) {
                    dsmax = tmpds;
                    imax = i;
                    jmax = j;
                }
                snm1[idx] = s[idx];
            }
        }

        //********************Check convergence********************
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << k << "\t" << f << "\t" << g << "\t" << td << "\t" << ti << "\t" << dsmax  << " (" << imax << "," << jmax << ")\t" << qf << "\t" << qg << "\t" << mnd << std::endl;
        }
        if ((k > 1 && dsmax < 0.01 && g <= 0.0) || k == nk) {
            if (MyRank == 0) {
                std::cout << "----------Convergence/Last step----------" << std::endl;
            }

            VTKXMLExport file(pf, "result/nsadncopt");
            file.AddPointData(pf, "k", [&](int _i, int _j, int _k) { return diffusivity[pg.Index(_i, _j)]; });
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
            file.AddPointData(pf, "ss", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j)];    });
            file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) { return dfds[pf.Index(_i, _j)];    });
            file.AddPointData(pf, "dkds", [&](int _i, int _j, int _k) { return dkds[pf.Index(_i, _j)];    });
            file.AddPointData(pf, "dads", [&](int _i, int _j, int _k) { return dads[pf.Index(_i, _j)];    });
            
            for (int i = 0; i < nx; ++i) {
                for (int j = 0; j < ny; ++j) {
                    int idx = pf.Index(i, j);
                    ss[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && ss[idx] < 0.1) ? 0.0 : 1.0;
                }
            }
            file.AddPointData(pf, "ss2", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j)];    });
            break;
        }
    }
    
    delete[] rho, ux, uy, uxp, uyp, tem, qx, qy, qxp, qyp, diffusivity, alpha, dkds, dads, irho, iux, iuy, iuxp, iuyp, imx, imy, item, iqx, iqy, iqxp, iqyp, gi, igi;
}