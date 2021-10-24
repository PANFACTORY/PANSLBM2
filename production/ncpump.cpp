#define _USE_MATH_DEFINES
//#define _USE_MPI_DEFINES
#define _USE_AVX_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/vtkxmlexport.h"
#include "../src/utility/residual.h"
#include "../src/utility/mma.h"
#include "../src/utility/densityfilter.h"
#include "../src/utility/normalize.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
#ifdef _USE_MPI_DEFINES
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int mx = atoi(argv[1]), my = atoi(argv[2]);
    assert(mx*my == PeTot);
#else
    int MyRank = 0, mx = 1, my = 1; 
#endif

    //--------------------Set parameters--------------------
    int lx = 51, ly = 101, nt0 = 0, dt = 1000, nt = 50000, period = 5000, nk = 2000, nb = 100;
    double viscosity = 0.1, diff_fluid = viscosity/6.0, Th = 1.0, Tl = 0.0, gx = 0.0, gy = 1e-4;
    double alphamax = 1e4, diff_solid = diff_fluid*10.0, qf = 1e-2, qg = 1e0, weightlimit = 0.5, movelimit = 0.2, R = 1.5, eps = 1e-5;
    D2Q9<double> pf(lx, ly, MyRank, mx, my), pg(lx, ly, MyRank, mx, my);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pf.nxyz], *qy = new double[pf.nxyz], *qxp = new double[pf.nxyz], *qyp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqxp = new double[pg.nxyz], *iqyp = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0;  ux[idx] = 0.0;  uy[idx] = 0.0;  uxp[idx] = 0.0;  uyp[idx] = 0.0;  tem[idx] = Tl;  qx[idx] = 0.0;  qy[idx] = 0.0;   qxp[idx] = 0.0; qyp[idx] = 0.0;
        irho[idx] = 0.0; iux[idx] = 0.0; iuy[idx] = 0.0; iuxp[idx] = 0.0; iuyp[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0; iqxp[idx] = 0.0; iqyp[idx] = 0.0;
    }
    double *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dads = new double[pf.nxyz], *dkds = new double[pf.nxyz];
    double *gi = new double[pg.nxyz*pg.nc], *igi = new double[pg.nxyz*pg.nc];
    double *directionx = new double[pf.nxyz], *directiony = new double[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            directionx[idx] = ((i + pf.offsetx) == lx/2 && (j + pf.offsety) > 9*ly/10) ? -1.0 : 0.0;
            directiony[idx] = 0.0;
        }
    }

    std::vector<double> s(pf.nxyz, 1.0), snm1(pf.nxyz, 1.0);
    MMA<double> optimizer(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), 
        std::vector<double>(s.size(), 1.0)
    );
    optimizer.move = movelimit;

    for (int k = 1; k <= nk; k++) {
        if (k%nb == 0) {
            qf = std::min(1e7, qf*10.0);
        }

        //********************Filter variables********************
        std::vector<double> ss = DensityFilter::GetFilteredValue(pf, R, s);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                int idx = pf.Index(i, j);
                ss[idx] = (j + pf.offsety) < ly/2 ? ss[idx] : 1.0;
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
                if ((j + pf.offsety) < ly/2) {
                    int idx = pf.Index(i, j);
                    g_buffer += ss[idx]/(weightlimit*lx*ly/2);
                    dgdss[idx] = 1.0/(weightlimit*lx*ly/2); 
                }
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&g_buffer, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        g = g_buffer;
#endif
        g -= 1.0;

        //********************Direct analyze********************
        int td = nt;
        if (MyRank == 0) {
            std::cout << "Direct analyse t = 0";
        }
        NS::InitialCondition(pf, rho, ux, uy);
        AD::InitialCondition(pg, tem, ux, uy);
        for (int t = 1; t <= nt0 + nt; ++t) {
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                if (Residual(ux, uy, uxp, uyp, pf.nxyz) < eps && Residual(qx, qy, qxp, qyp, pf.nxyz) < eps) {
                    td = t;
                    break;
                }
            }
            AD::MacroBrinkmanCollideNaturalConvection(
                pf, rho, ux, uy, alpha, viscosity, 
                pg, tem, qx, qy, diffusivity, gx, gy, Tl, true, gi
            );

            pf.Stream();
            pg.Stream();
            pf.BoundaryCondition([=](int _i, int _j) { return 1; });
            pg.BoundaryCondition([=](int _i, int _j) { return 0; });
            AD::BoundaryConditionSetT(pg, 
                [=](int _i, int _j) { return _i == 0 ? Th : Tl; }, 
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

            std::swap(ux, uxp);
            std::swap(uy, uyp);
            std::swap(qx, qxp);
            std::swap(qy, qyp);
        }

        //********************Inverse analyze********************
        int ti = nt;
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
        AAD::InitialCondition(pg, ux, uy, item, iqx, iqy);
        for (int t = 1; t <= nt0 + nt; ++t) {
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
                if (Residual(iux, iuy, iuxp, iuyp, pf.nxyz) < eps && Residual(iqx, iqy, iqxp, iqyp, pf.nxyz) < eps) {
                    ti = t;
                    break;
                }
            }
            AAD::MacroBrinkmanCollideNaturalConvectionMassFlow(
                pf, rho, ux, uy, irho, iux, iuy, imx, imy, alpha, viscosity,
                pg, tem, item, iqx, iqy, diffusivity, gx, gy, 
                directionx, directiony, true, igi
            );

            pf.iStream();
            pg.iStream();
            pf.iBoundaryCondition([=](int _i, int _j) { return 1; });
            pg.iBoundaryCondition([=](int _i, int _j) { return 0; });
            AAD::iBoundaryConditionSetT(pg, ux, uy, [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); });
            AAD::iBoundaryConditionSetQ(pg, ux, uy, [=](int _i, int _j) { return (_i == 0 && ly/2 <= _j) || (_i == lx - 1 && ly/2 <= _j) || _j == 0 || _j == ly - 1; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            pf.iBoundaryConditionAlongXEdge(lx/5, 1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            pf.iBoundaryConditionAlongYEdge(ly/2, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            pf.iBoundaryConditionAlongXEdge(4*lx/5, -1, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            pf.iBoundaryConditionAlongYEdge(9*ly/10, 1, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            pf.SmoothCornerAt(lx/5, ly/2, -1, -1);
            pf.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
            pf.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
            pf.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
            AAD::iBoundaryConditionSetQAlongXEdge(pg, lx/5, 1, ux, uy, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AAD::iBoundaryConditionSetQAlongYEdge(pg, ly/2, 1, ux, uy, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            AAD::iBoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, ux, uy, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AAD::iBoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, ux, uy, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
            pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
            pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
            pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);

            std::swap(iux, iuxp);
            std::swap(iuy, iuyp);
            std::swap(iqx, iqxp);
            std::swap(iqy, iqyp);
        }

        //********************Get sensitivity********************
        double f_buffer = 0.0, f;
        for (int j = 0; j < pf.ny; ++j) {
            int i = lx/2 - pf.offsetx;
            if (0 <= i && i < pf.nx && (j + pf.offsety) > 9*ly/10) {
                int idx = pf.Index(i, j);
                f_buffer += ux[idx]*directionx[idx] + uy[idx]*directiony[idx];
            }
        }
        #ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        f /= (double)(ly/10);
        std::vector<double> dfdss(s.size(), 0.0);
        AAD::SensitivityBrinkmanDiffusivity(pg, dfdss.data(), ux, uy, imx, imy, dads, tem, item, iqx, iqy, gi, igi, diffusivity, dkds);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                int idx = pf.Index(i, j);
                dfdss[idx] = (j + pf.offsety) < ly/2 ? dfdss[idx] : 0.0;
            }
        }
        Normalize(dfdss.data(), pg.nxyz);

        //********************Filter sensitivities********************
        std::vector<double> dfds = DensityFilter::GetFilteredValue(pf, R, dfdss);
        std::vector<double> dgds = DensityFilter::GetFilteredValue(pf, R, dgdss);

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
        mnd *= 4.0/(double)(lx*ly/2)*100.0;

        //********************Update design variable********************
        if (MyRank == 0) {
            std::cout << "\rUpdate design variable" << std::string(10, ' ');
        }
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
        double dsmax_buffer = 0.0, dsmax;
        int imax_buffer = 0, imax, jmax_buffer = 0, jmax;
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                int idx = pf.Index(i, j);
                s[idx] = (j + pf.offsety) < ly/2 ? s[idx] : 1.0;
                double tmpds = fabs(s[idx] - snm1[idx]);
                if (dsmax_buffer < tmpds) {
                    dsmax_buffer = tmpds;
                    imax_buffer = i;
                    jmax_buffer = j;
                }
                snm1[idx] = s[idx];
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&dsmax_buffer, &dsmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&imax_buffer, &imax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&jmax_buffer, &jmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
        dsmax = dsmax_buffer;
        imax = imax_buffer;
        jmax = jmax_buffer;
#endif

        //********************Check convergence********************
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << k << " " << f << " " << g << " " << td << " " << ti << " " << dsmax  << " (" << imax << "," << jmax << ") " << qf << " " << qg << " " << mnd << std::endl;
        }
        if ((dsmax < 0.01 && g <= 0.0) || k == nk) {
            if (MyRank == 0) {
                std::cout << "----------Convergence/Last step----------" << std::endl;
            }

            VTKXMLExport file(pf, "result/ncpump");
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
            file.AddPointData(pf, "dfdss", [&](int _i, int _j, int _k) { return dfdss[pf.Index(_i, _j)];    });
            break;
        }
    }

    delete[] rho;   delete[] ux;    delete[] uy;    delete[] uxp;   delete[] uyp;
    delete[] tem;   delete[] qx;    delete[] qy;    delete[] qxp;   delete[] qyp;
    delete[] diffusivity;   delete[] alpha; delete[] dkds;  delete[] dads;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] iuxp;  delete[] iuyp;  delete[] imx;   delete[] imy;
    delete[] item;  delete[] iqx;   delete[] iqy;   delete[] iqxp;  delete[] iqyp;
    delete[] gi;    delete[] igi;   delete[] directionx;    delete[] directiony;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}