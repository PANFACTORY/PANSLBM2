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
#include "../src/utility/residual.h"
#include "../src/utility/mma.h"
#include "../src/utility/densityfilter.h"
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
    int mx = atoi(argv[1]), my = atoi(argv[2]);
    assert(mx*my == PeTot);
#else
    int MyRank = 0, mx = 1, my = 1; 
#endif

    //--------------------Set parameters--------------------
    int lx = 101, ly = 201, dt = 1000, nt0 = 1000000, nt = 100000, period = 100000, nk = 2000, nb = 100, duty = 20;
    double viscosity = 0.1/6.0, diff_fluid = viscosity/1.0, Th = 1.0, Tl = 0.0, gx = 0.0, gy = 1000*pow(viscosity, 2)/(double)pow(lx - 1, 3);
    double alphamax = 1e5, diff_solid = diff_fluid*10.0, qf = 1e-6, qfmax = 1e-2, qg = 1e-4, weightlimit = 0.5, movelimit = 0.2, R = 0.5, eps = 1e-5, ratio = 0.5;
    D2Q9<double> pf(lx, ly, MyRank, mx, my), pg(lx, ly, MyRank, mx, my);
    double **rho = new double*[nt], **ux = new double*[nt], **uy = new double*[nt];
    double **tem = new double*[nt], *qx = new double[pg.nxyz], *qy = new double[pg.nxyz];
    double **gi = new double*[nt];
    double *f = new double[nt];
    f[0] = 0.0;
    for (int t = 0; t < nt; ++t) {
        rho[t] = new double[pf.nxyz];   ux[t] = new double[pf.nxyz];    uy[t] = new double[pf.nxyz];
        tem[t] = new double[pg.nxyz];   gi[t] = new double[pg.nxyz*pg.nc];
    }
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[0][idx] = 1.0; ux[0][idx] = 0.0; uy[0][idx] = 0.0; tem[0][idx] = 0.5*(Tl + Th); qx[idx] = 0.0; qy[idx] = 0.0;
        irho[idx] = 1.0; iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0;
    }
    double *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dads = new double[pf.nxyz], *dkds = new double[pf.nxyz];
    double *igi = new double[pg.nxyz*pg.nc];
    double *directionx = new double[pf.nxyz], *directiony = new double[pf.nxyz], *directionxt = new double[pf.nxyz], *directionyt = new double[pf.nxyz];
    for (int i = 0; i < pf.nx; ++i) {
        for (int j = 0; j < pf.ny; ++j) {
            int idx = pf.Index(i, j);
            directionx[idx] = ((i + pf.offsetx) == lx/2 && (j + pf.offsety) > 9*ly/10) ? -1.0 : 0.0;
            directiony[idx] = 0.0;
        }
    }

    auto tembc = [=](int _t) { return Th*(1 - cos(2*M_PI*_t/period)); };
    //auto tembc = [=](int _t) { return _t%period < period*duty/100.0 ? (Th - Tl)*100.0/(double)duty + Tl : Tl; };

    std::vector<double> s(pf.nxyz, 1.0), snm1(pf.nxyz, 1.0);
    MMA<double> optimizer(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), 
        std::vector<double>(s.size(), 1.0)
    );
    optimizer.move = movelimit;

    for (int k = 1, cnt = 1; k <= nk; k++) {
        if (cnt%nb == 0) {
            qf = std::min(qfmax, qf*10.0);
            cnt = 1;
        } else {
            cnt++;
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
        if (MyRank == 0) {
            std::cout << "Direct analyse t = 0";
        }
        if (k == 1) {
            NS::InitialCondition(pf, rho[0], ux[0], uy[0]);
            AD::InitialCondition(pg, tem[0], ux[0], uy[0]);
            for (int t = 1; t < nt0; ++t) {
                if (t%dt == 0 && MyRank == 0) {
                    std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                }
                AD::MacroBrinkmanCollideNaturalConvection(
                    pf, rho[0], ux[0], uy[0], alpha, viscosity, 
                    pg, tem[0], qx, qy, diffusivity, gx, gy, 0.5*(Th + Tl), true, gi[0]
                );

                pf.Stream();
                pg.Stream();
                pf.BoundaryCondition([=](int _i, int _j) { return 1; });
                pg.BoundaryCondition([=](int _i, int _j) { return 0; });
                AD::BoundaryConditionSetT(pg, 
                    [=](int _i, int _j) { return _i == 0 ? tembc(t) : Tl; }, 
                    ux[0], uy[0], 
                    [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); }
                );
                AD::BoundaryConditionSetQ(pg, 
                    [=](int _i, int _j) { return 0.0; },
                    ux[0], uy[0], diffusivity, 
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
                AD::BoundaryConditionSetQAlongXEdge(pg, lx/5, 1, [=](int _i, int _j) { return 0.0; }, ux[0], uy[0], diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                AD::BoundaryConditionSetQAlongYEdge(pg, ly/2, 1, [=](int _i, int _j) { return 0.0; }, ux[0], uy[0], diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                AD::BoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, [=](int _i, int _j) { return 0.0; }, ux[0], uy[0], diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
                AD::BoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, [=](int _i, int _j) { return 0.0; }, ux[0], uy[0], diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
                pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
                pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
                pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
                pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
            }
        } else {
            for (int idx = 0; idx < pf.nxyz; idx++) {
                rho[0][idx] = rho[nt - 1][idx]; ux[0][idx] = ux[nt - 1][idx]; uy[0][idx] = uy[nt - 1][idx]; tem[0][idx] = tem[nt - 1][idx];
            }
        }
        double faverage_buffer = 0.0, fsquare_buffer = 0.0;
        NS::InitialCondition(pf, rho[0], ux[0], uy[0]);
        AD::InitialCondition(pg, tem[0], ux[0], uy[0]);
        for (int t = 1; t < nt; ++t) {
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
            }
            AD::MacroBrinkmanCollideNaturalConvection(
                pf, rho[t], ux[t], uy[t], alpha, viscosity, 
                pg, tem[t], qx, qy, diffusivity, gx, gy, 0.5*(Th + Tl), true, gi[t]
            );

            pf.Stream();
            pg.Stream();
            pf.BoundaryCondition([=](int _i, int _j) { return 1; });
            pg.BoundaryCondition([=](int _i, int _j) { return 0; });
            AD::BoundaryConditionSetT(pg, 
                [=](int _i, int _j) { return _i == 0 ? tembc(t) : Tl; }, 
                ux[t], uy[t], 
                [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); }
            );
            AD::BoundaryConditionSetQ(pg, 
                [=](int _i, int _j) { return 0.0; },
                ux[t], uy[t], diffusivity, 
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
            AD::BoundaryConditionSetQAlongXEdge(pg, lx/5, 1, [=](int _i, int _j) { return 0.0; }, ux[t], uy[t], diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AD::BoundaryConditionSetQAlongYEdge(pg, ly/2, 1, [=](int _i, int _j) { return 0.0; }, ux[t], uy[t], diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            AD::BoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, [=](int _i, int _j) { return 0.0; }, ux[t], uy[t], diffusivity, [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AD::BoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, [=](int _i, int _j) { return 0.0; }, ux[t], uy[t], diffusivity, [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
            pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
            pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
            pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);

            f[t] = 0.0;
            for (int j = 0; j < pf.ny; ++j) {
                int i = lx/2 - pf.offsetx;
                if (0 <= i && i < pf.nx && (j + pf.offsety) > 9*ly/10) {
                    int idx = pf.Index(i, j);
                    f[t] += ux[t][idx]*directionx[idx] + uy[t][idx]*directiony[idx];
                }
            }
            f[t] /= (double)((ly - 1)/10);
            faverage_buffer += f[t];
            fsquare_buffer += pow(f[t], 2.0);
        }
        double faverage, fsquare;    
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&faverage_buffer, &faverage, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&fsquare_buffer, &fsquare, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        faverage = faverage_buffer;
        fsquare = fsquare_buffer;
#endif
        faverage /= (double)nt;
        double variance = fsquare/(double)nt - pow(faverage, 2.0), coef = (1.0 - ratio)/sqrt(variance);
        double F = ratio*faverage + (1.0 - ratio)*sqrt(variance);

        //********************Inverse analyze********************
        std::vector<double> dfdss(s.size(), 0.0);
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux[nt - 1], uy[nt - 1], irho, iux, iuy);
        AAD::InitialCondition(pg, ux[nt - 1], uy[nt - 1], item, iqx, iqy);
        for (int t = nt - 2; t >= 0; --t) {
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
            }
            for (int idx = 0; idx < pf.nxyz; ++idx) {
                directionxt[idx] = (ratio + coef*(f[t] - faverage))*directionx[idx];
                directionyt[idx] = (ratio + coef*(f[t] - faverage))*directiony[idx];
            }

            AAD::MacroBrinkmanCollideNaturalConvectionMassFlow(
                pf, rho[t], ux[t], uy[t], irho, iux, iuy, imx, imy, alpha, viscosity,
                pg, tem[t], item, iqx, iqy, diffusivity, gx, gy, 
                directionxt, directionyt, true, igi
            );

            AAD::SensitivityBrinkmanDiffusivity(pg, dfdss.data(), ux[t], uy[t], imx, imy, dads, tem[t], item, iqx, iqy, gi[t], igi, diffusivity, dkds);

            pf.iStream();
            pg.iStream();
            pf.iBoundaryCondition([=](int _i, int _j) { return 1; });
            pg.iBoundaryCondition([=](int _i, int _j) { return 0; });
            AAD::iBoundaryConditionSetT(pg, ux[t], uy[t], [=](int _i, int _j) { return (_i == 0 && _j < ly/2) || (_i == lx - 1 && _j < ly/2); });
            AAD::iBoundaryConditionSetQ(pg, ux[t], uy[t], [=](int _i, int _j) { return (_i == 0 && ly/2 <= _j) || (_i == lx - 1 && ly/2 <= _j) || _j == 0 || _j == ly - 1; });
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
            AAD::iBoundaryConditionSetQAlongXEdge(pg, lx/5, 1, ux[t], uy[t], [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AAD::iBoundaryConditionSetQAlongYEdge(pg, ly/2, 1, ux[t], uy[t], [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            AAD::iBoundaryConditionSetQAlongXEdge(pg, 4*lx/5, -1, ux[t], uy[t], [=](int _i, int _j) { return ly/2 <= _j && _j < 9*ly/10; });
            AAD::iBoundaryConditionSetQAlongYEdge(pg, 9*ly/10, 1, ux[t], uy[t], [=](int _i, int _j) { return lx/5 <= _i && _i < 4*lx/5; });
            pg.SmoothCornerAt(lx/5, ly/2, -1, -1);
            pg.SmoothCornerAt(4*lx/5, ly/2, 1, -1);
            pg.SmoothCornerAt(4*lx/5, 9*ly/10, 1, 1);
            pg.SmoothCornerAt(lx/5, 9*ly/10, -1, 1);
        }
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
        optimizer.UpdateVariables(s, F, dfds, { g }, { dgds });
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
            std::cout << "\r" << k << std::scientific << std::setprecision(6) << " " << F << " " << faverage << " " << variance << std::fixed << std::setprecision(6) << " " << g << " " << dsmax  << " (" << imax << "," << jmax << ") " << qf << " " << qg << " " << mnd << std::endl;
        }
        if (dsmax < 0.01 || k == nk) {
            if (qf < qfmax && k != nk) {
                cnt = 0;
            } else {
                if (MyRank == 0) {
                    std::cout << "----------Convergence/Last step----------" << std::endl;
                }
#ifdef _USE_MPI_DEFINES
                VTKXMLExport file(pf, "result/ncpump_periodic");
                file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[nt - 1][pf.Index(_i, _j)]; });
                file.AddPointData(pf, "u", 
                    [&](int _i, int _j, int _k) { return ux[nt - 1][pf.Index(_i, _j)]; },
                    [&](int _i, int _j, int _k) { return uy[nt - 1][pf.Index(_i, _j)]; },
                    [](int _i, int _j, int _k) { return 0.0; }
                );
                file.AddPointData(pf, "T", [&](int _i, int _j, int _k) { return tem[nt - 1][pg.Index(_i, _j)]; });
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
#else
                VTKExport file("result/ncpump_periodic.vtk", lx, ly);
                file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[nt - 1][pf.Index(_i, _j)]; });
                file.AddPointVector("u", 
                    [&](int _i, int _j, int _k) { return ux[nt - 1][pf.Index(_i, _j)]; },
                    [&](int _i, int _j, int _k) { return uy[nt - 1][pf.Index(_i, _j)]; },
                    [](int _i, int _j, int _k) { return 0.0; }
                );
                file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[nt - 1][pg.Index(_i, _j)]; });
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
                file.AddPointScaler("ss", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j)];    });
                file.AddPointScaler("dfdss", [&](int _i, int _j, int _k) { return dfdss[pf.Index(_i, _j)];    });
#endif
                break;
            }
        }
    }

    for (int t = 0; t < nt; ++t) {
        delete[] rho[t]; delete[] ux[t]; delete[] uy[t]; delete[] tem[t]; delete[] gi[t];
    }
    delete[] rho; delete[] ux; delete[] uy; delete[] tem; delete[] qx; delete[] qy; delete[] gi;
    delete[] diffusivity; delete[] alpha; delete[] dkds; delete[] dads;
    delete[] irho; delete[] iux; delete[] iuy; delete[] imx; delete[] imy; delete[] item; delete[] iqx; delete[] iqy;
    delete[] igi; delete[] directionx; delete[] directiony; delete[] f; delete[] directionxt; delete[] directionyt;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
}