//#define _USE_MPI_DEFINES
#define _USE_AVX_DEFINES
#include <iostream>
#include <chrono>
#include <vector>
#include <utility>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d3q15.h"
#include "../src/equation/advection.h"
#include "../src/equation/adjointadvection.h"
#include "../src/utility/vtkxmlexport.h"
#include "../src/utility/residual.h"
#include "../src/utility/mma.h"
#include "../src/utility/densityfilter.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
#ifdef _USE_MPI_DEFINES
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int nPEx = atoi(argv[1]), nPEy = atoi(argv[2]), nPEz = atoi(argv[3]);
    assert(nPEx*nPEy*nPEz == PeTot);
#else
    int MyRank = 0, nPEx = 1, nPEy = 1, nPEz = 1;
#endif

    //********************Parameters********************
    int lx = 41, ly = 81, lz = 41, mx = 31, my = 61, mz = 31, nt = 30000, dt = 100, nitr = 2000, nb = 100;
    double Pr = 6.0, Ra = 2.5e3, nu = 0.1, L = 8.0, tem0 = 0.0, qn0 = 1.0e-2, alphamax = 1.0e4;
    double qf = 1e-2, qg = 1e0, movelimit = 0.2, weightlimit = 0.05, R = 1.5, eps = 1.0e-4, s0 = 0.05;

    double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1), gz = 0.0;
    D3Q15<double> pf(lx, ly, lz, MyRank, nPEx, nPEy, nPEz), pg(lx, ly, lz, MyRank, nPEx, nPEy, nPEz);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uz = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz], *uzp = new double[pf.nxyz];
    double *tem = new double[pg.nxyz], *qx = new double[pf.nxyz], *qy = new double[pf.nxyz], *qz = new double[pf.nxyz], *qxp = new double[pf.nxyz], *qyp = new double[pf.nxyz], *qzp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *iuz = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *imz = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz], *iuzp = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqz = new double[pg.nxyz], *iqxp = new double[pg.nxyz], *iqyp = new double[pg.nxyz], *iqzp = new double[pg.nxyz];
    for (int idx = 0; idx < pf.nxyz; idx++) {
        rho[idx] = 1.0;  ux[idx] = 0.0;  uy[idx] = 0.0;  uz[idx] = 0.0;  uxp[idx] = 0.0;  uyp[idx] = 0.0;  uzp[idx] = 0.0;  tem[idx] = 0.0; qx[idx] = 0.0;  qy[idx] = 0.0;  qz[idx] = 0.0;   qxp[idx] = 0.0; qyp[idx] = 0.0; qzp[idx] = 0.0;
        irho[idx] = 0.0; iux[idx] = 0.0; iuy[idx] = 0.0; iuz[idx] = 0.0; iuxp[idx] = 0.0; iuyp[idx] = 0.0; iuzp[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; imz[idx] = 0.0; item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0; iqz[idx] = 0.0; iqxp[idx] = 0.0; iqyp[idx] = 0.0; iqzp[idx] = 0.0;
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
            for (int k = 0; k < pf.nz; ++k) {
                if ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) {
                    int idx = pf.Index(i, j, k);
                    s[idx] = s0;
                    snm1[idx] = s0;
                }
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

    for (int itr = 1; itr <= nitr; itr++) {
        if (itr%nb == 0) {
            qf = std::min(1e7, qf*10.0);
        }

        //********************Filter variables********************
        std::vector<double> ss = DensityFilter::GetFilteredValue(pf, R, s);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                for (int k = 0; k < pf.nz; ++k) {
                    int idx = pf.Index(i, j, k);
                    ss[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) ? ss[idx] : 1.0;
                }
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
                for (int k = 0; k < pf.nz; ++k) {
                    if ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) {
                        int idx = pf.Index(i, j, k);
                        g_buffer += (1.0 - ss[idx])/(weightlimit*mx*my*mz);
                        dgdss[idx] = -1.0/(weightlimit*mx*my*mz); 
                    }
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
        NS::InitialCondition(pf, rho, ux, uy, uz);
        AD::InitialCondition(pg, tem, ux, uy, uz);
        for (int t = 1; t <= nt; t++) {
            AD::MacroBrinkmanCollideNaturalConvection(pf, rho, ux, uy, uz, alpha, nu, pg, tem, qx, qy, qz, diffusivity, gx, gy, gz, tem0, true, gi);
            if (t%dt == 0) {
                if (MyRank == 0) {
                    std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
                }
                if (Residual(ux, uy, uz, uxp, uyp, uzp, pf.nxyz) < eps && Residual(qx, qy, qz, qxp, qyp, qzp, pf.nxyz) < eps) {
                    td = t;
                    break;
                }
            }
            pf.Stream();
            pg.Stream();
            pf.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 1; });
            AD::BoundaryConditionSetT(pg, 
                [=](int _i, int _j, int _k) { return tem0; }, 
                ux, uy, uz,
                [=](int _i, int _j, int _k) { return _i == lx - 1 || _j == ly - 1 || _k == lz - 1; }
            );
            AD::BoundaryConditionSetQ(pg, 
                [=](int _i, int _j, int _k) { return (_j == 0 && _i < L && _k < L) ? qn0 : 0.0; }, 
                ux, uy, uz, diffusivity,
                [=](int _i, int _j, int _k) { return _j == 0; }
            );
            pg.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 0; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            std::swap(ux, uxp);
            std::swap(uy, uyp);
            std::swap(uz, uzp);
            std::swap(qx, qxp);
            std::swap(qy, qyp);
            std::swap(qz, qzp);
        }
               
        //  Inverse analyse
        int ti = nt;
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux, uy, uz, irho, iux, iuy, iuz);
        AAD::InitialCondition(pg, ux, uy, uz, item, iqx, iqy, iqz);
        for (int t = 1; t <= nt; t++) {
            AAD::MacroBrinkmanCollideNaturalConvection(
                pf, rho, ux, uy, uz, irho, iux, iuy, iuz, imx, imy, imz, alpha, nu, 
                pg, tem, item, iqx, iqy, iqz, diffusivity,
                gx, gy, gz, true, igi
            );
            if (t%dt == 0) {
                if (MyRank == 0) {
                    std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
                }
                if (Residual(iux, iuy, iuz, iuxp, iuyp, iuzp, pf.nxyz) < eps && Residual(iqx, iqy, iqz, iqxp, iqyp, iqzp, pf.nxyz) < eps) {
                    ti = t;
                    break;
                }
            }
            pf.iStream();
            pg.iStream();
            AAD::iBoundaryConditionSetT(pg, ux, uy, uz, [=](int _i, int _j, int _k) { return _i == lx - 1 || _j == ly - 1 || _k == lz - 1; });
            AAD::iBoundaryConditionSetQ(pg, ux, uy, uz, [=](int _i, int _j, int _k) { return _j == 0; });
            AAD::iBoundaryConditionSetQ(pg, ux, uy, uz, [=](int _i, int _j, int _k) { return _j == 0 && _i < L && _k < L; }, 1.0);
            pg.iBoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 0; });
            pf.iBoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 1; });
            pf.SmoothCorner();
            pg.SmoothCorner();

            std::swap(iux, iuxp);
            std::swap(iuy, iuyp);
            std::swap(iuz, iuzp);
            std::swap(iqx, iqxp);
            std::swap(iqy, iqyp);
            std::swap(iqz, iqzp);
        }
               
        //  Get sensitivity
        double f_buffer = 0.0, f;
        for (int k = 0; k < pf.nz; ++k) {
            for (int i = 0; i < pf.nx; ++i) {
                if ((i + pf.offsetx) < L && (k + pf.offsetz) < L && pf.PEy == 0) {
                    f_buffer += tem[pf.Index(i, 0, k)];
                }
            }
        }
        #ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        f /= L*L;
        std::vector<double> dfdss(s.size(), 0.0);
        AAD::SensitivityTemperatureAtHeatSource(
            ux, uy, uz, imx, imy, imz, pg, tem, item, iqx, iqy, iqz, gi, igi, dfdss.data(), diffusivity, dads, dkds,
            [=](int _i, int _j, int _k) { return (_j == 0 && _i < L && _k < L) ? qn0 : 0.0; }, 
            [=](int _i, int _j, int _k) { return _j == 0 && _i < L && _k < L; }
        );
        double dfdssmax_buffer = 0.0, dfdssmax;
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            if (dfdssmax_buffer < fabs(dfdss[idx])) {
                dfdssmax_buffer = fabs(dfdss[idx]);
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&dfdssmax_buffer, &dfdssmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
        dfdssmax = dfdssmax_buffer;
#endif
        for (int idx = 0; idx < pf.nxyz; ++idx) {
            dfdss[idx] /= dfdssmax;
        }
        
        //********************Filter sensitivities********************
        std::vector<double> dfds = DensityFilter::GetFilteredValue(pf, R, dfdss);
        std::vector<double> dgds = DensityFilter::GetFilteredValue(pf, R, dgdss);
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                for (int k = 0; k < pf.nz; ++k) {
                    int idx = pf.Index(i, j, k);
                    dfds[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) ? dfds[idx] : 0.0;
                    dgds[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) ? dgds[idx] : 0.0;
                }
            }
        }

        //********************Check grayscale********************
        double mnd_buffer = 0.0, mnd;
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                for (int k = 0; k < pf.nz; ++k) {
                    mnd_buffer += ss[pf.Index(i, j, k)]*(1.0 - ss[pf.Index(i, j, k)]);
                }
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&mnd_buffer, &mnd, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        mnd = mnd_buffer;
#endif
        mnd *= 4.0/(double)(mx*my*mz)*100.0;

        //********************Update design variable********************
        if (MyRank == 0) {
            std::cout << "\rUpdate design variable" << std::string(10, ' ');
        }
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
        double dsmax_buffer = 0.0, dsmax;
        int imax_buffer = 0, imax, jmax_buffer = 0, jmax, kmax_buffer = 0, kmax;
        for (int i = 0; i < pf.nx; ++i) {
            for (int j = 0; j < pf.ny; ++j) {
                for (int k = 0; k < pf.nz; ++k) {
                    int idx = pf.Index(i, j, k);
                    s[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz) ? s[idx] : 1.0;
                    double tmpds = fabs(s[idx] - snm1[idx]);
                    if (dsmax_buffer < tmpds) {
                        dsmax_buffer = tmpds;
                        imax_buffer = i;
                        jmax_buffer = j;
                        kmax_buffer = k;
                    }
                    snm1[idx] = s[idx];
                }
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&dsmax_buffer, &dsmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&imax_buffer, &imax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&jmax_buffer, &jmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&kmax_buffer, &kmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else
        dsmax = dsmax_buffer;
        imax = imax_buffer;
        jmax = jmax_buffer;
        kmax = kmax_buffer;
#endif

        //********************Check convergence********************
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << itr << " " << f << " " << g << " " << td << " " << ti << " " << dsmax  << " (" << imax << "," << jmax << "," << kmax << ") " << qf << " " << qg << " " << mnd << std::endl;
        }
        if ((itr > 1 && dsmax < 0.01 && g <= 0.0) || itr == nitr) {
            if (MyRank == 0) {
                std::cout << "----------Convergence/Last step----------" << std::endl;
            }

            VTKXMLExport file(pf, "result/nsadncopt3D");
            file.AddPointData(pf, "k", [&](int _i, int _j, int _k) { return diffusivity[pg.Index(_i, _j, _k)]; });
            file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j, _k)]; });
            file.AddPointData(pf, "u", 
                [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return uz[pf.Index(_i, _j, _k)]; }
            );
            file.AddPointData(pf, "T", [&](int _i, int _j, int _k) { return tem[pg.Index(_i, _j, _k)]; });
            file.AddPointData(pf, "q", 
                [&](int _i, int _j, int _k) { return qx[pg.Index(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return qy[pg.Index(_i, _j, _k)]; },
                [&](int _i, int _j, int _k) { return qz[pg.Index(_i, _j, _k)]; }
            );
            file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) {   return irho[pf.Index(_i, _j, _k)];  });
            file.AddPointData(pf, "iu", 
                [&](int _i, int _j, int _k) {   return iux[pf.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return iuy[pf.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return iuz[pf.Index(_i, _j, _k)];   }
            );
            file.AddPointData(pf, "im", 
                [&](int _i, int _j, int _k) {   return imx[pf.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return imy[pf.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return imz[pf.Index(_i, _j, _k)];   }
            );
            file.AddPointData(pf, "iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j, _k)];  });
            file.AddPointData(pf, "iq", 
                [&](int _i, int _j, int _k) {   return iqx[pg.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return iqy[pg.Index(_i, _j, _k)];   },
                [&](int _i, int _j, int _k) {   return iqz[pg.Index(_i, _j, _k)];   }
            );
            file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j, _k)];    });
            file.AddPointData(pf, "ss", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j, _k)];    });
            file.AddPointData(pf, "dfds", [&](int _i, int _j, int _k) { return dfds[pf.Index(_i, _j, _k)];    });
            file.AddPointData(pf, "dkds", [&](int _i, int _j, int _k) { return dkds[pf.Index(_i, _j, _k)];    });
            file.AddPointData(pf, "dads", [&](int _i, int _j, int _k) { return dads[pf.Index(_i, _j, _k)];    });
            
            for (int i = 0; i < pf.nx; ++i) {
                for (int j = 0; j < pf.ny; ++j) {
                    for (int k = 0; k < pf.nz; ++k) {
                        int idx = pf.Index(i, j, k);
                        ss[idx] = ((i + pf.offsetx) < mx && (j + pf.offsety) < my && (k + pf.offsetz) < mz && ss[idx] < 0.1) ? 0.0 : 1.0;
                    }
                }
            }
            file.AddPointData(pf, "ss2", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j, _k)];    });
            break;
        }
    }
    
    delete[] rho;   delete[] ux;    delete[] uy;    delete[] uz;    delete[] uxp;   delete[] uyp;   delete[] uzp;
    delete[] tem;   delete[] qx;    delete[] qy;    delete[] qz;    delete[] qxp;   delete[] qyp;   delete[] qzp;
    delete[] diffusivity;   delete[] alpha; delete[] dkds;  delete[] dads;
    delete[] irho;  delete[] iux;   delete[] iuy;   delete[] iuz;   delete[] iuxp;  delete[] iuyp;  delete[] iuzp;  delete[] imx;   delete[] imy;   delete[] imz;
    delete[] item;  delete[] iqx;   delete[] iqy;   delete[] iqz;   delete[] iqxp;  delete[] iqyp;  delete[] iqzp;
    delete[] gi;    delete[] igi;
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    if (MyRank == 0) {
        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
    }
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}