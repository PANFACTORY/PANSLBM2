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

    assert(argc == 4);
    int nPEx = atoi(argv[1]), nPEy = atoi(argv[2]), nPEz = atoi(argv[3]);
    assert(nPEx*nPEy*nPEz == PeTot);
#else
    int MyRank = 0, nPEx = 1, nPEy = 1, nPEz = 1;
#endif

    //********************Parameters********************
    int lx = 81, ly = 161, lz = 81, nt = 90000, dt = 100, nitr = 2000, nb = 100;
    double Pr = 6.0, Ra = 2.5e4, nu = 0.1, L = (lx - 1)/10, tem0 = 0.0, qn = 1.0e-2, alphamax = 1.0e4;
    double qf = 1e-2, qfmax = 1e2, qg = 1e0, movelimit = 0.2, weightlimit = 0.05, R = 1.5, eps = 1.0e-5, s0 = 0.5;
    
    int mx = 3*(lx - 1)/4 + 1, my = 3*(ly - 1)/4 + 1, mz = 3*(lz - 1)/4 + 1;
    double U = nu*sqrt(Ra/Pr)/(double)(ly - 1), diff_fluid = nu/Pr, diff_solid = diff_fluid*10.0, gx = 0.0, gy = U*U/(double)(ly - 1), gz = 0.0;
    D3Q15<double> pf(lx, ly, lz, MyRank, nPEx, nPEy, nPEz), pg(lx, ly, lz, MyRank, nPEx, nPEy, nPEz);
    double **rho = new double*[nt], **ux = new double*[nt], **uy = new double*[nt], **uz = new double*[nt];
    double **tem = new double*[nt], **qx = new double*[nt], **qy = new double*[nt], **qz = new double*[nt];
    double **gi = new double*[nt];
    for (int t = 0; t < nt; ++t) {
        rho[t] = new double[pf.nxyz]; ux[t] = new double[pf.nxyz]; uy[t] = new double[pf.nxyz]; uz[t] = new double[pf.nxyz];
        tem[t] = new double[pg.nxyz]; qx[t] = new double[pg.nxyz]; qy[t] = new double[pg.nxyz]; qz[t] = new double[pg.nxyz];
        gi[t] = new double[pg.nxyz*pg.nc];
    }
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *iuz = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *imz = new double[pf.nxyz];
    double *item = new double[pg.nxyz], *iqx = new double[pg.nxyz], *iqy = new double[pg.nxyz], *iqz = new double[pg.nxyz];
    double *alpha = new double[pf.nxyz], *diffusivity = new double[pf.nxyz], *dads = new double[pf.nxyz], *dkds = new double[pf.nxyz];
    double *igi = new double[pg.nxyz*pg.nc];
    
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

    for (int itr = 1, cnt = 1; itr <= nitr; itr++) {
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
        if (MyRank == 0) {
            std::cout << "Direct analyse t = 0";
        }
        for (int idx = 0; idx < pf.nxyz; idx++) {
            rho[0][idx] = 1.0; ux[0][idx] = 0.0; uy[0][idx] = 0.0; uz[0][idx] = 0.0; 
            tem[0][idx] = 0.0; qx[0][idx] = 0.0; qy[0][idx] = 0.0; qz[0][idx] = 0.0;
        }   
        NS::InitialCondition(pf, rho[0], ux[0], uy[0], uz[0]);
        AD::InitialCondition(pg, tem[0], ux[0], uy[0], uz[0]);
        for (int t = 1; t < nt; ++t) {
            AD::MacroBrinkmanCollideNaturalConvection(
                pf, rho[t], ux[t], uy[t], uz[t], alpha, nu, 
                pg, tem[t], qx[t], qy[t], qz[t], diffusivity, 
                gx, gy, gz, tem0, true, gi[t]
            );
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rDirect analyse t = " << t << std::string(10, ' ');
            }
            pf.Stream();
            pg.Stream();
            pf.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 1; });
            AD::BoundaryConditionSetT(pg, 
                [=](int _i, int _j, int _k) { return tem0; }, 
                ux[t], uy[t], uz[t],
                [=](int _i, int _j, int _k) { return _i == lx - 1 || _j == ly - 1 || _k == lz - 1; }
            );
            AD::BoundaryConditionSetQ(pg, 
                [=](int _i, int _j, int _k) { return (_j == 0 && _i < L && _k < L) ? qn : 0.0; }, 
                ux[t], uy[t], uz[t], diffusivity,
                [=](int _i, int _j, int _k) { return _j == 0; }
            );
            pg.BoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 0; });
            pf.SmoothCorner();
            pg.SmoothCorner();
        }
                
        //  Inverse analyse
        std::vector<double> dfdss(s.size(), 0.0);
        for (int idx = 0; idx < pf.nxyz; idx++) {
            irho[idx] = 0.0; iux[idx] = 0.0; iuy[idx] = 0.0; iuz[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; imz[idx] = 0.0;
            item[idx] = 0.0; iqx[idx] = 0.0; iqy[idx] = 0.0; iqz[idx] = 0.0;
        }
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux[nt - 1], uy[nt - 1], uz[nt - 1], irho, iux, iuy, iuz);
        AAD::InitialCondition(pg, ux[nt - 1], uy[nt - 1], uz[nt - 1], item, iqx, iqy, iqz);
        for (int t = nt - 2; t >= 0; --t) {
            AAD::MacroBrinkmanCollideNaturalConvection(
                pf, rho[t], ux[t], uy[t], uz[t], irho, iux, iuy, iuz, imx, imy, imz, alpha, nu, 
                pg, tem[t], item, iqx, iqy, iqz, diffusivity,
                gx, gy, gz, true, igi
            );
            AAD::SensitivityTemperatureAtHeatSource(
                pg, dfdss.data(), ux[t], uy[t], uz[t], imx, imy, imz, dads, tem[t], item, iqx, iqy, iqz, gi[t], igi, diffusivity, dkds,
                [=](int _i, int _j, int _k) { return (_j == 0 && _i < L && _k < L) ? qn : 0.0; }, 
                [=](int _i, int _j, int _k) { return _j == 0 && _i < L && _k < L; }
            );
            if (t%dt == 0 && MyRank == 0) {
                std::cout << "\rInverse analyse t = " << t << std::string(10, ' ');
            }
            pf.iStream();
            pg.iStream();
            AAD::iBoundaryConditionSetT(pg, ux[t], uy[t], uz[t], [=](int _i, int _j, int _k) { return _i == lx - 1 || _j == ly - 1 || _k == lz - 1; });
            AAD::iBoundaryConditionSetQ(pg, ux[t], uy[t], uz[t], [=](int _i, int _j, int _k) { return _j == 0; });
            AAD::iBoundaryConditionSetQ(pg, ux[t], uy[t], uz[t], [=](int _i, int _j, int _k) { return _j == 0 && _i < L && _k < L; }, 1.0);
            pg.iBoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 0; });
            pf.iBoundaryCondition([=](int _i, int _j, int _k) { return (_i == 0 || _k == 0) ? 2 : 1; });
            pf.SmoothCorner();
            pg.SmoothCorner();
        }
               
        //  Get sensitivity
        double f_buffer = 0.0, f;
        for (int t = 0; t < nt; ++t) {
            for (int i = 0; i < pf.nx; ++i) {
                for (int k = 0; k < pf.nz; ++k) {
                    if ((i + pf.offsetx) < L && (k + pf.offsetz) < L && pf.PEy == 0) {
                        f_buffer += tem[t][pf.Index(i, 0, k)];
                    }
                }
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        f /= (double)(L*L*nt);
        
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
        Normalize(dfds.data(), pf.nxyz);

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
                        imax_buffer = i + pf.offsetx;
                        jmax_buffer = j + pf.offsety;
                        kmax_buffer = k + pf.offsetz;
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
            std::cout << "\r" << std::fixed << std::setprecision(6) << itr << " " << f << " " << g << " " << dsmax  << " (" << imax << "," << jmax << "," << kmax << ") " << qf << " " << qg << " " << mnd << std::endl;
        }
        if (dsmax < 0.01 || itr == nitr) {
            if (qf < qfmax && itr != nitr) {
                cnt = 0;
            } else {
                if (MyRank == 0) {
                    std::cout << "----------Convergence/Last step----------" << std::endl;
                    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
                    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
                }
#ifdef _USE_MPI_DEFINES
                VTKXMLExport file(pf, "result/heatsink3D_transient");
                file.AddPointData(pf, "rho", [&](int _i, int _j, int _k) { return rho[nt - 1][pf.Index(_i, _j, _k)]; });
                file.AddPointData(pf, "u", 
                    [&](int _i, int _j, int _k) { return ux[nt - 1][pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return uy[nt - 1][pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return uz[nt - 1][pf.Index(_i, _j, _k)]; }
                );
                file.AddPointData(pf, "T", [&](int _i, int _j, int _k) { return tem[nt - 1][pg.Index(_i, _j, _k)]; });
                file.AddPointData(pf, "q", 
                    [&](int _i, int _j, int _k) { return qx[nt - 1][pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return qy[nt - 1][pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return qz[nt - 1][pg.Index(_i, _j, _k)]; }
                );
                file.AddPointData(pf, "irho", [&](int _i, int _j, int _k) {   return irho[pf.Index(_i, _j, _k)];  });
                file.AddPointData(pf, "iu", 
                    [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iuz[pf.Index(_i, _j, _k)]; }
                );
                file.AddPointData(pf, "im", 
                    [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return imz[pf.Index(_i, _j, _k)]; }
                );
                file.AddPointData(pf, "iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j, _k)]; });
                file.AddPointData(pf, "iq", 
                    [&](int _i, int _j, int _k) { return iqx[pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iqy[pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iqz[pg.Index(_i, _j, _k)]; }
                );
                file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j, _k)]; });
                file.AddPointData(pf, "ss", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j, _k)]; });
                file.AddPointData(pf, "dfdss", [&](int _i, int _j, int _k) { return dfdss[pf.Index(_i, _j, _k)]; });
#else
                VTKExport file("result/heatsink3D_transient.vtk", lx, ly, lz);
                file.AddPointScaler("rho", [&](int _i, int _j, int _k) { return rho[nt - 1][pf.Index(_i, _j, _k)]; });
                file.AddPointVector("u", 
                    [&](int _i, int _j, int _k) { return ux[nt - 1][pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return uy[nt - 1][pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return uz[nt - 1][pf.Index(_i, _j, _k)]; }
                );
                file.AddPointScaler("T", [&](int _i, int _j, int _k) { return tem[nt - 1][pg.Index(_i, _j, _k)]; });
                file.AddPointVector("q", 
                    [&](int _i, int _j, int _k) { return qx[nt - 1][pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return qy[nt - 1][pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return qz[nt - 1][pg.Index(_i, _j, _k)]; }
                );
                file.AddPointScaler("irho", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j, _k)]; });
                file.AddPointVector("iu", 
                    [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iuz[pf.Index(_i, _j, _k)]; }
                );
                file.AddPointVector("im", 
                    [&](int _i, int _j, int _k) { return imx[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return imy[pf.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return imz[pf.Index(_i, _j, _k)]; }
                );
                file.AddPointScaler("iT", [&](int _i, int _j, int _k) { return item[pg.Index(_i, _j, _k)]; });
                file.AddPointVector("iq", 
                    [&](int _i, int _j, int _k) { return iqx[pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iqy[pg.Index(_i, _j, _k)]; },
                    [&](int _i, int _j, int _k) { return iqz[pg.Index(_i, _j, _k)]; }
                );
                file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j, _k)]; });
                file.AddPointScaler("ss", [&](int _i, int _j, int _k) { return ss[pf.Index(_i, _j, _k)]; });
                file.AddPointScaler("dfdss", [&](int _i, int _j, int _k) { return dfdss[pf.Index(_i, _j, _k)]; });
#endif
                break;
            }
        }
    }
    
    for (int t = 0; t < nt; ++t) {
        delete[] rho[t]; delete[] ux[t]; delete[] uy[t]; delete[] uz[t];
        delete[] tem[t]; delete[] qx[t]; delete[] qy[t]; delete[] qz[t];
        delete[] gi[t];    
    }
    delete[] rho; delete[] ux; delete[] uy; delete[] uz; 
    delete[] tem; delete[] qx; delete[] qy; delete[] qz; delete[] gi;
    delete[] diffusivity; delete[] alpha; delete[] dkds; delete[] dads;
    delete[] irho; delete[] iux; delete[] iuy; delete[] iuz; delete[] imx; delete[] imy; delete[] imz; 
    delete[] item; delete[] iqx; delete[] iqy; delete[] iqz; delete[] igi;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}