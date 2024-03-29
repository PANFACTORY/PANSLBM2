//#define _USE_MPI_DEFINES
//#define _USE_AVX_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/utility/mma.h"
#include "../src/utility/residual.h"
#include "../src/utility/vtkxmlexport.h"
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

    //********************Set parameters********************
    int lx = 201, ly = 201, nt = 10000, dt = 100, nk = 100;
    double nu = 0.1, u0 = 0.01, rho0 = 1.0, q = 0.01, amax = 2e2, scale0 = 1.0e0, weightlimit = 0.25, movelimit = 0.5, epsu = 1.0e-4;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    double *rho = new double[pf.nxyz], *ux = new double[pf.nxyz], *uy = new double[pf.nxyz], *uxp = new double[pf.nxyz], *uyp = new double[pf.nxyz];
    double *irho = new double[pf.nxyz], *iux = new double[pf.nxyz], *iuy = new double[pf.nxyz], *imx = new double[pf.nxyz], *imy = new double[pf.nxyz], *iuxp = new double[pf.nxyz], *iuyp = new double[pf.nxyz];
    double *alpha = new double[pf.nxyz], *dads = new double[pf.nxyz];
    for (int idx = 0; idx < pf.nxyz; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;  uxp[idx] = 0.0; uyp[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0; iuxp[idx] = 0.0;    iuyp[idx] = 0.0;
    }

    std::vector<double> s(pf.nxyz, 1.0);
    MMA<double> optimizer(s.size(), 1, 1.0,
		std::vector<double>(1, 0.0),
		std::vector<double>(1, 10000.0),
		std::vector<double>(1, 0.0), 
		std::vector<double>(s.size(), 0.0), std::vector<double>(s.size(), 1.0));
    optimizer.move = movelimit;
	
    for (int k = 0; k < nk; k++) {
        //********************Set parameter********************
        for (int idx = 0; idx < s.size(); idx++) {
            alpha[idx] = amax/(double)lx*q*(1.0 - s[idx])/(s[idx] + q);
            dads[idx] = -amax/(double)lx*q*(q + 1.0)/pow(q + s[idx], 2.0);
        }

        //********************Get weight********************
        double g_buffer = 0.0, g;
        std::vector<double> dgds(s.size(), 0.0);
        for(int idx = 0; idx < pf.nxyz; idx++){
            g_buffer += s[idx]/(weightlimit*(double)(lx*ly));
            dgds[idx] = 1.0/(weightlimit*(double)(lx*ly)); 
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&g_buffer, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        g = g_buffer;
#endif
        g -= 1.0;

        //********************Direct Analyse********************
        if (MyRank == 0) {
            std::cout << "Direct analyse t = 0";
        }
        NS::InitialCondition(pf, rho, ux, uy);
        int td;
        for (td = 1; td <= nt; ++td) {
            NS::MacroBrinkmanCollide(pf, rho, ux, uy, nu, alpha, true);
            if (td%dt == 0 && MyRank == 0) {
                std::cout << "\rDirect analyse t = " << td << std::string(10, ' ');
                if (Residual(ux, uy, uxp, uyp, pf.nxyz) < epsu) {
                    break;
                }
            }
            pf.Stream();
            pf.BoundaryCondition([=](int _i, int _j) { return ((_i == 0 && 0.7*ly < _j && _j < 0.9*ly) || (_j == 0 && 0.7*lx < _i && _i < 0.9*lx)) ? 0 : 1; });
            NS::BoundaryConditionSetU(pf, 
                [=](int _i, int _j) { return -u0*(_j - 0.7*ly)*(_j - 0.9*ly)/(0.1*ly*0.1*ly); }, 
                [=](int _i, int _j) { return 0.0; }, 
                [=](int _i, int _j) { return _i == 0 && 0.7*ly < _j && _j < 0.9*ly; }
            );
            NS::BoundaryConditionSetRho(pf, 
                [=](int _i, int _j) { return 1.0; }, 
                [=](int _i, int _j) { return 0.0; }, 
                [=](int _i, int _j) { return _j == 0 && 0.7*lx < _i && _i < 0.9*lx; }
            );
            pf.SmoothCorner();

            std::swap(ux, uxp);
            std::swap(uy, uyp);
        }

        //********************Invert analyze********************
        if (MyRank == 0) {
            std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
        int ti;
        for (ti = 1; ti <= nt; ++ti) {
            ANS::MacroBrinkmanCollide(pf, rho, ux, uy, irho, iux, iuy, imx, imy, nu, alpha, true);
            if (ti%dt == 0 && MyRank == 0) {
                std::cout << "\rInverse analyse t = " << ti << std::string(10, ' ');
                if (Residual(iux, iuy, iuxp, iuyp, pf.nxyz) < epsu) {
                    break;
                }
            }
            pf.iStream();
            pf.iBoundaryCondition([=](int _i, int _j) { return ((_i == 0 && 0.7*ly < _j && _j < 0.9*ly) || (_j == 0 && 0.7*lx < _i && _i < 0.9*lx)) ? 0 : 1; });
            ANS::iBoundaryConditionSetU(pf, 
                [=](int _i, int _j) { return -u0*(_j - 0.7*ly)*(_j - 0.9*ly)/(0.1*ly*0.1*ly); }, 
                [=](int _i, int _j) { return 0.0; }, 
                [=](int _i, int _j) { return _i == 0 && 0.7*ly < _j && _j < 0.9*ly; },
                1.0
            );
            ANS::iBoundaryConditionSetRho2D(pf, [=](int _i, int _j) { return _j == 0 && 0.7*lx < _i && _i < 0.9*lx; });
            pf.SmoothCorner();

            std::swap(iux, iuxp);
            std::swap(iuy, iuyp);
        }

        //********************Get sensitivity********************
        double f_buffer = 0.0, f;
        for (int j = 0; j < pf.ny; ++j) {
            if (0.7*ly < j + pf.offsety && j + pf.offsety < 0.9*ly && pf.PEx == 0) {
                f_buffer += rho[pf.Index(0, j)]/3.0;
            }
        }
        for (int i = 0; i < pf.nx; i++) {
            if (0.7*lx < i + pf.offsetx && i + pf.offsetx < 0.9*lx && pf.PEy == 0) {
                f_buffer -= rho[pf.Index(i, 0)]/3.0;
            }
        }
#ifdef _USE_MPI_DEFINES
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
        f = f_buffer;
#endif
        std::vector<double> dfds(s.size(), 0.0);
        ANS::SensitivityBrinkman(pf, dfds.data(), ux, uy, imx, imy, dads);
        Normalize(dfds.data(), pf.nxyz);
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << k << " " << f << " " << g << " " << td << " " << ti << std::endl;
        }

        //********************Check convergence********************
        if(g < 0.0 && optimizer.IsConvergence(f)){
            if (MyRank == 0) {
                std::cout << std::endl << "-----Optimized-----" << std::endl;
            }
            break;
        }

        //********************Update variable********************
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });
    }

    //********************Export result********************
    VTKXMLExport file(pf, "result/pipebend");
    file.AddPointData(pf, "p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointData(pf, "u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    //file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j)];  });
    file.AddPointData(pf, "ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointData(pf, "alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointData(pf, "s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });

    delete[] rho; delete[] ux; delete[] uy; delete[] uxp; delete[] uyp; delete[] irho; delete[] iux; delete[] iuy; delete[] imx; delete[] imy; delete[] iuxp; delete[] iuyp; delete[] alpha; delete[] dads;
#ifdef _USE_MPI_DEFINES
    MPI_Finalize();
#endif
    return 0;
}