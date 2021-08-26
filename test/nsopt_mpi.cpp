#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "mpi.h"

#define _USE_MPI_DEFINES
#include "../src/particle/d2q9.h"
#include "../src/equation/navierstokes.h"
#include "../src/equation/adjointnavierstokes.h"
#include "../src/utility/mma.h"
#include "../src/utility/vtkxmlexport.h"

using namespace PANSLBM2;

int main(int argc, char** argv) {
    int PeTot, MyRank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &PeTot);
    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    assert(argc == 3);
    int mx = atoi(argv[1]), my = atoi(argv[2]);
    assert(mx*my == PeTot);

    //********************Set parameters********************
    int lx = 100, ly = 100, nt = 10000, dt = 100, nk = 1;
    double nu = 0.1, u0 = 0.01, rho0 = 1.0, q = 0.01, amax = 50.0, scale0 = 1.0e0, weightlimit = 0.25, movelimit = 0.5, epsu = 1.0e-4;
    D2Q9<double> pf(lx, ly, MyRank, mx, my);
    double *rho = new double[pf.nxy], *ux = new double[pf.nxy], *uy = new double[pf.nxy], *uxm1 = new double[pf.nxy], *uym1 = new double[pf.nxy];
    double *irho = new double[pf.nxy], *iux = new double[pf.nxy], *iuy = new double[pf.nxy], *imx = new double[pf.nxy], *imy = new double[pf.nxy];
    double *alpha = new double[pf.nxy], *dads = new double[pf.nxy];
    for (int idx = 0; idx < pf.nxy; ++idx) {
        rho[idx] = 1.0; ux[idx] = 0.0;  uy[idx] = 0.0;  uxm1[idx] = 0.0;    uym1[idx] = 0.0;
        irho[idx] = 0.0;    iux[idx] = 0.0; iuy[idx] = 0.0; imx[idx] = 0.0; imy[idx] = 0.0;
    }

    pf.SetBoundary([&](int _i, int _j) {
        if ((_i == 0 && 0.7*ly < _j && _j < 0.9*ly) || (_j == 0 && 0.7*lx < _i && _i < 0.9*lx)) {
            return 0;
        } else {
            return 1;
        }
    });
    int *boundaryup = new int[pf.nbc];
    pf.SetBoundary(boundaryup, [&](int _i, int _j) {
        if (_i == 0 && 0.7*ly < _j && _j < 0.9*ly) {
            return 1;
        } else if (_j == 0 && 0.7*lx < _i && _i < 0.9*lx) {
            return 2;
        } else {
            return 0;
        }
    });
    double *uxbc = new double[pf.nbc], *uybc = new double[pf.nbc], *rhobc = new double[pf.nbc], *usbc = new double[pf.nbc];
    for (int j = 0; j < pf.ny; ++j) {
        if (0.7*ly < j + pf.offsety && j + pf.offsety < 0.9*ly) {
            uxbc[j + pf.offsetxmin] = -u0*((j + pf.offsety) - 0.7*pf.ly)*((j + pf.offsety) - 0.9*pf.ly)/(0.1*pf.ly*0.1*pf.ly);
        } else {
            uxbc[j + pf.offsetxmin] = 0.0;
        }
        uybc[j + pf.offsetxmin] = 0.0;  rhobc[j + pf.offsetxmin] = 1.0; usbc[j + pf.offsetxmin] = 0.0;
        uxbc[j + pf.offsetxmax] = 0.0;  uybc[j + pf.offsetxmax] = 0.0;  rhobc[j + pf.offsetxmax] = 1.0; usbc[j + pf.offsetxmax] = 0.0;
    }
    for (int i = 0; i < pf.nx; ++i) {
        uxbc[i + pf.offsetymin] = 0.0;  uybc[i + pf.offsetymin] = 0.0;  rhobc[i + pf.offsetymin] = 1.0; usbc[i + pf.offsetymin] = 0.0;
        uxbc[i + pf.offsetymax] = 0.0;  uybc[i + pf.offsetymax] = 0.0;  rhobc[i + pf.offsetymax] = 1.0; usbc[i + pf.offsetymax] = 0.0;
    }

    std::vector<double> s(pf.nxy, 1.0);
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
        for(int idx = 0; idx < pf.nxy; idx++){
            g_buffer += s[idx]/(weightlimit*(double)(lx*ly));
            dgds[idx] = 1.0/(weightlimit*(double)(lx*ly)); 
        }
        MPI_Allreduce(&g_buffer, &g, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        g -= 1.0;

        //********************Direct Analyse********************
        if (MyRank == 0) {
            //std::cout << "Direct analyse t = 0";
        }
        NS::InitialCondition(pf, rho, ux, uy);
        for (int t = 1; t <= nt; ++t) {
            if (MyRank == 0 && t%dt == 0) {
                std::cout << "Direct analyse t = " << t << std::endl;
            }
            NS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, nu, alpha, true);
            pf.Swap();
            pf.Synchronize();
            pf.BoundaryCondition();
            NS::BoundaryConditionSetU(pf, uxbc, uybc, boundaryup);
            NS::BoundaryConditionSetRho(pf, rhobc, usbc, boundaryup);
            pf.SmoothCorner();
        }

        //********************Invert analyze********************
        if (MyRank == 0) {
            //std::cout << "\rInverse analyse t = 0" << std::string(10, ' ');
        }
        ANS::InitialCondition(pf, ux, uy, irho, iux, iuy);
        for (int t = 1; t <= nt; ++t) {
            if (MyRank == 0 && t%dt == 0) {
                std::cout << "Inverse analyse t = " << t << std::endl;
            }
            ANS::Macro_Brinkman_Collide_Stream(pf, rho, ux, uy, irho, iux, iuy, imx, imy, nu, alpha, true);
            pf.Swap();
            pf.iSynchronize();
            pf.iBoundaryCondition();
            ANS::BoundaryConditionSetiU(pf, uxbc, uybc, boundaryup);
            ANS::BoundaryConditionSetiRho<double>(pf, boundaryup);
            pf.SmoothCorner();
        }

        //********************Get sensitivity********************
        double f_buffer = 0.0, f;
        for (int j = 0; j < pf.ny; ++j) {
            if (0.7*ly < j + pf.offsety && j + pf.offsety < 0.9*ly) {
                f_buffer += rho[pf.Index(0, j)]/3.0;
            }
        }
        for (int i = 0; i < pf.nx; i++) {
            if (0.7*lx < i + pf.offsetx && i + pf.offsetx < 0.9*lx) {
                f_buffer -= rho[pf.Index(i, 0)]/3.0;
            }
        }
        MPI_Allreduce(&f_buffer, &f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        double dfdsmax_buffer = 0.0, dfdsmax;
        std::vector<double> dfds(s.size(), 0.0);
        for (int idx = 0; idx < pf.nxy; ++idx) {
            dfds[idx] = 3.0*(imx[idx]*ux[idx] + imy[idx]*uy[idx])*dads[idx];
            if (dfdsmax_buffer < fabs(dfds[idx])) {
                dfdsmax_buffer = fabs(dfds[idx]);
            }
        }
        MPI_Allreduce(&dfdsmax_buffer, &dfdsmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        for (int idx = 0; idx < pf.nxy; ++idx) {  
            dfds[idx] /= dfdsmax;
        }
        if (MyRank == 0) {
            std::cout << "\r" << std::fixed << std::setprecision(6) << k << "\t" << f << "\t" << g << std::endl;
        }

        //********************Update variable********************
        /*if (MyRank == 0) {
            std::cout << "Update design variable" << std::endl;
        }
        optimizer.UpdateVariables(s, f, dfds, { g }, { dgds });

        //********************Check convergence********************
        if(g < 0.0 && optimizer.IsConvergence(f)){
            if (MyRank == 0) {
                std::cout << std::endl << "-----Optimized-----" << std::endl;
            }
            break;
        }*/
    }

    //********************Export result********************
    VTKXMLExport file("result/pipebend", MyRank, lx, ly, 1, mx, my, 1);
    file.AddPointScaler("p", [&](int _i, int _j, int _k) { return rho[pf.Index(_i, _j)]/3.0; });
    file.AddPointVector("u", 
        [&](int _i, int _j, int _k) { return ux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return uy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    //file.AddPointScaler("dfds", [&](int _i, int _j, int _k) {   return dfds[pf.Index(_i, _j)];  });
    file.AddPointScaler("ip", [&](int _i, int _j, int _k) { return irho[pf.Index(_i, _j)]; });
    file.AddPointVector("iu", 
        [&](int _i, int _j, int _k) { return iux[pf.Index(_i, _j)]; },
        [&](int _i, int _j, int _k) { return iuy[pf.Index(_i, _j)]; },
        [](int _i, int _j, int _k) { return 0.0; }
    );
    file.AddPointScaler("alpha", [&](int _i, int _j, int _k) { return alpha[pf.Index(_i, _j)]; });
    file.AddPointScaler("s", [&](int _i, int _j, int _k) { return s[pf.Index(_i, _j)]; });

    delete[] rho, ux, uy, uxm1, uym1, irho, iux, iuy, imx, imy, alpha, dads, boundaryup, uxbc, uybc, rhobc, usbc;

    MPI_Finalize();

    return 0;
}