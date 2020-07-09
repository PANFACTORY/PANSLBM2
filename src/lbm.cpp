//*****************************************************************************
//  Title       :   src/lbm.cpp
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/07
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <chrono>


#include <xmmintrin.h>


int main() {
    //--------------------Parameters--------------------
    int nx = 200;
    int ny = 80;
    int tmax = 30000;

    double dx = 1.0;
    double dt = 1.0;
    double viscosity = 0.02;
    double ux0 = 0.1;

    double omega = 1.0/(3.0*viscosity*dt/pow(dx, 2.0) + 0.5);
    double t0 = 4.0/9.0, t1 = 1.0/9.0, t2 = 1.0/36.0;


    //--------------------Variables--------------------
    std::vector<std::vector<double> > f0t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t0*(1.0 - 1.5*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f1t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t1*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f2t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t1*(1.0 - 1.5*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f3t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t1*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f4t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t1*(1.0 - 1.5*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f5t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t2*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f6t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t2*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f7t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t2*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0))));
    std::vector<std::vector<double> > f8t = std::vector<std::vector<double> >(nx, std::vector<double>(ny, t2*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0))));

    std::vector<std::vector<double> > f0tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f1tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f2tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f3tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f4tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f5tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f6tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f7tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > f8tp1 = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));

    std::vector<std::vector<double> > rho = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > ux = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double> > uy = std::vector<std::vector<double> >(nx, std::vector<double>(ny, 0.0));

    std::vector<std::vector<bool> > barrier0 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier1 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier2 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier3 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier4 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier5 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier6 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier7 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));
    std::vector<std::vector<bool> > barrier8 = std::vector<std::vector<bool> >(nx, std::vector<bool>(ny, false));


    //--------------------Set barrier--------------------
    for (int j = ny/2 - 8; j <= ny/2 + 8; j++) {
        barrier0[ny/2][j] = true;
    }

    for (int i = 1; i < nx - 1; i++) {
        for (int j = 1; j < ny - 1; j++) {
            barrier1[i][j] = barrier0[i - 1][j];
            barrier2[i][j] = barrier0[i][j - 1];
            barrier3[i][j] = barrier0[i + 1][j];
            barrier4[i][j] = barrier0[i][j + 1];
            barrier5[i][j] = barrier0[i - 1][j - 1];
            barrier6[i][j] = barrier0[i + 1][j - 1];
            barrier7[i][j] = barrier0[i + 1][j + 1];
            barrier8[i][j] = barrier0[i - 1][j + 1];
        }
    }


    //--------------------Loop for time step--------------------
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now();

    for (int t = 0; t < tmax; t++) {
        //..........Stream..........
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                f0tp1[i][j] = f0t[i][j];
                f1tp1[i][j] = (i == 0) ? f1t[nx - 1][j] : f1t[i - 1][j];
                f2tp1[i][j] = (j == 0) ? f2t[i][ny - 1] :f2t[i][j - 1];
                f3tp1[i][j] = (i == nx - 1) ? f3t[0][j] : f3t[i + 1][j];
                f4tp1[i][j] = (j == ny - 1) ? f4t[i][0] : f4t[i][j + 1];
                f5tp1[i][j] = (i == 0 || j == 0) ? f5t[nx - 1][ny - 1] : f5t[i - 1][j - 1];
                f6tp1[i][j] = (i == nx - 1 || j == 0) ? f6t[0][ny - 1] : f6t[i + 1][j - 1];
                f7tp1[i][j] = (i == nx - 1 || j == ny - 1) ? f7t[0][0] : f7t[i + 1][j + 1];
                f8tp1[i][j] = (i == 0 || j == ny - 1) ? f8t[nx - 1][0] : f8t[i - 1][j + 1];
            }
        }


        //..........Boundary condition (Bouns-Back)..........
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                if (barrier1[i][j]) {
                    f1tp1[i][j] = f3t[i][j];
                }
                if (barrier2[i][j]) {
                    f2tp1[i][j] = f4t[i][j];
                }
                if (barrier3[i][j]) {
                    f3tp1[i][j] = f1t[i][j];
                }
                if (barrier4[i][j]) {
                    f4tp1[i][j] = f2t[i][j];
                }
                if (barrier5[i][j]) {
                    f5tp1[i][j] = f7t[i][j];
                }
                if (barrier6[i][j]) {
                    f6tp1[i][j] = f8t[i][j];
                }
                if (barrier7[i][j]) {
                    f7tp1[i][j] = f5t[i][j];
                }
                if (barrier8[i][j]) {
                    f8tp1[i][j] = f6t[i][j];
                }
            }
        }


        //..........Update macroscopic values..........
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                rho[i][j] = f0tp1[i][j] + f1tp1[i][j] + f2tp1[i][j] + f3tp1[i][j] + f4tp1[i][j] + f5tp1[i][j] + f6tp1[i][j] + f7tp1[i][j] + f8tp1[i][j];
                ux[i][j] = (f1tp1[i][j] - f3tp1[i][j] + f5tp1[i][j] - f6tp1[i][j] - f7tp1[i][j] + f8tp1[i][j])/rho[i][j];
                uy[i][j] = (f2tp1[i][j] - f4tp1[i][j] + f5tp1[i][j] + f6tp1[i][j] - f7tp1[i][j] - f8tp1[i][j])/rho[i][j];
            }
        } 


        //..........Export result..........
        /*if (t%100 == 0) {
            std::cout << t << std::endl;

            std::ofstream fout("result/result" + std::to_string(t/100) + ".vtk");
            fout << "# vtk DataFile Version 3.0" << std::endl;
            fout << "2D flow" << std::endl;
            fout << "ASCII" << std::endl;
            fout << "DATASET\tSTRUCTURED_GRID" << std::endl;
            fout << "DIMENSIONS\t" << nx << "\t" << ny << "\t" << 1 << std::endl;
            
            fout << "POINTS\t" << nx*ny << "\t" << "float" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    fout << dx*i << "\t" << dx*j << "\t" << 0.0 << std::endl;
                }
            }

            fout << "POINT_DATA\t" << nx*ny << std::endl;
            fout << "SCALARS\tcurl\tfloat" << std::endl;
            fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1) {
                        fout << 0.0 << std::endl;
                    } else {
                        fout << uy[i + 1][j] - uy[i - 1][j] - ux[i][j + 1] + ux[i][j - 1] << std::endl;
                    }
                }
            }
        }*/


        //..........Collision..........
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) {
                double ux2 = pow(ux[i][j], 2.0);
                double uy2 = pow(uy[i][j], 2.0);
                double u2 = ux2 + uy2;
                double uxuy = ux[i][j]*uy[i][j];
                double omu215 = 1.0 - 1.5*u2;

                f0t[i][j] = (1.0 - omega)*f0tp1[i][j] + omega*t0*rho[i][j]*omu215;
                f1t[i][j] = (1.0 - omega)*f1tp1[i][j] + omega*t1*rho[i][j]*(omu215 + 3.0*ux[i][j] + 4.5*ux2);
                f2t[i][j] = (1.0 - omega)*f2tp1[i][j] + omega*t1*rho[i][j]*(omu215 + 3.0*uy[i][j] + 4.5*uy2);
                f3t[i][j] = (1.0 - omega)*f3tp1[i][j] + omega*t1*rho[i][j]*(omu215 - 3.0*ux[i][j] + 4.5*ux2);
                f4t[i][j] = (1.0 - omega)*f4tp1[i][j] + omega*t1*rho[i][j]*(omu215 - 3.0*uy[i][j] + 4.5*uy2);
                f5t[i][j] = (1.0 - omega)*f5tp1[i][j] + omega*t2*rho[i][j]*(omu215 + 3.0*(ux[i][j] + uy[i][j]) + 4.5*(u2 + 2.0*uxuy));
                f6t[i][j] = (1.0 - omega)*f6tp1[i][j] + omega*t2*rho[i][j]*(omu215 - 3.0*(ux[i][j] - uy[i][j]) + 4.5*(u2 - 2.0*uxuy));
                f7t[i][j] = (1.0 - omega)*f7tp1[i][j] + omega*t2*rho[i][j]*(omu215 - 3.0*(ux[i][j] + uy[i][j]) + 4.5*(u2 + 2.0*uxuy));
                f8t[i][j] = (1.0 - omega)*f8tp1[i][j] + omega*t2*rho[i][j]*(omu215 + 3.0*(ux[i][j] - uy[i][j]) + 4.5*(u2 - 2.0*uxuy));
            }
        }


        //..........Boundary condition (inlet)..........
        for (int j = 0; j < ny; j++) {
            f1t[0][j] = t1*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0));
            f3t[0][j] = t1*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0));
            f5t[0][j] = t2*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0));
            f6t[0][j] = t2*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0));
            f7t[0][j] = t2*(1.0 - 3.0*ux0 + 3.0*pow(ux0, 2.0));
            f8t[0][j] = t2*(1.0 + 3.0*ux0 + 3.0*pow(ux0, 2.0));
        }
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
}