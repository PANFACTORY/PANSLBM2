//*****************************************************************************
//  Title       :   src/lbm.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/07/10
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************


namespace PANSLBM2 {
    template<class T>
    class LBM {
public:
        LBM(int _nx, int _ny, T _viscosity);
        ~LBM();

        void Boundary();
        void Stream();
        void UpdateMacro();
        void Collision();

        const int nx, ny;

        std::vector<std::vector<T> > rho;
        std::vector<std::vector<T> > u;
        std::vector<std::vector<T> > v;

        const static int INLET = 0;
        const static int OUTLET = 1;
        const static int NOSRIP = 2;
        const static int SLIP = 3;


private:
        T dx, dt, viscosity, omega, t0, t1, t2;

        std::vector<std::vector<T> > f0t;
        std::vector<std::vector<T> > f1t;
        std::vector<std::vector<T> > f2t;
        std::vector<std::vector<T> > f3t;
        std::vector<std::vector<T> > f4t;
        std::vector<std::vector<T> > f5t;
        std::vector<std::vector<T> > f6t;
        std::vector<std::vector<T> > f7t;
        std::vector<std::vector<T> > f8t;

        std::vector<std::vector<T> > f0tp1;
        std::vector<std::vector<T> > f1tp1;
        std::vector<std::vector<T> > f2tp1;
        std::vector<std::vector<T> > f3tp1;
        std::vector<std::vector<T> > f4tp1;
        std::vector<std::vector<T> > f5tp1;
        std::vector<std::vector<T> > f6tp1;
        std::vector<std::vector<T> > f7tp1;
        std::vector<std::vector<T> > f8tp1;
    };


    template<class T>
    LBM<T>::LBM(int _nx, int _ny, T _viscosity) : nx(_nx), ny(_ny) {
        this->dx = 1.0;
        this->dt = 1.0;
        this->viscosity = _viscosity;
        this->omega = 1.0/(3.0*this->viscosity*this->dt/pow(this->dx, 2.0) + 0.5);
        this->t0 = 4.0/9.0;
        this->t1 = 1.0/9.0;
        this->t2 = 1.0/36.0;

        this->f0t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t0));
        this->f1t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f2t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f3t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f4t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t1));
        this->f5t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f6t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f7t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));
        this->f8t = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, t2));

        this->f0tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f1tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f2tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f3tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f4tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f5tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f6tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f7tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->f8tp1 = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));

        this->rho = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->u = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
        this->v = std::vector<std::vector<T> >(this->nx, std::vector<T>(this->ny, T()));
    }


    template<class T>
    LBM<T>::~LBM() {}


    template<class T>
    void LBM<T>::Stream() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->f0tp1[i][j] = this->f0t[i][j];
                this->f1tp1[i][j] = (i == 0) ? this->f1t[i][j] : this->f1t[i - 1][j];
                this->f2tp1[i][j] = (j == 0) ? this->f2t[i][j] : this->f2t[i][j - 1];
                this->f3tp1[i][j] = (i == this->nx - 1) ? this->f3t[i][j] : this->f3t[i + 1][j];
                this->f4tp1[i][j] = (j == this->ny - 1) ? this->f4t[i][j] : this->f4t[i][j + 1];
                this->f5tp1[i][j] = (i == 0 || j == 0) ? this->f5t[i][j] : this->f5t[i - 1][j - 1];
                this->f6tp1[i][j] = (i == this->nx - 1 || j == 0) ? this->f6t[i][j] : this->f6t[i + 1][j - 1];
                this->f7tp1[i][j] = (i == this->nx - 1 || j == this->ny - 1) ? this->f7t[i][j] : this->f7t[i + 1][j + 1];
                this->f8tp1[i][j] = (i == 0 || j == this->ny - 1) ? this->f8t[i][j] : this->f8t[i - 1][j + 1];
            }
        }
    }


    template<class T>
    void LBM<T>::UpdateMacro() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                this->rho[i][j] = this->f0tp1[i][j] + this->f1tp1[i][j] + this->f2tp1[i][j] + this->f3tp1[i][j] + this->f4tp1[i][j] + this->f5tp1[i][j] + this->f6tp1[i][j] + this->f7tp1[i][j] + this->f8tp1[i][j];
                this->u[i][j] = (this->f1tp1[i][j] - this->f3tp1[i][j] + this->f5tp1[i][j] - this->f6tp1[i][j] - this->f7tp1[i][j] + this->f8tp1[i][j])/this->rho[i][j];
                this->v[i][j] = (this->f2tp1[i][j] - this->f4tp1[i][j] + this->f5tp1[i][j] + this->f6tp1[i][j] - this->f7tp1[i][j] - this->f8tp1[i][j])/this->rho[i][j];
            }
        }
    }


    template<class T>
    void LBM<T>::Collision() {
        for (int i = 0; i < this->nx; i++) {
            for (int j = 0; j < this->ny; j++) {
                T u2 = pow(this->u[i][j], 2.0);
                T v2 = pow(this->v[i][j], 2.0);
                T u2v2 = u2 + v2;
                T uv = this->u[i][j]*this->v[i][j];
                T omu215 = 1.0 - 1.5*u2v2;

                this->f0t[i][j] = (1.0 - this->omega)*this->f0tp1[i][j] + this->omega*this->t0*this->rho[i][j]*omu215;
                this->f1t[i][j] = (1.0 - this->omega)*this->f1tp1[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 + 3.0*this->u[i][j] + 4.5*ux2);
                this->f2t[i][j] = (1.0 - this->omega)*this->f2tp1[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 + 3.0*this->v[i][j] + 4.5*uy2);
                this->f3t[i][j] = (1.0 - this->omega)*this->f3tp1[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 - 3.0*this->u[i][j] + 4.5*ux2);
                this->f4t[i][j] = (1.0 - this->omega)*this->f4tp1[i][j] + this->omega*this->t1*this->rho[i][j]*(omu215 - 3.0*this->v[i][j] + 4.5*uy2);
                this->f5t[i][j] = (1.0 - this->omega)*this->f5tp1[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 + 3.0*(this->u[i][j] + this->v[i][j]) + 4.5*(u2v2 + 2.0*uv));
                this->f6t[i][j] = (1.0 - this->omega)*this->f6tp1[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 - 3.0*(this->u[i][j] - this->v[i][j]) + 4.5*(u2v2 - 2.0*uv));
                this->f7t[i][j] = (1.0 - this->omega)*this->f7tp1[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 - 3.0*(this->u[i][j] + this->v[i][j]) + 4.5*(u2v2 + 2.0*uv));
                this->f8t[i][j] = (1.0 - this->omega)*this->f8tp1[i][j] + this->omega*this->t2*this->rho[i][j]*(omu215 + 3.0*(this->u[i][j] - this->v[i][j]) + 4.5*(u2v2 - 2.0*uv));
            }
        }
    }
}