#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>

namespace {
    template<class T, class F>
    std::vector<T> solvels(std::vector<T> _A, std::vector<T> _b, F _idx){
        std::vector<T> x(_b.size());
        for(int i = 0; i < _b.size() - 1; i++){
            //----------Get pivot----------
            T pivot = fabs(_A[_idx(i, i)]);
            int pivoti = i;
            for(int j = i + 1; j < _b.size(); j++){
                if(pivot < fabs(_A[_idx(j, i)])){
                    pivot = fabs(_A[_idx(j, i)]);
                    pivoti = j;
                }
            }
            
            //----------Exchange pivot----------
            if(pivoti != i){
                T tmp = _b[i];
                _b[i] = _b[pivoti];
                _b[pivoti] = tmp;
                for(int j = i; j < _b.size(); j++){
                    tmp = _A[_idx(i, j)];
                    _A[_idx(i, j)] = _A[_idx(pivoti, j)];
                    _A[_idx(pivoti, j)] = tmp;
                }
            }
            
            //----------Forward erase----------
            for(int j = i + 1; j < _b.size(); j++){
                for(int k = i + 1; k < _b.size(); k++){
                    _A[_idx(j, k)] -= _A[_idx(i, k)]*_A[_idx(j, i)]/_A[_idx(i, i)];
                }
                _b[j] -= _b[i]*_A[_idx(j, i)]/_A[_idx(i, i)];
            }
        }
        
        //----------Back substitution----------
        for(int i = _b.size() - 1; i >= 0; i--){
            x[i] = _b[i];
            for(int j = _b.size() - 1; j > i; j--){
                x[i] -= x[j]*_A[_idx(i, j)];
            }
            x[i] /= _A[_idx(i, i)];
        }
        return x;
    }
}

namespace PANSLBM2 {
    //********************Optimizational solver with MMA********************
    template<class T>
    class MMA{
public:
        MMA() = delete;
        MMA(int _n, int _m, T _a0, std::vector<T> _a, std::vector<T> _c, std::vector<T> _d, const std::vector<T>& _xmin, const std::vector<T>& _xmax) : n(_n), m(_m) {
            //----------Initialize solver parameter----------
            this->k = 0;
            this->previousvalue = T();
            this->epsvalue = 1.0e-6;
            this->xmin = _xmin;
            this->xmax = _xmax;
            this->xkm2 = std::vector<T>(this->n, T());
            this->xkm1 = std::vector<T>(this->n, T());
            this->a0 = _a0;
            this->a = _a;
            this->c = _c;
            this->d = _d;

            //----------Set default MMA parameters----------
            this->raa0 = 1.0e-5;
            this->albefa = 0.1;
            this->move = 0.5;
            this->asyinit = 0.5;
            this->asydecr = 0.7;
            this->asyincr = 1.2;
            this->L = std::vector<T>(this->n, T());
            this->U = std::vector<T>(this->n, T());
        }
        MMA(const MMA<T>&) = delete;
        ~MMA() {}

        void ImportSetting(std::string _fname);
        void ExportSetting(std::string _fname);

        bool IsConvergence(T _currentf0) {
            return fabs((_currentf0 - this->previousvalue)/(_currentf0 + this->previousvalue)) < this->epsvalue;
        }
        void UpdateVariables(std::vector<T>& _xk, T _f, const std::vector<T>& _dfdx, const std::vector<T>& _g, const std::vector<std::vector<T> >& _dgdx);

        //----------Parameters for solver----------
        T epsvalue;             //  Convergence parameter of objective value

        //----------Parameters for MMA----------
        T raa0;                 //  Used in equation(3.3) and (3.4)
        T albefa;               //  Used in equation(3.6) and (3.7)
        T move;                 //  Used in equation(3.6) and (3.7)
        T asyinit;              //  Used in equation(3.11)
        T asydecr;              //  Used in equation(3.13)
        T asyincr;              //  Used in equation(3.13)
        
private:
        //----------Parameters for solver----------
		int k;                  //  Counter for outer loop
        const int n;            //  Number of design variables
        const int m;            //  Number of constraint  
        T previousvalue;        //  Previous function value
        std::vector<T> xmin;    //  Minimum value of design variable
        std::vector<T> xmax;    //  Maximum value of design variable
        std::vector<T> xkm2;    //  k-2 th value of design variable
        std::vector<T> xkm1;    //  k-1 th value of design variable

        //----------Parameters for MMA----------
        T a0;                   //  Used in equation(3.1)
        std::vector<T> a;       //  Used in equation(3.1)
        std::vector<T> c;       //  Used in equation(3.1)
        std::vector<T> d;       //  Used in equation(3.1)
        std::vector<T> L;       //  Parameter of asymptotes
        std::vector<T> U;       //  Parameter of asymptotes

        int IdxAm(int _i, int _j) const {
            return (this->m + 1)*_i + _j;
        }
        int IdxAn(int _i, int _j) const {
            return (this->n + 1)*_i + _j;
        }
        int IdxG(int _i, int _j) const {
            return this->n*_i + _j;
        }

		T KKTNorm(const std::vector<T>& _x, const std::vector<T>& _y, T _z, const std::vector<T>& _lambda, const std::vector<T>& _gsi, const std::vector<T>& _ita, const std::vector<T>& _mu, T _zeta, const std::vector<T>& _s,
            T _eps, const std::vector<std::vector<T> >& _p, const std::vector<std::vector<T> >& _q, const std::vector<T>& _p0, const std::vector<T>&  _q0, const std::vector<T>& _alpha, const std::vector<T>& _beta, const std::vector<T>& _b);
	};

    template<class T>
    void MMA<T>::ImportSetting(std::string _fname) {
        std::ifstream fin(_fname);
        
        fin >> this->k;
        fin >> this->previousvalue;
        fin >> this->epsvalue;
        
        fin >> this->raa0;
        fin >> this->albefa;
        fin >> this->move;
        fin >> this->asyinit;
        fin >> this->asydecr;
        fin >> this->asyincr;

        fin >> this->a0;

        for (int i = 0; i < this->m; ++i) {
            fin >> this->a[i] >> this->c[i] >> this->d[i];
        }

        for (int j = 0; j < this->n; ++j) {
            fin >> this->xmin[j] >> this->xmax[j] >> this->xkm2[j] >> this->xkm1[j] >> this->L[j] >> this->U[j];
        }
    }

    template<class T>
    void MMA<T>::ExportSetting(std::string _fname) {
        std::ofstream fout(_fname);

        fout << this->k << std::endl;
        fout << this->previousvalue << std::endl;
        fout << this->epsvalue << std::endl;
        
        fout << this->raa0 << std::endl;
        fout << this->albefa << std::endl;
        fout << this->move << std::endl;
        fout << this->asyinit << std::endl;
        fout << this->asydecr << std::endl;
        fout << this->asyincr << std::endl;

        fout << this->a0 << std::endl;

        for (int i = 0; i < this->m; ++i) {
            fout << this->a[i] << "\t" << this->c[i] << "\t" << this->d[i] << std::endl;
        }

        for (int j = 0; j < this->n; ++j) {
            fout << this->xmin[j] << "\t" << this->xmax[j] << "\t" << this->xkm2[j] << "\t" << this->xkm1[j] << "\t" << this->L[j] << "\t" << this->U[j] << std::endl;
        }
    }

    template<class T>
    void MMA<T>::UpdateVariables(std::vector<T>& _xk, T _f, const std::vector<T>& _dfdx, const std::vector<T>& _g, const std::vector<std::vector<T> >& _dgdx){       
		//----------Get asymptotes parameter L and U----------
		if(this->k < 2){
			for(int j = 0; j < this->n; j++){
				this->L[j] = _xk[j] - this->asyinit*(this->xmax[j] - this->xmin[j]);
				this->U[j] = _xk[j] + this->asyinit*(this->xmax[j] - this->xmin[j]);
                this->L[j] = std::min(std::max(_xk[j] - 10.0*(this->xmax[j] - this->xmin[j]), this->L[j]), _xk[j] - 0.01*(this->xmax[j] - this->xmin[j]));
				this->U[j] = std::min(std::max(_xk[j] + 0.01*(this->xmax[j] - this->xmin[j]), this->U[j]), _xk[j] + 10.0*(this->xmax[j] - this->xmin[j]));
            }
		} else {
			for(int j = 0; j < this->n; j++){
				T tmp = (_xk[j] - this->xkm1[j])*(this->xkm1[j] - this->xkm2[j]);
				if(tmp < T()){
					this->L[j] = _xk[j] - this->asydecr*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->asydecr*(this->U[j] - this->xkm1[j]);
				} else if(tmp > T()){
					this->L[j] = _xk[j] - this->asyincr*(this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + this->asyincr*(this->U[j] - this->xkm1[j]);
				} else{
					this->L[j] = _xk[j] - (this->xkm1[j] - this->L[j]);
					this->U[j] = _xk[j] + (this->U[j] - this->xkm1[j]);
				}
                this->L[j] = std::min(std::max(_xk[j] - 10.0*(this->xmax[j] - this->xmin[j]), this->L[j]), _xk[j] - 0.01*(this->xmax[j] - this->xmin[j]));
				this->U[j] = std::min(std::max(_xk[j] + 0.01*(this->xmax[j] - this->xmin[j]), this->U[j]), _xk[j] + 10.0*(this->xmax[j] - this->xmin[j]));
			}
		}

		//----------Get movelimit at k----------
		std::vector<T> alpha(this->n);
        std::vector<T> beta(this->n);
		for(int j = 0; j < this->n; j++){
			alpha[j] = std::max({this->xmin[j], this->L[j] + this->albefa*(_xk[j] - this->L[j]), _xk[j] - this->move*(this->xmax[j] - this->xmin[j])});
			beta[j] = std::min({this->xmax[j], this->U[j] - this->albefa*(this->U[j] - _xk[j]), _xk[j] + this->move*(this->xmax[j] - this->xmin[j])});
		}

		//----------Get p0 and q0----------
		std::vector<T> p0(this->n);
        std::vector<T> q0(this->n);
		for(int j = 0; j < this->n; j++){
			T dfdxp = std::max(_dfdx[j], T());
			T dfdxm = std::max(-_dfdx[j], T());
			p0[j] = pow(this->U[j] - _xk[j], 2.0)*(1.001*dfdxp + 0.001*dfdxm + this->raa0/(this->xmax[j] - this->xmin[j]));
			q0[j] = pow(_xk[j] - this->L[j], 2.0)*(0.001*dfdxp + 1.001*dfdxm + this->raa0/(this->xmax[j] - this->xmin[j]));
		}

		//----------Get p, q and b----------
		std::vector<std::vector<T> > p(this->m, std::vector<T>(this->n));
        std::vector<std::vector<T> > q(this->m, std::vector<T>(this->n));
        std::vector<T> b(this->m), b_buffer(this->m, T());
		for(int i = 0; i < this->m; i++){
			for(int j = 0; j < this->n; j++){
				T dfdxp = std::max(_dgdx[i][j], T());
				T dfdxm = std::max(-_dgdx[i][j], T());
				p[i][j] = pow(this->U[j] - _xk[j], 2.0)*(1.001*dfdxp + 0.001*dfdxm + this->raa0/(this->xmax[j] - this->xmin[j]));
				q[i][j] = pow(_xk[j] - this->L[j], 2.0)*(0.001*dfdxp + 1.001*dfdxm + this->raa0/(this->xmax[j] - this->xmin[j]));
				b_buffer[i] += p[i][j]/(this->U[j] - _xk[j]) + q[i][j]/(_xk[j] - this->L[j]);
			}
		}
        MPI_Allreduce(b_buffer.data(), b.data(), this->m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for(int i = 0; i < this->m; i++){
			b[i] -= _g[i];
        }
     
        //----------Inner loop----------
        T eps = 1.0;
        std::vector<T> x(this->n);
        std::vector<T> y(this->m, 1.0);
        T z = 1.0;
        T zeta = 1.0;
        std::vector<T> lambda(this->m, 1.0);
        std::vector<T> s(this->m, 1.0);
        std::vector<T> gsi(this->n);
        std::vector<T> ita(this->n);
        std::vector<T> mu(this->m);

        for(int i = 0; i < this->m; i++){
            mu[i] = std::max(1.0, 0.5*this->c[i]);
        }

        for(int j = 0; j < this->n; j++){
            x[j] = 0.5*(alpha[j] + beta[j]);
            gsi[j] = std::max(1.0, 1.0/(x[j] - alpha[j]));
            ita[j] = std::max(1.0, 1.0/(beta[j] - x[j]));
        }

        std::vector<T> plambda(this->n);
        std::vector<T> qlambda(this->n);          
        std::vector<T> G(this->m*this->n);
        std::vector<T> Dx(this->n);
        std::vector<T> deltilx(this->n);
        std::vector<T> Dy(this->m);
        std::vector<T> Dlambda(this->m);
        std::vector<T> deltily(this->m);
        std::vector<T> deltillambda(this->m), deltillambda_buffer(this->m);;
        std::vector<T> Dlambday(this->m);
        std::vector<T> deltillambday(this->m);
        
        std::vector<T> dx(this->n);
        std::vector<T> dy(this->m);
        T dz = T();
        std::vector<T> dlambda(this->m), dlambda_buffer(this->m);
        std::vector<T> dgsi(this->n);
        std::vector<T> dita(this->n);
        std::vector<T> dmu(this->m);
        T dzeta = T();
        std::vector<T> ds(this->m);

        std::vector<T> xpdx(this->n);
        std::vector<T> ypdy(this->m);
        T zpdz;
        std::vector<T> lambdapdlambda(this->m);
        std::vector<T> gsipdgsi(this->n);
        std::vector<T> itapdita(this->n);
        std::vector<T> mupdmu(this->m);
        T zetapdzeta;
        std::vector<T> spds(this->m);
            
        for(int l = 0; eps > 1.0e-7; l++){
            //.....Get coefficients.....
            for(int j = 0; j < this->n; j++){
                plambda[j] = p0[j];
                qlambda[j] = q0[j];
                for(int i = 0; i < this->m; i++){
                    plambda[j] += lambda[i]*p[i][j];
                    qlambda[j] += lambda[i]*q[i][j];
                }
            }

            for(int i = 0; i < this->m; i++){
                for(int j = 0; j < this->n; j++){
                    G[this->IdxG(i, j)] = p[i][j]/pow(this->U[j] - x[j], 2.0) - q[i][j]/pow(x[j] - this->L[j], 2.0);
                }
            }

            for(int j = 0; j < this->n; j++){
                Dx[j] = 2.0*plambda[j]/pow(this->U[j] - x[j], 3.0) + 2.0*qlambda[j]/pow(x[j] - this->L[j], 3.0) + gsi[j]/(x[j] - alpha[j]) + ita[j]/(beta[j] - x[j]);
                deltilx[j] = plambda[j]/pow(this->U[j] - x[j], 2.0) - qlambda[j]/pow(x[j] - this->L[j], 2.0) - eps/(x[j] - alpha[j]) + eps/(beta[j] - x[j]);
            }

            for(int i = 0; i < this->m; i++){
                Dy[i] = d[i] + mu[i]/y[i];
                Dlambda[i] = s[i]/lambda[i];
                deltily[i] = c[i] + this->d[i]*y[i] - lambda[i] - eps/y[i];
            }

            T deltilz = this->a0 - eps/z - std::inner_product(lambda.begin(), lambda.end(), this->a.begin(), T());
            
            for (int i = 0; i < this->m; ++i) {
                deltillambda_buffer[i] = T();
                for (int j = 0; j < this->n; j++) {
                    deltillambda_buffer[i] += p[i][j]/(this->U[j] - x[j]) + q[i][j]/(x[j] - this->L[j]);
                }
            }
            MPI_Allreduce(deltillambda_buffer.data(), deltillambda.data(), this->m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            for (int i = 0; i < this->m; i++) {
                deltillambda[i] += -this->a[i]*z - y[i] - b[i] + eps/lambda[i];
                Dlambday[i] = Dlambda[i] + 1.0/Dy[i];
                deltillambday[i] = deltillambda[i] + deltily[i]/Dy[i];
            }

            //.....Get Newton direction.....
            if(this->n > this->m){
                std::vector<T> A((this->m + 1)*(this->m + 1), T()), A_buffer((this->m + 1)*(this->m + 1), T());
                for(int ii = 0; ii < this->m; ii++){
                    for(int jj = 0; jj < this->m; jj++){
                        for(int kk= 0; kk < this->n; kk++){
                            A_buffer[this->IdxAm(ii, jj)] += G[this->IdxG(ii, kk)]*G[this->IdxG(jj, kk)]/Dx[kk];
                        }
                    }
                }
                MPI_Allreduce(A_buffer.data(), A.data(), (this->m + 1)*(this->m + 1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int ii = 0; ii < this->m; ii++){
                    A[this->IdxAm(ii, ii)] += Dlambday[ii];
                    A[this->IdxAm(ii, this->m)] = this->a[ii];
                    A[this->IdxAm(this->m, ii)] = this->a[ii];
                }
                A[this->IdxAm(this->m, this->m)] = -zeta/z;
                std::vector<T> B(this->m + 1), B_buffer(this->m, T());
                for(int ii = 0; ii < this->m; ii++){
                    for(int jj = 0; jj < this->n; jj++){
                        B_buffer[ii] -= G[this->IdxG(ii, jj)]*deltilx[jj]/Dx[jj];
                    }
                }
                MPI_Allreduce(B_buffer.data(), B.data(), this->m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int ii = 0; ii < this->m; ii++){
                    B[ii] += deltillambday[ii];
                }
                B[this->m] = deltilz;
                std::vector<T> dlambdaz = solvels<T>(A, B, [&](int _i, int _j) { return this->IdxAm(_i, _j); });
                for(int i = 0; i < this->m; i++){
                    dlambda[i] = dlambdaz[i];
                }
                dz = dlambdaz[this->m];
                for(int j = 0; j < this->n; j++){
                    dx[j] = -deltilx[j]/Dx[j];
                    for(int i = 0; i < this->m; i++){
                        dx[j] -= G[this->IdxG(i, j)]*dlambda[i]/Dx[j];
                    }
                }
            } else {
                std::vector<T> A((this->n + 1)*(this->n + 1), T());             
                for(int ii = 0; ii < this->n; ii++){
                    for(int jj = 0; jj < this->n; jj++){
                        for(int kk = 0; kk < this->m; kk++){
                            A[this->IdxAn(ii, jj)] += G[this->IdxG(kk, ii)]*G[this->IdxG(kk, jj)]/Dlambday[kk];
                        }
                    }
                    A[this->IdxAn(ii, ii)] += Dx[ii];
                    for(int jj = 0; jj < this->m; jj++){
                        A[this->IdxAn(ii, this->n)] -= G[this->IdxG(jj, ii)]*this->a[jj]/Dlambday[jj];
                        A[this->IdxAn(this->n, ii)] -= G[this->IdxG(jj, ii)]*this->a[jj]/Dlambday[jj];
                        A[this->IdxAn(this->n, this->n)] += this->a[jj]*this->a[jj]/Dlambday[jj]; 
                    }
                }
                A[this->IdxAn(this->n, this->n)] += zeta/z;
                std::vector<T> B(this->n + 1);
                for(int ii = 0; ii < this->n; ii++){
                    B[ii] = -deltilx[ii];
                    for(int jj = 0; jj < this->m; jj++){
                        B[ii] -= G[this->IdxG(jj, ii)]*deltillambday[jj]/Dlambday[jj];
                    }
                }
                B[this->n] = -deltilz;
                for(int jj = 0; jj < this->m; jj++){
                    B[this->n] += this->a[jj]*deltillambday[jj]/Dlambday[jj];
                }
                std::vector<T> dxz = solvels<T>(A, B, [&](int _i, int _j) { return this->IdxAn(_i, _j); });
                for(int j = 0; j < this->n; j++){
                    dx[j] = dxz[j];
                }
                dz = dxz[this->n];
                for(int i = 0; i < this->m; i++){
                    dlambda_buffer[i] = T();
                    for(int j = 0; j < this->n; j++){
                        dlambda_buffer[i] += G[this->IdxG(i, j)]*dx[j]/Dlambday[i];
                    }
                }
                MPI_Allreduce(dlambda_buffer.data(), dlambda.data(), this->m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                for(int i = 0; i < this->m; i++){
                    dlambda[i] += -this->a[i]*dz/Dlambday[i] + deltillambday[i]/Dlambday[i];
                }
            }
            
            for(int i = 0; i < this->m; i++){
                dy[i] = dlambda[i]/Dy[i] - deltily[i]/Dy[i];
                dmu[i] = -mu[i]*dy[i]/y[i] - mu[i] + eps/y[i];
                ds[i] = -s[i]*dlambda[i]/lambda[i] - s[i] + eps/lambda[i];
            }

            for(int j = 0; j < this->n; j++){
                dgsi[j] = -gsi[j]*dx[j]/(x[j] - alpha[j]) - gsi[j] + eps/(x[j] - alpha[j]);
                dita[j] = ita[j]*dx[j]/(beta[j] - x[j]) - ita[j] + eps/(beta[j] - x[j]);
            }

            dzeta = -zeta*dz/z - zeta + eps/z;

            //.....Get step size.....
            T txmax, txmax_buffer = T();
            for(int j = 0; j < this->n; j++){
                T txmaxj = std::max({-1.01*dx[j]/(x[j] - alpha[j]), 1.01*dx[j]/(beta[j] - x[j]), -1.01*dgsi[j]/gsi[j], -1.01*dita[j]/ita[j]});
                if(txmax_buffer < txmaxj){
                    txmax_buffer = txmaxj;
                }
            }
            MPI_Allreduce(&txmax_buffer, &txmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            T tymax = T();
            for(int i = 0; i < this->m; i++){
                T tymaxi = std::max({-1.01*dy[i]/y[i], -1.01*dlambda[i]/lambda[i], -1.01*dmu[i]/mu[i], -1.01*ds[i]/s[i]});
                if(tymax < tymaxi){
                    tymax = tymaxi;
                }
            }
            T tau = 1.0/std::max({1.0, txmax, tymax, -1.01*dz/z, -1.01*dzeta/zeta});
            T deltawl = this->KKTNorm(x, y, z, lambda, gsi, ita, mu, zeta, s, eps, p, q, p0, q0, alpha, beta, b);
            for(int ll = 0; ll < 50; ll++){
                for(int j = 0; j < this->n; j++){
                    xpdx[j] = x[j] + tau*dx[j];
                    gsipdgsi[j] = gsi[j] + tau*dgsi[j];
                    itapdita[j] = ita[j] + tau*dita[j];
                }
                for(int i = 0; i < this->m; i++){
                    ypdy[i] = y[i] + tau*dy[i];
                    lambdapdlambda[i] = lambda[i] + tau*dlambda[i];
                    mupdmu[i] = mu[i] + tau*dmu[i];
                    spds[i] = s[i] + tau*ds[i];
                }
                zpdz = z + tau*dz;
                zetapdzeta = zeta + tau*dzeta;
                
                T deltawlp1 = this->KKTNorm(xpdx, ypdy, zpdz, lambdapdlambda, gsipdgsi, itapdita, mupdmu, zetapdzeta, spds, eps, p, q, p0, q0, alpha, beta, b);
                if(deltawlp1 < deltawl){
                    break;
                }
                tau *= 0.5;
            }

            //.....Update w.....
            x = xpdx;
            y = ypdy;
            z = zpdz;
            lambda = lambdapdlambda;
            gsi = gsipdgsi;
            ita = itapdita;
            mu = mupdmu;
            zeta = zetapdzeta;
            s = spds;

            //.....Update epsl.....
            T deltawlp1 = this->KKTNorm(x, y, z, lambda, gsi, ita, mu, zeta, s, eps, p, q, p0, q0, alpha, beta, b);
            eps *= deltawlp1 < 0.9*eps ? 0.1 : 1.0;
        }
        
        //----------Update outer loop counter k----------
        this->previousvalue = _f;
        this->k++;
        this->xkm2 = this->xkm1;
        this->xkm1 = _xk;
        _xk = x; 
    }

    template<class T>
    T MMA<T>::KKTNorm(
        const std::vector<T>& _x, const std::vector<T>& _y, T _z, 
        const std::vector<T>& _lambda, const std::vector<T>& _gsi, const std::vector<T>& _ita, const std::vector<T>& _mu, T _zeta, const std::vector<T>& _s,
        T _eps, const std::vector<std::vector<T> >& _p, const std::vector<std::vector<T> >& _q, const std::vector<T>& _p0, const std::vector<T>&  _q0, 
        const std::vector<T>& _alpha, const std::vector<T>& _beta, const std::vector<T>& _b
    ){
        T norm = T(), norm_buffer = T();

        //----------Get parameters----------
        std::vector<T> plambda(this->n);
        std::vector<T> qlambda(this->n);
        std::vector<T> g(this->m, T()), g_buffer(this->m, T());
        for(int j = 0; j < this->n; j++){
            plambda[j] = _p0[j];
            qlambda[j] = _q0[j];
            for(int i = 0; i < this->m; i++){
                plambda[j] += _lambda[i]*_p[i][j];
                qlambda[j] += _lambda[i]*_q[i][j];
                g_buffer[i] += _p[i][j]/(this->U[j] - _x[j]) + _q[i][j]/(_x[j] - this->L[j]);
            }
        }
        MPI_Allreduce(g_buffer.data(), g.data(), this->m, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //----------Equation(5.9a)(5.9e)(5.9f)----------
        for(int j = 0; j < this->n; j++){
            norm_buffer += pow(plambda[j]/pow(this->U[j] - _x[j], 2.0) - qlambda[j]/pow(_x[j] - this->L[j], 2.0) - _gsi[j] + _ita[j], 2.0); //  Equation(5.9a)
            norm_buffer += pow(_gsi[j]*(_x[j] - _alpha[j]) - _eps, 2.0);    //  Equation(5.9e)
            norm_buffer += pow(_ita[j]*(_beta[j] - _x[j]) - _eps, 2.0);     //  Equation(5.9f)
        }

        //----------Equation(5.9b)(5.9d)(5.9g)(5.9i)----------
        for(int i = 0; i < this->m; i++){
            norm_buffer += pow(this->c[i] + this->d[i]*_y[i] - _lambda[i] - _mu[i], 2.0);   //  Equation(5.9b)
            norm_buffer += pow(g[i] - this->a[i]*_z - _y[i] + _s[i] - _b[i], 2.0);          //  Equation(5.9d)
            norm_buffer += pow(_mu[i]*_y[i] - _eps, 2.0);                                   //  Equation(5.9g)
            norm_buffer += pow(_lambda[i]*_s[i] - _eps, 2.0);                               //  Equation(5.9i)
        }

        //----------Equation(5.9c)(5.9h)----------
        norm_buffer += pow(this->a0 - _zeta - std::inner_product(_lambda.begin(), _lambda.end(), this->a.begin(), T()), 2.0);   //  Equation(5.9c)
        norm_buffer += pow(_zeta*_z - _eps, 2.0);                                                                               //  Equation(5.9h)

        MPI_Allreduce(&norm_buffer, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(norm);
    }
}