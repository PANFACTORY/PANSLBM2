//*****************************************************************************
//  Title       :   src/particle/d2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/21
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>

namespace PANSLBM2 {
    namespace {
        const int BARRIER = 1;
        const int MIRROR = 2;
    }
    
    template<class T>
    class D2Q9 {
public:
        D2Q9() = delete;
        D2Q9(int _nx, int _ny) : 
            nx(_nx), ny(_ny), nxy(_nx*_ny), nbc(2*(_nx + _ny)),
            offsetxmin(0), offsetxmax(_ny), offsetymin(2*_ny), offsetymax(2*_ny + _nx)
        {
            assert(0 < _nx && 0 < _ny);
            this->f = new T[this->nxy*D2Q9<T>::nc];
            this->fnext = new T[this->nxy*D2Q9<T>::nc];
            this->bctype = new int[this->nbc];
        }
        D2Q9(const D2Q9<T>& _p) = delete;
        ~D2Q9() {
            delete[] this->f, this->fnext, this->bctype;
        }

        template<class F>
        void SetBoundary(int *_bctype, F _func);
        template<class F>
        void SetBoundary(F _func) {
            D2Q9<T>::SetBoundary(this->bctype, _func);
        }
        //void SetStress(int _i, int _j, T _tx, T _ty);                       //  Set boundary condition for Elastic
        //void SetiStress(int _i, int _j, T _rho, T _tx, T _ty);              //  Set boundary condition for Adjoint of Elastic

        int Index(int _i, int _j) const {
            return _i + this->nx*_j;
        }
        int IndexStream(int _i, int _j, int _c) const {
            int ip1 = _i + D2Q9<T>::cx[_c] == -1 ? this->nx - 1 : (_i + D2Q9<T>::cx[_c] == this->nx ? 0 : _i + D2Q9<T>::cx[_c]);
            int jp1 = _j + D2Q9<T>::cy[_c] == -1 ? this->ny - 1 : (_j + D2Q9<T>::cy[_c] == this->ny ? 0 : _j + D2Q9<T>::cy[_c]);
            return this->Index(ip1, jp1);
        }
        int IndexiStream(int _i, int _j, int _c) const {
            int ip1 = _i - D2Q9<T>::cx[_c] == -1 ? this->nx - 1 : (_i - D2Q9<T>::cx[_c] == this->nx ? 0 : _i - D2Q9<T>::cx[_c]);
            int jp1 = _j - D2Q9<T>::cy[_c] == -1 ? this->ny - 1 : (_j - D2Q9<T>::cy[_c] == this->ny ? 0 : _j - D2Q9<T>::cy[_c]);
            return this->Index(ip1, jp1);
        }
        static int IndexF(int _idx, int _c) const {
            return D2Q9<T>::nc*_idx + _c;
        }

        void Swap();
        void BoundaryCondition();
        void iBoundaryCondition();
        void SmoothCorner();
        
        const int nx, ny, nxy, nbc, offsetxmin, offsetxmax, offsetymin, offsetymax;
        static const int nc = 9, nd = 2, cx[nc], cy[nc];
        static const T ei[nc];
        T *f, *fnext;

private:
        int *bctype;
    };

    template<class T>const int D2Q9<T>::cx[D2Q9<T>::nc] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    template<class T>const int D2Q9<T>::cy[D2Q9<T>::nc] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    template<class T>const T D2Q9<T>::ei[D2Q9<T>::nc] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };

    template<class T>
    template<class F>
    void D2Q9<T>::SetBoundary(int *_bctype, F _func) {
        for (int j = 0; j < this->ny; ++j) {
            this->bctype[j + this->offsetxmin] = _func(0, j);
            this->bctype[j + this->offsetxmax] = _func(this->nx - 1, j);
        }
        for (int i = 0; i < this->nx; ++i) {
            this->bctype[i + this->offsetymin] = _func(i, 0);
            this->bctype[i + this->offsetymax] = _func(i, this->ny - 1);
        }
    }

    template<class T>
    void D2Q9<T>::Swap() {
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;
    }

    template<class T>
    void D2Q9<T>::BoundaryCondition() {
        for (int j = 0; j < this->ny; ++j) {
            //  On xmin
            if (this->bctype[j + this->offsetxmin] == BARRIER) {
                int idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[j + this->offsetxmin] == MIRROR) {
                int idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 7)];
            }

            //  On xmax
            if (this->bctype[j + this->offsetxmax] == BARRIER) {
                int idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[j + this->offsetxmax] == MIRROR) {
                int idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 5)];
            }
        }

        for (int i = 0; i < this->nx; ++i) {
            //  On ymin
            if (this->bctype[i + this->offsetymin] == BARRIER) {
                int idx = this->Index(j, 0);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[i + this->offsetymin] == MIRROR) {
                int idx = this->Index(j, 0);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 7)];
            }

            //  On ymax
            if (this->bctype[i + this->offsetymax] == BARRIER) {
                int idx = this->Index(j, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[i + this->offsetymax] == MIRROR) {
                int idx = this->Index(j, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 5)];
            }
        }
    }

    template<class T>
    void D2Q9<T>::iBoundaryCondition() {
        for (int j = 0; j < this->ny; ++j) {
            //  On xmin
            if (this->bctype[j + this->offsetxmin] == BARRIER) {
                int idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[j + this->offsetxmin] == MIRROR) {
                int idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            }

            //  On xmax
            if (this->bctype[j + this->offsetxmax] == BARRIER) {
                int idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[j + this->offsetxmax] == MIRROR) {
                int idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            }
        }

        for (int i = 0; i < this->nx; ++i) {
            //  On ymin
            if (this->bctype[i + this->offsetymin] == BARRIER) {
                int idx = this->Index(j, 0);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[i + this->offsetymin] == MIRROR) {
                int idx = this->Index(j, 0);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            }

            //  On ymax
            if (this->bctype[i + this->offsetymax] == BARRIER) {
                int idx = this->Index(j, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[i + this->offsetymax] == MIRROR) {
                int idx = this->Index(j, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            }
        }
    }

    template<class T>
    void D2Q9<T>::SmoothCorner() {
        for (int c = 0; c < D2Q9<T>::nc; ++c) {
            int idx, idxx, idxy;

            //  Corner at xmin, ymin
            idx = this->Index(0, 0);
            idxx = this->Index(1, 0);
            idxy = this->Index(0, 1);
            this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);

            //  Corner at xmin, ymax
            idx = this->Index(0, this->ny - 1);
            idxx = this->Index(1, this->ny - 1);
            idxy = this->Index(0,this->ny - 2);
            this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);    

            //  Corner at xmax, ymin
            idx = this->Index(this->nx - 1, 0);
            idxx = this->Index(this->nx - 2, 0);
            idxy = this->Index(this->nx - 1, 1);
            this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  

            //  Corner at xmax, ymax
            idx = this->Index(this->nx - 1, this->ny - 1);
            idxx = this->Index(this->nx - 2, this->ny - 1);
            idxy = this->Index(this->nx - 1, this->ny - 2);
            this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  
        }
    }
/*    
    template<class T>
    void D2Q9<T>::SetStress(int _i, int _j, T _tx, T _ty) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            this->ft[1][ij] = this->ft[3][ij] - 4.0*(this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + 2.0*_tx/3.0;
            this->ft[5][ij] = this->ft[6][ij] - (this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + (_tx + 3.0*_ty)/6.0;
            this->ft[8][ij] = this->ft[7][ij] - (this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 + (_tx - 3.0*_ty)/6.0;
        } else if (_i == this->nx - 1) {
            this->ft[3][ij] = this->ft[1][ij] - 4.0*(this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - 2.0*_tx/3.0;
            this->ft[6][ij] = this->ft[5][ij] - (this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - (_tx - 3.0*_ty)/6.0;
            this->ft[7][ij] = this->ft[8][ij] - (this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 - (_tx + 3.0*_ty)/6.0;
        } else if (_j == 0) {
            this->ft[2][ij] = this->ft[4][ij] - 4.0*(this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + 2.0*_ty/3.0;
            this->ft[5][ij] = this->ft[8][ij] - (this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + (_ty + 3.0*_tx)/6.0;
            this->ft[6][ij] = this->ft[7][ij] - (this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 + (_ty - 3.0*_tx)/6.0;
        } else if (_j == this->ny - 1) {
            this->ft[4][ij] = this->ft[2][ij] - 4.0*(this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - 2.0*_ty/3.0;
            this->ft[7][ij] = this->ft[6][ij] - (this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - (_ty + 3.0*_tx)/6.0;
            this->ft[8][ij] = this->ft[5][ij] - (this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 - (_ty - 3.0*_tx)/6.0;
        } else {
            //  境界に沿っていないことを警告する
        }
    }


    template<class T>
    void D2Q9<T>::SetiStress(int _i, int _j, T _rho, T _tx, T _ty) {
        int ij = this->GetIndex(_i, _j);
        if (_i == 0) {
            this->ft[3][ij] = this->ft[1][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*_tx/_rho;
            this->ft[6][ij] = this->ft[5][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*(_tx - _ty)/_rho;
            this->ft[7][ij] = this->ft[8][ij] - (4.0*this->ft[1][ij] + this->ft[5][ij] + this->ft[8][ij])/3.0 + 2.0*(_tx + _ty)/_rho;
        } else if (_i == this->nx - 1) {
            this->ft[1][ij] = this->ft[3][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*_tx/_rho;
            this->ft[5][ij] = this->ft[6][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*(_tx + _ty)/_rho;
            this->ft[8][ij] = this->ft[7][ij] - (4.0*this->ft[3][ij] + this->ft[6][ij] + this->ft[7][ij])/3.0 - 2.0*(_tx - _ty)/_rho;
        } else if (_j == 0) {
            this->ft[4][ij] = this->ft[2][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*_ty/_rho;
            this->ft[7][ij] = this->ft[6][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*(_ty + _tx)/_rho;
            this->ft[8][ij] = this->ft[5][ij] - (4.0*this->ft[2][ij] + this->ft[5][ij] + this->ft[6][ij])/3.0 + 2.0*(_ty - _tx)/_rho;
        } else if (_j == this->ny - 1) {
            this->ft[2][ij] = this->ft[4][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*_ty/_rho;
            this->ft[5][ij] = this->ft[8][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*(_ty + _tx)/_rho;
            this->ft[6][ij] = this->ft[7][ij] - (4.0*this->ft[4][ij] + this->ft[7][ij] + this->ft[8][ij])/3.0 - 2.0*(_ty - _tx)/_rho;
        } else {
            //  境界に沿っていないことを警告する
        }
    }
*/
}