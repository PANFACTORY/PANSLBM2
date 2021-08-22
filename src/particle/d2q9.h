//*****************************************************************************
//  Title       :   src/particle/d2q9mpi.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/8/14
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include <cmath>
#include <cassert>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

namespace PANSLBM2 {
    namespace {
        const int BARRIER = 1;
        const int MIRROR = 2;
    }
    
    template<class T>
    class D2Q9 {
public:
        D2Q9() = delete;
        D2Q9(int _lx, int _ly, int _PEid = 0, int _mx = 1, int _my = 1) :
            lx(_lx), ly(_ly), PEid(_PEid), mx(_mx), my(_my), 
            PEx(this->PEid%this->mx), PEy(this->PEid/this->mx),
            nx((this->lx + this->PEx)/this->mx), ny((this->ly + this->PEy)/this->my),
            nxy(this->nx*this->ny), nbc(2*(this->nx + this->ny)),
            offsetxmin(0), offsetxmax(this->ny), offsetymin(2*this->ny), offsetymax(2*this->ny + this->nx),
            offsetx(this->mx - this->PEx > this->lx%this->mx ? this->PEx*this->nx : this->lx - (this->mx - this->PEx)*this->nx),
            offsety(this->my - this->PEy > this->ly%this->my ? this->PEy*this->ny : this->ly - (this->my - this->PEy)*this->ny)
        {
            assert(0 < _lx && 0 < _ly && 0 <= _PEid && 0 < _mx && 0 < _my);

            this->f = new T[this->nxy*D2Q9<T>::nc];
            this->fnext = new T[this->nxy*D2Q9<T>::nc];
            this->fsend = new T[this->nbc*3 + 4];
            this->frecv = new T[this->nbc*3 + 4];
            this->bctype = new int[this->nbc];
            for (int idx = 0; idx < this->nbc; ++idx) {
                this->bctype[idx] = 0;
            }

#ifdef _USE_MPI_DEFINES
            this->StatSend = new MPI_Status[8];
            this->StatRecv = new MPI_Status[8];
            this->ReqSend = new MPI_Request[8];
            this->ReqRecv = new MPI_Request[8];
#endif
        }
        D2Q9(const D2Q9<T>& _p) = delete;
        ~D2Q9() {
            delete[] this->f, this->fnext, this->fsend, this->frecv, this->bctype;
#ifdef _USE_MPI_DEFINES
            delete[] this->StatSend, this->StatRecv, this->ReqSend, this->ReqRecv;
#endif
        }

        template<class F>
        void SetBoundary(int *_bctype, F _func);
        template<class F>
        void SetBoundary(F _func) {
            this->SetBoundary(this->bctype, _func);
        }
        
        int Index(int _i, int _j) const {
            int i = _i == -1 ? this->nx - 1 : (_i == this->nx ? 0 : _i);
            int j = _j == -1 ? this->ny - 1 : (_j == this->ny ? 0 : _j);
            return i + this->nx*j;
        }
        static int IndexF(int _idx, int _c) {
            return D2Q9<T>::nc*_idx + _c;
        }
        int IndexPE(int _i, int _j) const {
            int i = _i == -1 ? this->mx - 1 : (_i == this->mx ? 0 : _i);
            int j = _j == -1 ? this->my - 1 : (_j == this->my ? 0 : _j);
            return i + this->mx*j;
        }

        void Swap();
        void BoundaryCondition();
        void iBoundaryCondition();
        void SmoothCorner();
        
        void Synchronize();
        void iSynchronize();
        
        const int lx, ly, PEid, mx, my, PEx, PEy, nx, ny, nxy, nbc, offsetxmin, offsetxmax, offsetymin, offsetymax, offsetx, offsety;
        static const int nc = 9, nd = 2, cx[nc], cy[nc];
        static const T ei[nc];
        T *f, *fnext, *fsend, *frecv;

private:
        int *bctype;
#ifdef _USE_MPI_DEFINES
        MPI_Status *StatSend, *StatRecv;
        MPI_Request *ReqSend, *ReqRecv;
#endif
    };

    template<class T>const int D2Q9<T>::cx[D2Q9<T>::nc] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    template<class T>const int D2Q9<T>::cy[D2Q9<T>::nc] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    template<class T>const T D2Q9<T>::ei[D2Q9<T>::nc] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };

    template<class T>
    template<class F>
    void D2Q9<T>::SetBoundary(int *_bctype, F _func) {
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                _bctype[j + this->offsetxmin] = _func(0 + this->offsetx, j + this->offsety);
            }
        }
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                _bctype[j + this->offsetxmax] = _func((this->nx - 1) + this->offsetx, j + this->offsety);
            }
        }
        if (this->PEy == 0) {
            for (int i = 0; i < this->nx; ++i) {
                _bctype[i + this->offsetymin] = _func(i + this->offsetx, 0 + this->offsety);
            }
        }
        if (this->PEy == this->my - 1) {
            for (int i = 0; i < this->nx; ++i) {
                _bctype[i + this->offsetymax] = _func(i + this->offsetx, (this->ny - 1) + this->offsety);
            }
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
                int idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[i + this->offsetymin] == MIRROR) {
                int idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 7)];
            }

            //  On ymax
            if (this->bctype[i + this->offsetymax] == BARRIER) {
                int idx = this->Index(i, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[i + this->offsetymax] == MIRROR) {
                int idx = this->Index(i, this->ny - 1);
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
                int idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            } else if (this->bctype[i + this->offsetymin] == MIRROR) {
                int idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 6)];
            }

            //  On ymax
            if (this->bctype[i + this->offsetymax] == BARRIER) {
                int idx = this->Index(i, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            } else if (this->bctype[i + this->offsetymax] == MIRROR) {
                int idx = this->Index(i, this->ny - 1);
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
            if (this->PEx == 0 && this->PEy == 0) {
                idx = this->Index(0, 0);
                idxx = this->Index(1, 0);
                idxy = this->Index(0, 1);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);
            }
            
            //  Corner at xmin, ymax
            if (this->PEx == 0 && this->PEy == this->my - 1) {
               idx = this->Index(0, this->ny - 1);
                idxx = this->Index(1, this->ny - 1);
                idxy = this->Index(0,this->ny - 2);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);    
            }
            
            //  Corner at xmax, ymin
            if (this->PEx == this->mx - 1 && this->PEy == 0) {
                idx = this->Index(this->nx - 1, 0);
                idxx = this->Index(this->nx - 2, 0);
                idxy = this->Index(this->nx - 1, 1);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  
            }
            
            //  Corner at xmax, ymax
            if (this->PEx == this->mx - 1 && this->PEy == this->my - 1) {
                idx = this->Index(this->nx - 1, this->ny - 1);
                idxx = this->Index(this->nx - 2, this->ny - 1);
                idxy = this->Index(this->nx - 1, this->ny - 2);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  
            }
        }
    }

    template<class T>
    void D2Q9<T>::Synchronize() {
#ifdef _USE_MPI_DEFINES
        int idx, idxedge;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin  
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmin;
                this->fsend[idxedge*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->fsend[idxedge*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->fsend[idxedge*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 7)];

                //  Edge along xmax
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmax;
                this->fsend[idxedge*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->fsend[idxedge*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[idxedge*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)]; 
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymin;
                this->fsend[idxedge*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->fsend[idxedge*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->fsend[idxedge*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)]; 

                //  Edge along ymax
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymax;
                this->fsend[idxedge*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->fsend[idxedge*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[idxedge*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 6)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->fsend[this->nbc*3 + 0] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 7)];    //  Corner at xmin and ymin
            this->fsend[this->nbc*3 + 1] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 6)];               //  Corner at xmin and ymax
            this->fsend[this->nbc*3 + 2] = this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 8)];               //  Corner at xmax and ymin
            this->fsend[this->nbc*3 + 3] = this->f[D2Q9<T>::IndexF(this->Index(0, 0), 5)];                          //  Corner at xmax and ymax
        }
        
        //  Communicate with other PE
        int neib = 0;
        if (this->mx != 1) {
            //  Left
            MPI_Isend(&fsend[this->offsetxmin*3], this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 0, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->offsetxmax*3], this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 0, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Right
            MPI_Isend(&fsend[this->offsetxmax*3], this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 1, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->offsetxmin*3], this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 1, MPI_COMM_WORLD, &this->ReqRecv[neib++]);
        }
        if (this->my != 1) {
            //  Bottom
            MPI_Isend(&fsend[this->offsetymin*3], this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 2, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->offsetymax*3], this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 2, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Top
            MPI_Isend(&fsend[this->offsetymax*3], this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 3, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->offsetymin*3], this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 3, MPI_COMM_WORLD, &this->ReqRecv[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  Left bottom
            MPI_Isend(&fsend[this->nbc*3 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 4, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->nbc*3 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 4, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Left top
            MPI_Isend(&fsend[this->nbc*3 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 5, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[this->nbc*3 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 5, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Right bottom
            MPI_Irecv(&frecv[this->nbc*3 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 6, MPI_COMM_WORLD, &this->ReqRecv[neib]);
            MPI_Isend(&fsend[this->nbc*3 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 6, MPI_COMM_WORLD, &this->ReqSend[neib++]);

            //  Right top
            MPI_Irecv(&frecv[this->nbc*3 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 7, MPI_COMM_WORLD, &this->ReqRecv[neib]);
            MPI_Isend(&fsend[this->nbc*3 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 7, MPI_COMM_WORLD, &this->ReqSend[neib++]);      
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->ReqSend, this->StatSend);
            MPI_Waitall(neib, this->ReqRecv, this->StatRecv);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmin;
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->frecv[idxedge*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[idxedge*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[idxedge*3 + 2]; 
            
                //  Edge along xmax
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmax;
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->frecv[idxedge*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[idxedge*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[idxedge*3 + 2];
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymin;
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->frecv[idxedge*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[idxedge*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[idxedge*3 + 2];
            
                //  Edge along ymax
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymax;
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->frecv[idxedge*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[idxedge*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[idxedge*3 + 2];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->f[D2Q9<T>::IndexF(this->Index(0, 0), 5)] = this->frecv[this->nbc*3 + 0];                          //  Corner at xmin and ymin
            this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 8)] = this->frecv[this->nbc*3 + 1];               //  Corner at xmin and ymax
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 6)] = this->frecv[this->nbc*3 + 2];               //  Corner at xmax and ymin
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 7)] = this->frecv[this->nbc*3 + 3];    //  Corner at xmax and ymax
        }
#endif
    }

    template<class T>
    void D2Q9<T>::iSynchronize() {
#ifdef _USE_MPI_DEFINES
        int idx, idxedge;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin  
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmin;
                this->fsend[D2Q9<T>::IndexF(idxedge, 1)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 5)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 8)] = this->f[D2Q9<T>::IndexF(idx, 8)];

                //  Edge along xmax
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmax;
                this->fsend[D2Q9<T>::IndexF(idxedge, 3)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 6)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 7)] = this->f[D2Q9<T>::IndexF(idx, 7)]; 
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymin;
                this->fsend[D2Q9<T>::IndexF(idxedge, 2)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 5)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 6)] = this->f[D2Q9<T>::IndexF(idx, 6)]; 

                //  Edge along ymax
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymax;
                this->fsend[D2Q9<T>::IndexF(idxedge, 4)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 7)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 8)] = this->f[D2Q9<T>::IndexF(idx, 8)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->fsend[D2Q9<T>::IndexF(this->nbc + 0, 5)] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 5)];  //  Corner at xmin and ymin
            this->fsend[D2Q9<T>::IndexF(this->nbc + 1, 8)] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 8)];             //  Corner at xmin and ymax
            this->fsend[D2Q9<T>::IndexF(this->nbc + 2, 6)] = this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 6)];             //  Corner at xmax and ymin
            this->fsend[D2Q9<T>::IndexF(this->nbc + 3, 7)] = this->f[D2Q9<T>::IndexF(this->Index(0, 0), 7)];                        //  Corner at xmax and ymax
        }
        
        //  Communicate with other PE
        int neib = 0;
        if (this->mx != 1) {
            //  Left
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetxmin, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 0, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetxmax, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 0, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Right
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetxmax, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 1, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetxmin, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 1, MPI_COMM_WORLD, &this->ReqRecv[neib++]);
        }
        if (this->my != 1) {
            //  Bottom
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetymin, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 2, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetymax, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 2, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Top
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetymax, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 3, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetymin, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 3, MPI_COMM_WORLD, &this->ReqRecv[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  Left bottom
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 0, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 4, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 3, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 4, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Left top
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 1, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 5, MPI_COMM_WORLD, &this->ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 2, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 5, MPI_COMM_WORLD, &this->ReqRecv[neib++]);

            //  Right bottom
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 1, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 6, MPI_COMM_WORLD, &this->ReqRecv[neib]);
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 2, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 6, MPI_COMM_WORLD, &this->ReqSend[neib++]);

            //  Right top
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 0, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 7, MPI_COMM_WORLD, &this->ReqRecv[neib]);
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 3, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 7, MPI_COMM_WORLD, &this->ReqSend[neib++]);      
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->ReqSend, this->StatSend);
            MPI_Waitall(neib, this->ReqRecv, this->StatRecv);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmin;
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->frecv[D2Q9<T>::IndexF(idxedge, 3)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[D2Q9<T>::IndexF(idxedge, 6)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[D2Q9<T>::IndexF(idxedge, 7)]; 
            
                //  Edge along xmax
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmax;
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->frecv[D2Q9<T>::IndexF(idxedge, 1)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[D2Q9<T>::IndexF(idxedge, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[D2Q9<T>::IndexF(idxedge, 8)];
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymin;
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->frecv[D2Q9<T>::IndexF(idxedge, 4)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[D2Q9<T>::IndexF(idxedge, 7)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[D2Q9<T>::IndexF(idxedge, 8)];
            
                //  Edge along ymax
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymax;
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->frecv[D2Q9<T>::IndexF(idxedge, 2)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[D2Q9<T>::IndexF(idxedge, 5)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[D2Q9<T>::IndexF(idxedge, 6)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->f[D2Q9<T>::IndexF(this->Index(0, 0), 7)] = this->frecv[D2Q9<T>::IndexF(this->nbc + 0, 7)];                        //  Corner at xmin and ymin
            this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 6)] = this->frecv[D2Q9<T>::IndexF(this->nbc + 1, 6)];             //  Corner at xmin and ymax
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 8)] = this->frecv[D2Q9<T>::IndexF(this->nbc + 2, 8)];             //  Corner at xmax and ymin
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 5)] = this->frecv[D2Q9<T>::IndexF(this->nbc + 3, 5)];  //  Corner at xmax and ymax
        }
#endif
    }
}