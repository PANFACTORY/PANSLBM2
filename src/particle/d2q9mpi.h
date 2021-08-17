//*****************************************************************************
//  Title       :   src/particle/d2q9mpi.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/8/14
//  Copyright   :   (C)2021 TanabeYuta
//*****************************************************************************

#pragma once
#include "mpi.h"
#include <cmath>
#include <cassert>

namespace PANSLBM2 {
    namespace {
        const int BARRIER = 1;
        const int MIRROR = 2;
    }
    
    template<class T>
    class D2Q9 {
public:
        D2Q9() = delete;
        D2Q9(int _nx, int _ny, int _mpiid = 0) : 
            nx(_nx), ny(_ny), nxy(_nx*_ny), nbc(2*(_nx + _ny)),
            offsetxmin(0), offsetxmax(_ny), offsetymin(2*_ny), offsetymax(2*_ny + _nx),
            mpiid(_mpiid)
        {
            assert(0 < _nx && 0 < _ny);
            this->f = new T[this->nxy*D2Q9<T>::nc];
            this->fnext = new T[this->nxy*D2Q9<T>::nc];
            this->fsend = new T[(this->nbc + 4)*D2Q9<T>::nc];
            this->frecv = new T[(this->nbc + 4)*D2Q9<T>::nc];
            this->bctype = new int[this->nbc];

            this->neighbornum = 0;
            this->StatSend = nullptr;
            this->StatRecv = nullptr;
            this->ReqSend = nullptr;
            this->ReqRecv = nullptr;
        }
        D2Q9(const D2Q9<T>& _p) = delete;
        ~D2Q9() {
            delete[] this->f, this->fnext, this->fsend, this->frecv, this->bctype;
            if (this->neighbornum) {
                delete[] this->StatSend, this->StatRecv, this->ReqSend, this->ReqRecv;
            }
        }

        template<class F>
        void SetBoundary(int *_bctype, F _func);
        template<class F>
        void SetBoundary(F _func) {
            D2Q9<T>::SetBoundary(this->bctype, _func);
        }
        
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
        static int IndexF(int _idx, int _c) {
            return D2Q9<T>::nc*_idx + _c;
        }

        void Swap();
        void BoundaryCondition();
        void iBoundaryCondition();
        void SmoothCorner();
        
        void SetNeighborId(int _leftId, int _rightId, int _bottomId, int _topId, int _leftbottomId, int _lefttopId, int _rightbottomId, int _righttopId);
        void Synchronize();
        
        const int nx, ny, nxy, nbc, offsetxmin, offsetxmax, offsetymin, offsetymax, mpiid;
        static const int nc = 9, nd = 2, cx[nc], cy[nc];
        static const T ei[nc];
        T *f, *fnext, *fsend, *frecv;

private:
        int *bctype;
        int neighborid[8] = {
            -1, -1, -1, -1, -1, -1, -1, -1
        };  //  Left, Right, Bottom, Top, LeftBottom, LeftTop, RightBottom, RightTop
        int neighbornum;
        MPI_Status *StatSend, *StatRecv;
        MPI_Request *ReqSend, *ReqRecv;
    };

    template<class T>const int D2Q9<T>::cx[D2Q9<T>::nc] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    template<class T>const int D2Q9<T>::cy[D2Q9<T>::nc] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    template<class T>const T D2Q9<T>::ei[D2Q9<T>::nc] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };

    template<class T>
    template<class F>
    void D2Q9<T>::SetBoundary(int *_bctype, F _func) {
        for (int j = 0; j < this->ny; ++j) {
            _bctype[j + this->offsetxmin] = _func(0, j);
            _bctype[j + this->offsetxmax] = _func(this->nx - 1, j);
        }
        for (int i = 0; i < this->nx; ++i) {
            _bctype[i + this->offsetymin] = _func(i, 0);
            _bctype[i + this->offsetymax] = _func(i, this->ny - 1);
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
            if (this->neighborid[0] == -1 && this->neighborid[2] == -1) {
                idx = this->Index(0, 0);
                idxx = this->Index(1, 0);
                idxy = this->Index(0, 1);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);
            }
            
            //  Corner at xmin, ymax
            if (this->neighborid[0] == -1 && this->neighborid[3] == -1) {
               idx = this->Index(0, this->ny - 1);
                idxx = this->Index(1, this->ny - 1);
                idxy = this->Index(0,this->ny - 2);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);    
            }
            
            //  Corner at xmax, ymin
            if (this->neighborid[1] == -1 && this->neighborid[2] == -1) {
                idx = this->Index(this->nx - 1, 0);
                idxx = this->Index(this->nx - 2, 0);
                idxy = this->Index(this->nx - 1, 1);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  
            }
            
            //  Corner at xmax, ymax
            if (this->neighborid[1] == -1 && this->neighborid[3] == -1) {
                idx = this->Index(this->nx - 1, this->ny - 1);
                idxx = this->Index(this->nx - 2, this->ny - 1);
                idxy = this->Index(this->nx - 1, this->ny - 2);
                this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);  
            }
        }
    }

    template<class T>
    void D2Q9<T>::SetNeighborId(int _leftId, int _rightId, int _bottomId, int _topId, int _leftbottomId, int _lefttopId, int _rightbottomId, int _righttopId) {
        if (this->neighbornum) {
            delete[] this->StatSend, this->StatRecv, this->ReqSend, this->ReqRecv;
        }
        
        this->neighborid[0] = _leftId;
        this->neighborid[1] = _rightId;
        this->neighborid[2] = _bottomId;
        this->neighborid[3] = _topId;
        this->neighborid[4] = _leftbottomId;
        this->neighborid[5] = _lefttopId;
        this->neighborid[6] = _rightbottomId;
        this->neighborid[7] = _righttopId;

        this->neighbornum = 0;
        for (int i = 0; i < 8; ++i) {
            if (this->neighborid[i] != -1) {
                this->neighbornum += 1;
            }
        }

        if (this->neighbornum) {
            this->StatSend = new MPI_Status[this->neighbornum];
            this->StatRecv = new MPI_Status[this->neighbornum];
            this->ReqSend = new MPI_Request[this->neighbornum];
            this->ReqRecv = new MPI_Request[this->neighbornum];
        } else {
            this->StatSend = nullptr;
            this->StatRecv = nullptr;
            this->ReqSend = nullptr;
            this->ReqRecv = nullptr;
        }
    }

    template<class T>
    void D2Q9<T>::Synchronize() {
        int idx, idxedge, idxcorner, neib = 0;

        //  Copy from f to fedge along xmin
        if (this->neighborid[0] != -1) {
            for (int j = 0; j < this->ny; ++j) {
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmin;
                this->fsend[D2Q9<T>::IndexF(idxedge, 3)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 6)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 7)] = this->f[D2Q9<T>::IndexF(idx, 7)]; 
            }
        }
        //  Copy from f to fedge along xmax
        if (this->neighborid[1] != -1) {
            for (int j = 0; j < this->ny; ++j) {
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmax;
                this->fsend[D2Q9<T>::IndexF(idxedge, 1)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 5)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 8)] = this->f[D2Q9<T>::IndexF(idx, 8)]; 
            }
        }
        //  Copy from f to fedge along ymin
        if (this->neighborid[2] != -1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymin;
                this->fsend[D2Q9<T>::IndexF(idxedge, 4)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 7)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 8)] = this->f[D2Q9<T>::IndexF(idx, 8)]; 
            }
        }
        //  Copy from f to fedge along ymax
        if (this->neighborid[3] != -1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymax;
                this->fsend[D2Q9<T>::IndexF(idxedge, 2)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 5)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend[D2Q9<T>::IndexF(idxedge, 6)] = this->f[D2Q9<T>::IndexF(idx, 6)]; 
            }
        }
        //  Copy from f to fcorner at xmin and ymin
        if (this->neighborid[4] != -1) {
            idx = this->Index(this->nx - 1, this->ny - 1);
            idxcorner = this->nbc + 0;
            this->fsend[D2Q9<T>::IndexF(idxcorner, 7)] = this->f[D2Q9<T>::IndexF(idx, 7)];
        }
        //  Copy from f to fcorner at xmin and ymax
        if (this->neighborid[5] != -1) {
            idx = this->Index(this->nx - 1, 0);
            idxcorner = this->nbc + 1;
            this->fsend[D2Q9<T>::IndexF(idxcorner, 6)] = this->f[D2Q9<T>::IndexF(idx, 6)];
        }
        //  Copy from f to fcorner at xmax and ymin
        if (this->neighborid[6] != -1) {
            idx = this->Index(0, this->ny - 1);
            idxcorner = this->nbc + 2;
            this->fsend[D2Q9<T>::IndexF(idxcorner, 8)] = this->f[D2Q9<T>::IndexF(idx, 8)];
        }
        //  Copy from f to fcorner at xmax and ymax
        if (this->neighborid[7] != -1) {
            idx = this->Index(0, 0);
            idxcorner = this->nbc + 3;
            this->fsend[D2Q9<T>::IndexF(idxcorner, 5)] = this->f[D2Q9<T>::IndexF(idx, 5)];
        }
        
        //  Communicate with other PE p164~参照
        if (this->neighborid[0] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetxmin, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[0], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetxmin, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[0], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[1] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetxmax, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[1], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetxmax, 0)], this->ny*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[1], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[2] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetymin, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[2], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&fsend[D2Q9<T>::IndexF(this->offsetymin, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[2], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[3] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->offsetymax, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[3], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->offsetymax, 0)], this->nx*D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[3], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[4] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 0, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[4], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 0, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[4], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[5] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 1, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[5], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 1, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[5], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[6] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 2, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[6], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 2, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[6], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        if (this->neighborid[7] != -1) {
            MPI_Isend(&fsend[D2Q9<T>::IndexF(this->nbc + 3, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[7], 0, MPI_COMM_WORLD, &ReqSend[neib]);
            MPI_Irecv(&frecv[D2Q9<T>::IndexF(this->nbc + 3, 0)], D2Q9<T>::nc, MPI_DOUBLE, this->neighborid[7], 0, MPI_COMM_WORLD, &ReqRecv[neib]);
            neib++;
        }
        MPI_Waitall(this->neighbornum, this->ReqSend, this->StatSend);
        MPI_Waitall(this->neighbornum, this->ReqRecv, this->StatRecv);

        //  Copy to f from fedge along xmin
        if (this->neighborid[0] != -1) {
            for (int j = 0; j < this->ny; ++j) {
                idx = this->Index(0, j);
                idxedge = j + this->offsetxmin;
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->frecv[D2Q9<T>::IndexF(idxedge, 1)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[D2Q9<T>::IndexF(idxedge, 5)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[D2Q9<T>::IndexF(idxedge, 8)]; 
            }
        }
        //  Copy to f from fedge along xmax
        if (this->neighborid[1] != -1) {
            for (int j = 0; j < this->ny; ++j) {
                idx = this->Index(this->nx - 1, j);
                idxedge = j + this->offsetxmax;
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->frecv[D2Q9<T>::IndexF(idxedge, 3)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[D2Q9<T>::IndexF(idxedge, 6)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[D2Q9<T>::IndexF(idxedge, 7)];
            }
        }
        //  Copy to f from fedge along ymin
        if (this->neighborid[2] != -1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, 0);
                idxedge = i + this->offsetymin;
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->frecv[D2Q9<T>::IndexF(idxedge, 2)];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[D2Q9<T>::IndexF(idxedge, 5)];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[D2Q9<T>::IndexF(idxedge, 6)];
            }
        }
        //  Copy to f from fedge along ymax
        if (this->neighborid[3] != -1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, this->ny - 1);
                idxedge = i + this->offsetymax;
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->frecv[D2Q9<T>::IndexF(idxedge, 4)];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[D2Q9<T>::IndexF(idxedge, 7)];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[D2Q9<T>::IndexF(idxedge, 8)];
            }
        }
        //  Copy to f from fcorner at xmin and ymin
        if (this->neighborid[4] != -1) {
            idx = this->Index(this->nx - 1, this->ny - 1);
            idxcorner = this->nbc + 0;
            this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv[D2Q9<T>::IndexF(idxcorner, 5)];
        }
        //  Copy to f from fcorner at xmin and ymax
        if (this->neighborid[5] != -1) {
            idx = this->Index(this->nx - 1, 0);
            idxcorner = this->nbc + 1;
            this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv[D2Q9<T>::IndexF(idxcorner, 8)];
        }
        //  Copy to f from fcorner at xmax and ymin
        if (this->neighborid[6] != -1) {
            idx = this->Index(0, this->ny - 1);
            idxcorner = this->nbc + 2;
            this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv[D2Q9<T>::IndexF(idxcorner, 6)];
        }
        //  Copy to f from fcorner at xmax and ymax
        if (this->neighborid[7] != -1) {
            idx = this->Index(0, 0);
            idxcorner = this->nbc + 3;
            this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv[D2Q9<T>::IndexF(idxcorner, 7)];
        }
    }
}