//*****************************************************************************
//  Title       :   src/particle/d3q15.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/29
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
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
    class D3Q15 {
public:
        D3Q15() = delete;
        D3Q15(int _lx, int _ly, int _lz, int _PEid = 0, int _mx = 1, int _my = 1, int _mz = 1) :
            lx(_lx), ly(_ly), lz(_lz), PEid(_PEid), mx(_mx), my(_my), mz(_mz), 
            PEx(this->PEid%this->mx), PEy((this->PEid/this->mx)%this->my), PEz(this->PEid/(this->mx*this->my)),
            nx((this->lx + this->PEx)/this->mx), ny((this->ly + this->PEy)/this->my), nz((this->lz + this->PEz)/this->mz),
            nxyz(this->nx*this->ny*this->nz), nbc(2*(this->nx*this->ny + this->ny*this->nz + this->nz*this->nx)),
            offsetxmin(0), offsetxmax(this->ny*this->nz), 
            offsetymin(2*this->ny*this->nz), offsetymax(2*this->ny*this->nz + this->nz*this->nx),
            offsetzmin(2*(this->ny*this->nz + this->nz*this->nx)), offsetzmax(2*(this->ny*this->nz + this->nz*this->nx) + this->nx*this->ny),
            offsetx(this->mx - this->PEx > this->lx%this->mx ? this->PEx*this->nx : this->lx - (this->mx - this->PEx)*this->nx),
            offsety(this->my - this->PEy > this->ly%this->my ? this->PEy*this->ny : this->ly - (this->my - this->PEy)*this->ny),
            offsetz(this->mz - this->PEz > this->lz%this->mz ? this->PEz*this->nz : this->lz - (this->mz - this->PEz)*this->nz)
        {
            assert(0 < _lx && 0 < _ly && 0 < _lz && 0 <= _PEid && 0 < _mx && 0 < _my && 0 < _mz);

            this->f = new T[this->nxyz*D3Q15<T>::nc];
            this->fnext = new T[this->nxyz*D3Q15<T>::nc];
            this->fsend = new T[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 8];
            this->frecv = new T[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 8];
            this->bctype = new int[this->nbc];
            for (int idx = 0; idx < this->nbc; ++idx) {
                this->bctype[idx] = 0;
            }
        }
        D3Q15(const D3Q15<T>& _p) = delete;
        ~D3Q15() {
            delete[] this->f, this->fnext, this->fsend, this->frecv, this->bctype;
        }

        template<class F>
        void SetBoundary(int *_bctype, F _func);
        template<class F>
        void SetBoundary(F _func) {
            this->SetBoundary(this->bctype, _func);
        }
        T GetBoundary(int *_bctype, int _i, int _j, int _k) const;
        T GetBoundary(int _i, int _j, int _k) const {
            return this->GetBoundary(this->bctype, _i, _j, _k);
        }

        int Index(int _i, int _j, int _k) const {
            int i = _i == -1 ? this->nx - 1 : (_i == this->nx ? 0 : _i);
            int j = _j == -1 ? this->ny - 1 : (_j == this->ny ? 0 : _j);
            int k = _k == -1 ? this->nz - 1 : (_k == this->nz ? 0 : _k);
            return i + this->nx*j + this->nx*this->ny*k;
        }
        static int IndexF(int _idx, int _c) {
            return D3Q15<T>::nc*_idx + _c;
        }
        int IndexPE(int _i, int _j, int _k) const {
            int i = _i == -1 ? this->mx - 1 : (_i == this->mx ? 0 : _i);
            int j = _j == -1 ? this->my - 1 : (_j == this->my ? 0 : _j);
            int k = _k == -1 ? this->mz - 1 : (_k == this->mz ? 0 : _k);
            return i + this->mx*j + this->mx*this->my*k;
        }
        int IndexBCx(int _j, int _k) const {
            return _j + this->ny*_k;
        }
        int IndexBCy(int _k, int _i) const {
            return _k + this->nz*_i;
        }
        int IndexBCz(int _i, int _j) const {
            return _i + this->nx*_j;
        }

        void Swap();
        void BoundaryCondition();
        void iBoundaryCondition();
        void SmoothCorner();

        void Synchronize();
        void iSynchronize();

        const int lx, ly, lz, PEid, mx, my, mz, PEx, PEy, PEz, nx, ny, nz, nxyz, nbc, offsetxmin, offsetxmax, offsetymin, offsetymax, offsetzmin, offsetzmax, offsetx, offsety, offsetz;
        static const int nc = 15, nd = 3, cx[nc], cy[nc], cz[nc];
        static const T ei[nc];
        T *f, *fnext, *fsend, *frecv;

private:
        int *bctype;
#ifdef _USE_MPI_DEFINES
        MPI_Status status[52];
        MPI_Request request[52];
#endif
    };

    template<class T>const int D3Q15<T>::cx[D3Q15<T>::nc] = { 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, -1, -1 };
    template<class T>const int D3Q15<T>::cy[D3Q15<T>::nc] = { 0, 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, -1, -1, 1, -1 };
    template<class T>const int D3Q15<T>::cz[D3Q15<T>::nc] = { 0, 0, 0, 1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1 };
    template<class T>const T D3Q15<T>::ei[D3Q15<T>::nc] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0 };

    template<class T>
    template<class F>
    void D3Q15<T>::SetBoundary(int *_bctype, F _func) {
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    _bctype[this->IndexBCx(j, k) + this->offsetxmin] = _func(0 + this->offsetx, j + this->offsety, k + this->offsetz);
                }
            }
        }
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    _bctype[this->IndexBCx(j, k) + this->offsetxmax] = _func((this->nx - 1) + this->offsetx, j + this->offsety, k + this->offsetz);
                }   
            }
        }
        if (this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    _bctype[this->IndexBCy(k, i) + this->offsetymin] = _func(i + this->offsetx, 0 + this->offsety, k + this->offsetz);
                }
            }
        }
        if (this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    _bctype[this->IndexBCy(k, i) + this->offsetymax] = _func(i + this->offsetx, (this->ny - 1) + this->offsety, k + this->offsetz);
                }
            }
        }
        if (this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    _bctype[this->IndexBCz(i, j) + this->offsetzmin] = _func(i + this->offsetx, j + this->offsety, 0 + this->offsetz);
                }
            }
        }
        if (this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    _bctype[this->IndexBCz(i, j) + this->offsetzmax] = _func(i + this->offsetx, j + this->offsety, (this->nz - 1) + this->offsetz);
                }
            }
        }
    }

    template<class T>
    T D3Q15<T>::GetBoundary(int *_bctype, int _i, int _j, int _k) const {
        if (_i == 0) {
            return T(_bctype[this->IndexBCx(_j, _k) + this->offsetxmin]);
        } else if (_i == this->nx - 1) {
            return T(_bctype[this->IndexBCx(_j, _k) + this->offsetxmax]);
        } else if (_j == 0) {
            return T(_bctype[this->IndexBCy(_k, _i) + this->offsetymin]);
        } else if (_j == this->ny - 1) {
            return T(_bctype[this->IndexBCy(_k, _i) + this->offsetymax]);
        } else if (_k == 0) {
            return T(_bctype[this->IndexBCz(_i, _j) + this->offsetzmin]);
        } else if (_k == this->nz - 1) {
            return T(_bctype[this->IndexBCz(_i, _j) + this->offsetzmax]);
        } else {
            return T();
        }
    }

    template<class T>
    void D3Q15<T>::Swap() {
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;
    }

    template<class T>
    void D3Q15<T>::BoundaryCondition() {
        int idx, idxbc;
        for (int j = 0; j < this->ny; ++j) {
            for (int k = 0; k < this->nz; ++k) {
                //  On xmin
                idx = this->Index(0, j, k);
                idxbc = this->IndexBCx(j, k) + this->offsetxmin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }

                //  On xmax
                idx = this->Index(this->nx - 1, j, k);
                idxbc = this->IndexBCx(j, k) + this->offsetxmax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                }
            }
        }
        for (int k = 0; k < this->nz; ++k) {
            for (int i = 0; i < this->nx; ++i) {
                //  On ymin
                idx = this->Index(i, 0, k);
                idxbc = this->IndexBCy(k, i) + this->offsetymin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }

                //  On ymax
                idx = this->Index(i, this->ny - 1, k);
                idxbc = this->IndexBCy(k, i) + this->offsetymax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                }
            }
        }
        for (int i = 0; i < this->nx; ++i) {
            for (int j = 0; j < this->ny; ++j) {
                //  On zmin
                idx = this->Index(i, j, 0);
                idxbc = this->IndexBCz(i, j) + this->offsetzmin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }

                //  On zmax
                idx = this->Index(i, j, this->nz - 1);
                idxbc = this->IndexBCz(i, j) + this->offsetzmax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                }
            }
        }
    }

    template<class T>
    void D3Q15<T>::iBoundaryCondition() {
        int idx, idxbc;
        for (int j = 0; j < this->ny; ++j) {
            for (int k = 0; k < this->nz; ++k) {
                //  On xmin
                idx = this->Index(0, j, k);
                idxbc = this->IndexBCx(j, k) + this->offsetxmin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                }

                //  On xmax
                idx = this->Index(this->nx - 1, j, k);
                idxbc = this->IndexBCx(j, k) + this->offsetxmax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }
            }
        }
        for (int k = 0; k < this->nz; ++k) {
            for (int i = 0; i < this->nx; ++i) {
                //  On ymin
                idx = this->Index(i, 0, k);
                idxbc = this->IndexBCy(k, i) + this->offsetymin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                }

                //  On ymax
                idx = this->Index(i, this->ny - 1, k);
                idxbc = this->IndexBCy(k, i) + this->offsetymax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }
            }
        }
        for (int i = 0; i < this->nx; ++i) {
            for (int j = 0; j < this->ny; ++j) {
                //  On zmin
                idx = this->Index(i, j, 0);
                idxbc = this->IndexBCz(i, j) + this->offsetzmin;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                }

                //  On zmax
                idx = this->Index(i, j, this->nz - 1);
                idxbc = this->IndexBCz(i, j) + this->offsetzmax;
                if (this->bctype[idxbc] == BARRIER) {
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                } else if (this->bctype[idxbc] == MIRROR) {
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                }
            }
        }
    }

    template<class T>
    void D3Q15<T>::SmoothCorner() {
        int idx, idxx, idxy, idxz;

        //  Along line ymin and zmin
        if (this->PEy == 0 && this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, 0, 0);
                idxy = this->Index(i, 1, 0);
                idxz = this->Index(i, 0, 1);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymax and zmin
        if (this->PEy == this->my - 1 && this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, this->ny - 1, 0);
                idxy = this->Index(i, this->ny - 2, 0);
                idxz = this->Index(i, this->ny - 1, 1);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymax and zmax
        if (this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                idxy = this->Index(i, this->ny - 2, this->nz - 1);
                idxz = this->Index(i, this->ny - 1, this->nz - 2);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymin and zmax
        if (this->PEy == 0 && this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                idx = this->Index(i, 0, this->nz - 1);
                idxy = this->Index(i, 1, this->nz - 1);
                idxz = this->Index(i, 0, this->nz - 2);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }

        //  Along line zmin and xmin
        if (this->PEz == 0 && this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                idx = this->Index(0, j, 0);
                idxz = this->Index(1, j, 0);
                idxx = this->Index(0, j, 1);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmax and xmin
        if (this->PEz == this->mz - 1 && this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {                
                idx = this->Index(0, j, this->nz - 1);
                idxz = this->Index(0, j, this->nz - 2);
                idxx = this->Index(1, j, this->nz - 1);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmax and xmax
        if (this->PEz == this->mz - 1 && this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {                                        
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                idxz = this->Index(this->nx - 1, j, this->nz - 2);
                idxx = this->Index(this->nx - 2, j, this->nz - 1);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmin and xmax
        if (this->PEz == 0 && this->PEy == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {                                        
                idx = this->Index(this->nx - 1, j, 0);
                idxz = this->Index(this->nx - 1, j, 1);
                idxx = this->Index(this->nx - 2, j, 0);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        
        //  Along line xmin and ymin
        if (this->PEx == 0 && this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {
                idx = this->Index(0, 0, k);
                idxx = this->Index(1, 0, k);
                idxy = this->Index(0, 1, k);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmax and ymin
        if (this->PEx == this->mx - 1 && this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {               
                idx = this->Index(this->nx - 1, 0, k);
                idxx = this->Index(this->nx - 2, 0, k);
                idxy = this->Index(this->nx - 1, 1, k);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmax and ymax
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {               
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                idxx = this->Index(this->nx - 2, this->ny - 1, k);
                idxy = this->Index(this->nx - 1, this->ny - 2, k);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmin and ymax
        if (this->PEx == 0 && this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {               
                idx = this->Index(0, this->ny - 1, k);
                idxx = this->Index(1, this->ny - 1, k);
                idxy = this->Index(0, this->ny - 2, k);
                for (int c = 0; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        
        //  Corner at xmin, ymin and zmin
        if (this->PEx == 0 && this->PEy == 0 && this->PEz == 0) {
            idx = this->Index(0, 0, 0);
            idxx = this->Index(1, 0, 0);
            idxy = this->Index(0, 1, 0);
            idxz = this->Index(0, 0, 1);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymin and zmin
        if (this->PEx == this->mx - 1 && this->PEy == 0 && this->PEz == 0) {
            idx = this->Index(this->nx - 1, 0, 0);
            idxx = this->Index(this->nx - 2, 0, 0);
            idxy = this->Index(this->nx - 1, 1, 0);
            idxz = this->Index(this->nx - 1, 0, 1);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymax and zmin
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1 && this->PEz == 0) {
            idx = this->Index(this->nx - 1, this->ny - 1, 0);
            idxx = this->Index(this->nx - 2, this->ny - 1, 0);
            idxy = this->Index(this->nx - 1, this->ny - 2, 0);
            idxz = this->Index(this->nx - 1, this->ny - 1, 1);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymax and zmin
        if (this->PEx == 0 && this->PEy == this->my - 1 && this->PEz == 0) {
            idx = this->Index(0, this->ny - 1, 0);
            idxx = this->Index(1, this->ny - 1, 0);
            idxy = this->Index(0, this->ny - 2, 0);
            idxz = this->Index(0, this->ny - 1, 1);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymin and zmax
        if (this->PEx == 0 && this->PEy == 0 && this->PEz == this->mz - 1) {
            idx = this->Index(0, 0, this->nz - 1);
            idxx = this->Index(1, 0, this->nz - 1);
            idxy = this->Index(0, 1, this->nz - 1);
            idxz = this->Index(0, 0, this->nz - 2);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymin and zmax
        if (this->PEx == this->mx - 1 && this->PEy == 0 && this->PEz == this->mz - 1) {
            idx = this->Index(this->nx - 1, 0, this->nz - 1);
            idxx = this->Index(this->nx - 2, 0, this->nz - 1);
            idxy = this->Index(this->nx - 1, 1, this->nz - 1);
            idxz = this->Index(this->nx - 1, 0, this->nz - 2);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymax and zmax
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            idx = this->Index(this->nx - 1, this->ny - 1, this->nz - 1);
            idxx = this->Index(this->nx - 2, this->ny - 1, this->nz - 1);
            idxy = this->Index(this->nx - 1, this->ny - 2, this->nz - 1);
            idxz = this->Index(this->nx - 1, this->ny - 1, this->nz - 2);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymax and zmax
        if (this->PEx == 0 && this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            idx = this->Index(0, this->ny - 1, this->nz - 1);
            idxx = this->Index(1, this->ny - 1, this->nz - 1);
            idxy = this->Index(0, this->ny - 2, this->nz - 1);
            idxz = this->Index(0, this->ny - 1, this->nz - 2);
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
    }

    template<class T>
    void D3Q15<T>::Synchronize() {
#ifdef _USE_MPI_DEFINES
        int idx, idxface, idxedge, neib = 0;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on xmax
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 12)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on ymax
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j) + this->offsetzmin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];

                    //  Face on zmax
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j) + this->offsetzmax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                idxedge = i + 0*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on ymin and zmax
                idx = this->Index(i, this->ny - 1, 0);
                idxedge = i + 1*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on ymax and zmin
                idx = this->Index(i, 0, this->nz - 1);
                idxedge = i + 2*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on ymax and zmax
                idx = this->Index(i, 0, 0);
                idxedge = i + 3*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                idxedge = j + 4*this->nx + 0*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on zmin and xmax
                idx = this->Index(0, j, this->nz - 1);
                idxedge = j + 4*this->nx + 1*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on zmax and xmin
                idx = this->Index(this->nx - 1, j, 0);
                idxedge = j + 4*this->nx + 2*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on zmax and xmax
                idx = this->Index(0, j, 0);
                idxedge = j + 4*this->nx + 3*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 0*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on xmin and ymax
                idx = this->Index(this->nx - 1, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 1*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on xmax and ymin
                idx = this->Index(0, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 2*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on xmax and ymax
                idx = this->Index(0, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 3*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)]; //  Corner at xmin, ymin and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)];            //  Corner at xmax, ymin and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)];            //  Corner at xmin, ymax and zmin 
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)];                       //  Corner at xmax, ymax and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)];            //  Corner at xmin, ymin and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)];                        //  Corner at xmax, ymin and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)];                        //  Corner at xmin, ymax and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)];                                   //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        if (this->mx != 1) {
            //  To xmin
            MPI_Isend(&fsend[this->offsetxmin*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetxmax*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  To xmax
            MPI_Isend(&fsend[this->offsetxmax*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetxmin*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  To ymin
            MPI_Isend(&fsend[this->offsetymin*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetymax*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  To ymax
            MPI_Isend(&fsend[this->offsetymax*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetymin*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1) {
            //  To zmin
            MPI_Isend(&fsend[this->offsetzmin*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetzmax*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  To zmax
            MPI_Isend(&fsend[this->offsetzmax*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetzmin*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1 || this->mz != 1) {
            //  To ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 0*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 3*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 1*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 2*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 2*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 1*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmax
            MPI_Isend(&fsend[this->nbc*5 + 3*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 0*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1 || this->mx != 1) {
            //  To zmin and xmin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 0*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 3*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To zmin and xmax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 1*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 2*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
                 
            //  To zmax and xmin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 2*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 1*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To zmax and xmax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 3*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 0*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  To xmin and ymin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 0*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 3*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To xmin and ymax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 1*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 2*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 2*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 1*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 3*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 0*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmax
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmin;
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 4];

                    //  Face on xmax
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmax;
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymin;
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 4];

                    //  Face on ymax
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymax;
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j) + this->offsetzmin;
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];

                    //  Face on zmax
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j) + this->offsetzmax;
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, 0, 0);
                idxedge = i + 0*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymin and zmax
                idx = this->Index(i, 0, this->nz - 1);
                idxedge = i + 1*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymax and zmin
                idx = this->Index(i, this->ny - 1, 0);
                idxedge = i + 2*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymax and zmax
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                idxedge = i + 3*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(0, j, 0);
                idxedge = j + 4*this->nx + 0*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmin and xmax
                idx = this->Index(this->nx - 1, j, 0);
                idxedge = j + 4*this->nx + 1*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmax and xmin
                idx = this->Index(0, j, this->nz - 1);
                idxedge = j + 4*this->nx + 2*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmax and xmax
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                idxedge = j + 4*this->nx + 3*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(0, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 0*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmin and ymax
                idx = this->Index(0, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 1*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmax and ymin
                idx = this->Index(this->nx - 1, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 2*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmax and ymax
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 3*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0];                                   //  Corner at xmin, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1];                        //  Corner at xmax, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2];                        //  Corner at xmin, ymax and zmin 
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3];            //  Corner at xmax, ymax and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4];                       //  Corner at xmin, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5];            //  Corner at xmax, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6];            //  Corner at xmin, ymax and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7]; //  Corner at xmax, ymax and zmax
        }
#endif
    }

    template<class T>
    void D3Q15<T>::iSynchronize() {
#ifdef _USE_MPI_DEFINES
        int idx, idxface, idxedge, neib = 0;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 12)];

                    //  Face on xmax
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];

                    //  Face on ymax
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j) + this->offsetzmin;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on zmax
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j) + this->offsetzmax;
                    this->fsend[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->fsend[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                idxedge = i + 0*this->nx;
                 this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];

                //  Edge on ymin and zmax
                idx = this->Index(i, this->ny - 1, 0);
                idxedge = i + 1*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on ymax and zmin
                idx = this->Index(i, 0, this->nz - 1);
                idxedge = i + 2*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on ymax and zmax
                idx = this->Index(i, 0, 0);
                idxedge = i + 3*this->nx;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                idxedge = j + 4*this->nx + 0*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];

                //  Edge on zmin and xmax
                idx = this->Index(0, j, this->nz - 1);
                idxedge = j + 4*this->nx + 1*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on zmax and xmin
                idx = this->Index(this->nx - 1, j, 0);
                idxedge = j + 4*this->nx + 2*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on zmax and xmax
                idx = this->Index(0, j, 0);
                idxedge = j + 4*this->nx + 3*this->ny;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 0*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];

                //  Edge on xmin and ymax
                idx = this->Index(this->nx - 1, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 1*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];
                
                //  Edge on xmax and ymin
                idx = this->Index(0, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 2*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on xmax and ymax
                idx = this->Index(0, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 3*this->nz;
                this->fsend[this->nbc*5 + idxedge*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend[this->nbc*5 + idxedge*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 7)];  //  Corner at xmin, ymin and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 8)];             //  Corner at xmax, ymin and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 9)];             //  Corner at xmin, ymax and zmin 
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 14)];                       //  Corner at xmax, ymax and zmin
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 10)];            //  Corner at xmin, ymin and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 13)];                       //  Corner at xmax, ymin and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 12)];                       //  Corner at xmin, ymax and zmax
            this->fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 11)];                                  //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        if (this->mx != 1) {
            //  To xmin
            MPI_Isend(&fsend[this->offsetxmin*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetxmax*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  To xmax
            MPI_Isend(&fsend[this->offsetxmax*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetxmin*5], this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  To ymin
            MPI_Isend(&fsend[this->offsetymin*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetymax*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  To ymax
            MPI_Isend(&fsend[this->offsetymax*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetymin*5], this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1) {
            //  To zmin
            MPI_Isend(&fsend[this->offsetzmin*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetzmax*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  To zmax
            MPI_Isend(&fsend[this->offsetzmax*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->offsetzmin*5], this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1 || this->mz != 1) {
            //  To ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 0*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 0*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 1*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 1*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 2*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 2*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmax
            MPI_Isend(&fsend[this->nbc*5 + 3*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 3*this->nx*2], this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1 || this->mx != 1) {
            //  To zmin and xmin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 0*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 0*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To zmin and xmax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 1*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 1*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
                 
            //  To zmax and xmin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 2*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 2*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To zmax and xmax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 3*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 3*this->ny*2], this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  To xmin and ymin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 0*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 0*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To xmin and ymax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 1*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 1*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymin
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 2*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 2*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymax
            MPI_Isend(&fsend[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 3*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*this->nx*2 + 4*this->ny*2 + 3*this->nz*2], this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmin
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmax
            MPI_Isend(&fsend[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmax
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmin;
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 4];

                    //  Face on xmax
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k) + this->offsetxmax;
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymin;
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 4];

                    //  Face on ymax
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i) + this->offsetymax;
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j) + this->offsetzmin;
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[idxface*5 + 4];

                    //  Face on zmax
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j) + this->offsetzmax;
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->frecv[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, 0, 0);
                idxedge = i + 0*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymin and zmax
                idx = this->Index(i, 0, this->nz - 1);
                idxedge = i + 1*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymax and zmin
                idx = this->Index(i, this->ny - 1, 0);
                idxedge = i + 2*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on ymax and zmax
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                idxedge = i + 3*this->nx;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(0, j, 0);
                idxedge = j + 4*this->nx + 0*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmin and xmax
                idx = this->Index(this->nx - 1, j, 0);
                idxedge = j + 4*this->nx + 1*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmax and xmin
                idx = this->Index(0, j, this->nz - 1);
                idxedge = j + 4*this->nx + 2*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on zmax and xmax
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                idxedge = j + 4*this->nx + 3*this->ny;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(0, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 0*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmin and ymax
                idx = this->Index(0, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 1*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmax and ymin
                idx = this->Index(this->nx - 1, 0, k);
                idxedge = k + 4*this->nx + 4*this->ny + 2*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv[this->nbc*5 + idxedge*2 + 1];

                //  Edge on xmax and ymax
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                idxedge = k + 4*this->nx + 4*this->ny + 3*this->nz;
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv[this->nbc*5 + idxedge*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv[this->nbc*5 + idxedge*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 0];                                   //  Corner at xmin, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 1];                        //  Corner at xmax, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 2];                        //  Corner at xmin, ymax and zmin 
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 3];            //  Corner at xmax, ymax and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 4];                       //  Corner at xmin, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 5];            //  Corner at xmax, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 6];            //  Corner at xmin, ymax and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)] = this->frecv[this->nbc*5 + 4*(this->nx + this->ny + this->nz)*2 + 7]; //  Corner at xmax, ymax and zmax
        }
#endif
    }
}