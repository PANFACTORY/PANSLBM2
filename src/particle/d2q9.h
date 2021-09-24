//*****************************************************************************
//  Title       :   src/particle/d2q9.h
//  Author      :   Tanabe Yuta
//  Date        :   2021/8/14
//  Copyright   :   (C)2021 TanabeYuta
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
    class D2Q9 {
public:
        D2Q9() = delete;
        D2Q9(int _lx, int _ly, int _PEid = 0, int _mx = 1, int _my = 1) :
            lx(_lx), ly(_ly), lz(1), PEid(_PEid), mx(_mx), my(_my), mz(1), 
            PEx(this->PEid%this->mx), PEy(this->PEid/this->mx), PEz(0),
            nx((this->lx + this->PEx)/this->mx), ny((this->ly + this->PEy)/this->my), nz(1),
            nxyz(this->nx*this->ny),
            offsetx(this->mx - this->PEx > this->lx%this->mx ? this->PEx*this->nx : this->lx - (this->mx - this->PEx)*this->nx),
            offsety(this->my - this->PEy > this->ly%this->my ? this->PEy*this->ny : this->ly - (this->my - this->PEy)*this->ny),
            offsetz(0)
        {
            assert(0 < _lx && 0 < _ly && 0 <= _PEid && 0 < _mx && 0 < _my);

            this->f = new T[this->nxyz*D2Q9<T>::nc];
            this->fnext = new T[this->nxyz*D2Q9<T>::nc];
            this->fsend_xmin = new T[this->ny*3];
            this->fsend_xmax = new T[this->ny*3];
            this->fsend_ymin = new T[this->nx*3];
            this->fsend_ymax = new T[this->nx*3];
            this->frecv_xmin = new T[this->ny*3];
            this->frecv_xmax = new T[this->ny*3];
            this->frecv_ymin = new T[this->nx*3];
            this->frecv_ymax = new T[this->nx*3];
        }
        D2Q9(const D2Q9<T>& _p) = delete;
        ~D2Q9() {
            delete[] this->f, this->fnext;
            delete[] this->fsend_xmin, this->fsend_xmax, this->fsend_ymin, this->fsend_ymax;
            delete[] this->frecv_xmin, this->frecv_xmax, this->frecv_ymin, this->frecv_ymax;
        }
        
        int Index(int _i, int _j) const {
            int i = _i == -1 ? this->nx - 1 : (_i == this->nx ? 0 : _i);
            int j = _j == -1 ? this->ny - 1 : (_j == this->ny ? 0 : _j);
            return i + this->nx*j;
        }
        int Index(int _i, int _j, int _k) const {
            return this->Index(_i, _j);
        }
        static int IndexF(int _idx, int _c) {
            return D2Q9<T>::nc*_idx + _c;
        }
        int IndexPE(int _i, int _j) const {
            int i = _i == -1 ? this->mx - 1 : (_i == this->mx ? 0 : _i);
            int j = _j == -1 ? this->my - 1 : (_j == this->my ? 0 : _j);
            return i + this->mx*j;
        }
        int IndexPE(int _i, int _j, int _k) const {
            return this->IndexPE(_i, _j);
        }

        void Stream();
        void Swap();
        template<class Ff>
        void BoundaryCondition(Ff _bctype);
        template<class Ff>
        void iBoundaryCondition(Ff _bctype);
        void SmoothCorner();
        
        void Synchronize();
        void iSynchronize();
        
        const int lx, ly, lz, PEid, mx, my, mz, PEx, PEy, PEz, nx, ny, nz, nxyz, offsetx, offsety, offsetz;
        static const int nc = 9, nd = 2, cx[nc], cy[nc], cz[nc];
        static const T ei[nc];
        T *f, *fnext;
        
private:
        T *fsend_xmin, *fsend_xmax, *fsend_ymin, *fsend_ymax, *frecv_xmin, *frecv_xmax, *frecv_ymin, *frecv_ymax;
        T fsend_corner[4], frecv_corner[4];
#ifdef _USE_MPI_DEFINES
        MPI_Status status[16];
        MPI_Request request[16];
#endif
    };

    template<class T>const int D2Q9<T>::cx[D2Q9<T>::nc] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
    template<class T>const int D2Q9<T>::cy[D2Q9<T>::nc] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
    template<class T>const int D2Q9<T>::cz[D2Q9<T>::nc] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    template<class T>const T D2Q9<T>::ei[D2Q9<T>::nc] = { 4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0 };

    template<class T>
    void D2Q9<T>::Stream() {
        //  Stream
#pragma omp parallel for
        for (int j = 0; j < this->ny; ++j) {
            for (int i = 0; i < this->nx; ++i) {
                for (int c = 0; c < D2Q9<T>::nc; ++c) {
                    int idx = this->Index(i, j), idxstream = this->Index(i - D2Q9<T>::cx[c], j - D2Q9<T>::cy[c]);
                    this->fnext[D2Q9<T>::IndexF(idx, c)] = this->f[D2Q9<T>::IndexF(idxstream, c)];
                }
            }
        }

        //  Swap
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;
    }

    template<class T>
    void D2Q9<T>::Swap() {
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;
    }

    template<class T>
    template<class Ff>
    void D2Q9<T>::BoundaryCondition(Ff _bctype) {
        //  On xmin
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                if (_bctype(0 + this->offsetx, j + this->offsety) == BARRIER) {
                    int idx = this->Index(0, j);
                    this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                } else if (_bctype(0 + this->offsetx, j + this->offsety) == MIRROR) {
                    int idx = this->Index(0, j);
                    this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                }
            }
        }
        //  On xmax
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety) == BARRIER) {
                    int idx = this->Index(this->nx - 1, j);
                    this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                } else if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety) == MIRROR) {
                    int idx = this->Index(this->nx - 1, j);
                    this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                }
            }
        }
        //  On ymin
        if (this->PEy == 0) {
            for (int i = 0; i < this->nx; ++i) {
                if (_bctype(i + this->offsetx, 0 + this->offsety) == BARRIER) {
                    int idx = this->Index(i, 0);
                    this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                } else if (_bctype(i + this->offsetx, 0 + this->offsety) == MIRROR) {
                    int idx = this->Index(i, 0);
                    this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                }
            }
        }
        //  On ymax
        if (this->PEy == this->my - 1) {
            for (int i = 0; i < this->nx; ++i) {
                if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety) == BARRIER) {
                    int idx = this->Index(i, this->ny - 1);
                    this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                } else if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety) == MIRROR) {
                    int idx = this->Index(i, this->ny - 1);
                    this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                }
            }
        }
    }

    template<class T>
    template<class Ff>
    void D2Q9<T>::iBoundaryCondition(Ff _bctype) {
        //  On xmin
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                if (_bctype(0 + this->offsetx, j + this->offsety) == BARRIER) {
                    int idx = this->Index(0, j);
                    this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                } else if (_bctype(0 + this->offsetx, j + this->offsety) == MIRROR) {
                    int idx = this->Index(0, j);
                    this->f[D2Q9<T>::IndexF(idx, 3)] = this->f[D2Q9<T>::IndexF(idx, 1)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                }
            }
        }
        //  On xmax
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety) == BARRIER) {
                    int idx = this->Index(this->nx - 1, j);
                    this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                } else if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety) == MIRROR) {
                    int idx = this->Index(this->nx - 1, j);
                    this->f[D2Q9<T>::IndexF(idx, 1)] = this->f[D2Q9<T>::IndexF(idx, 3)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                }
            }
        }
        //  On ymin
        if (this->PEy == 0) {
            for (int i = 0; i < this->nx; ++i) {
                if (_bctype(i + this->offsetx, 0 + this->offsety) == BARRIER) {
                    int idx = this->Index(i, 0);
                    this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                } else if (_bctype(i + this->offsetx, 0 + this->offsety) == MIRROR) {
                    int idx = this->Index(i, 0);
                    this->f[D2Q9<T>::IndexF(idx, 4)] = this->f[D2Q9<T>::IndexF(idx, 2)];
                    this->f[D2Q9<T>::IndexF(idx, 8)] = this->f[D2Q9<T>::IndexF(idx, 5)];
                    this->f[D2Q9<T>::IndexF(idx, 7)] = this->f[D2Q9<T>::IndexF(idx, 6)];
                }
            }
        }
        //  On ymax
        if ( this->PEy == this->my - 1) {
            for (int i = 0; i < this->nx; ++i) {
                if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety) == BARRIER) {
                    int idx = this->Index(i, this->ny - 1);
                    this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                } else if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety) == MIRROR) {
                    int idx = this->Index(i, this->ny - 1);
                    this->f[D2Q9<T>::IndexF(idx, 2)] = this->f[D2Q9<T>::IndexF(idx, 4)];
                    this->f[D2Q9<T>::IndexF(idx, 6)] = this->f[D2Q9<T>::IndexF(idx, 7)];
                    this->f[D2Q9<T>::IndexF(idx, 5)] = this->f[D2Q9<T>::IndexF(idx, 8)];
                }
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
        int idx, neib = 0;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin  
                idx = this->Index(this->nx - 1, j);
                this->fsend_xmin[j*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->fsend_xmin[j*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->fsend_xmin[j*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 7)];

                //  Edge along xmax
                idx = this->Index(0, j);
                this->fsend_xmax[j*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->fsend_xmax[j*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend_xmax[j*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)]; 
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, this->ny - 1);
                this->fsend_ymin[i*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->fsend_ymin[i*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->fsend_ymin[i*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)]; 

                //  Edge along ymax
                idx = this->Index(i, 0);
                this->fsend_ymax[i*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->fsend_ymax[i*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend_ymax[i*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 6)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->fsend_corner[0] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 7)];   //  Corner at xmin and ymin
            this->fsend_corner[1] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 6)];              //  Corner at xmin and ymax
            this->fsend_corner[2] = this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 8)];              //  Corner at xmax and ymin
            this->fsend_corner[3] = this->f[D2Q9<T>::IndexF(this->Index(0, 0), 5)];                         //  Corner at xmax and ymax
        }
        
        //  Communicate with other PE
        if (this->mx != 1) {
            //  Left
            MPI_Isend(this->fsend_xmin, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right
            MPI_Isend(this->fsend_xmax, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  Bottom
            MPI_Isend(this->fsend_ymin, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  Top
            MPI_Isend(this->fsend_ymax, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  Left bottom
            MPI_Isend(&this->fsend_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  Left top
            MPI_Isend(&this->fsend_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right bottom
            MPI_Irecv(&this->frecv_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Isend(&this->fsend_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right top
            MPI_Irecv(&this->frecv_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Isend(&this->fsend_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);      
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin
                idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->frecv_xmin[j*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv_xmin[j*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv_xmin[j*3 + 2]; 
            
                //  Edge along xmax
                idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->frecv_xmax[j*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv_xmax[j*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv_xmax[j*3 + 2];
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->frecv_ymin[i*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv_ymin[i*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv_ymin[i*3 + 2];
            
                //  Edge along ymax
                idx = this->Index(i, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->frecv_ymax[i*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv_ymax[i*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv_ymax[i*3 + 2];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->f[D2Q9<T>::IndexF(this->Index(0, 0), 5)] = this->frecv_corner[0];                         //  Corner at xmin and ymin
            this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 8)] = this->frecv_corner[1];              //  Corner at xmin and ymax
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 6)] = this->frecv_corner[2];              //  Corner at xmax and ymin
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 7)] = this->frecv_corner[3];   //  Corner at xmax and ymax
        }
#endif
    }

    template<class T>
    void D2Q9<T>::iSynchronize() {
#ifdef _USE_MPI_DEFINES
        int idx, neib = 0;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin  
                idx = this->Index(this->nx - 1, j);
                this->fsend_xmin[j*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 1)];
                this->fsend_xmin[j*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend_xmin[j*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)];

                //  Edge along xmax
                idx = this->Index(0, j);
                this->fsend_xmax[j*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 3)];
                this->fsend_xmax[j*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 6)];
                this->fsend_xmax[j*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 7)]; 
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, this->ny - 1);
                this->fsend_ymin[i*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 2)];
                this->fsend_ymin[i*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 5)];
                this->fsend_ymin[i*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 6)]; 

                //  Edge along ymax
                idx = this->Index(i, 0);
                this->fsend_ymax[i*3 + 0] = this->f[D2Q9<T>::IndexF(idx, 4)];
                this->fsend_ymax[i*3 + 1] = this->f[D2Q9<T>::IndexF(idx, 7)];
                this->fsend_ymax[i*3 + 2] = this->f[D2Q9<T>::IndexF(idx, 8)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->fsend_corner[0] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 5)];   //  Corner at xmin and ymin
            this->fsend_corner[1] = this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 8)];              //  Corner at xmin and ymax
            this->fsend_corner[2] = this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 6)];              //  Corner at xmax and ymin
            this->fsend_corner[3] = this->f[D2Q9<T>::IndexF(this->Index(0, 0), 7)];                         //  Corner at xmax and ymax
        }
        
        //  Communicate with other PE
        if (this->mx != 1) {
            //  Left
            MPI_Isend(this->fsend_xmin, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right
            MPI_Isend(this->fsend_xmax, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin, this->ny*3, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  Bottom
            MPI_Isend(this->fsend_ymin, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  Top
            MPI_Isend(this->fsend_ymax, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin, this->nx*3, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  Left bottom
            MPI_Isend(&this->fsend_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  Left top
            MPI_Isend(&this->fsend_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right bottom
            MPI_Irecv(&this->frecv_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Isend(&this->fsend_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);

            //  Right top
            MPI_Irecv(&this->frecv_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Isend(&this->fsend_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);      
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge along xmin
                idx = this->Index(0, j);
                this->f[D2Q9<T>::IndexF(idx, 3)] = this->frecv_xmin[j*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv_xmin[j*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv_xmin[j*3 + 2]; 
            
                //  Edge along xmax
                idx = this->Index(this->nx - 1, j);
                this->f[D2Q9<T>::IndexF(idx, 1)] = this->frecv_xmax[j*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv_xmax[j*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv_xmax[j*3 + 2];
            }
        }
        if (this->my != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge along ymin
                idx = this->Index(i, 0);
                this->f[D2Q9<T>::IndexF(idx, 4)] = this->frecv_ymin[i*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 7)] = this->frecv_ymin[i*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 8)] = this->frecv_ymin[i*3 + 2];
            
                //  Edge along ymax
                idx = this->Index(i, this->ny - 1);
                this->f[D2Q9<T>::IndexF(idx, 2)] = this->frecv_ymax[i*3 + 0];
                this->f[D2Q9<T>::IndexF(idx, 5)] = this->frecv_ymax[i*3 + 1];
                this->f[D2Q9<T>::IndexF(idx, 6)] = this->frecv_ymax[i*3 + 2];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            this->f[D2Q9<T>::IndexF(this->Index(0, 0), 7)] = this->frecv_corner[0];                         //  Corner at xmin and ymin
            this->f[D2Q9<T>::IndexF(this->Index(0, this->ny - 1), 6)] = this->frecv_corner[1];              //  Corner at xmin and ymax
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, 0), 8)] = this->frecv_corner[2];              //  Corner at xmax and ymin
            this->f[D2Q9<T>::IndexF(this->Index(this->nx - 1, this->ny - 1), 5)] = this->frecv_corner[3];   //  Corner at xmax and ymax
        }
#endif
    }
}