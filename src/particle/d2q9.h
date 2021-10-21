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
#ifdef _USE_AVX_DEFINES
    #include <immintrin.h>
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
#ifdef _USE_AVX_DEFINES
            this->f0 = (T*)_mm_malloc(sizeof(T)*this->nxyz, 32);
            this->f = (T*)_mm_malloc(sizeof(T)*this->nxyz*(D2Q9<T>::nc - 1), 32);
            this->fnext = (T*)_mm_malloc(sizeof(T)*this->nxyz*(D2Q9<T>::nc - 1), 32);
#else
            this->f0 = new T[this->nxyz];
            this->f = new T[this->nxyz*(D2Q9<T>::nc - 1)];
            this->fnext = new T[this->nxyz*(D2Q9<T>::nc - 1)];
#endif
            this->fsend_xmin = new T[this->ny*3];
            this->fsend_xmax = new T[this->ny*3];
            this->fsend_ymin = new T[this->nx*3];
            this->fsend_ymax = new T[this->nx*3];
            this->frecv_xmin = new T[this->ny*3];
            this->frecv_xmax = new T[this->ny*3];
            this->frecv_ymin = new T[this->nx*3];
            this->frecv_ymax = new T[this->nx*3];
#ifdef _USE_AVX_DEFINES
            D2Q9<T>::LoadCxCyCzEi();
#endif
        }
        D2Q9(const D2Q9<T>& _p) = delete;
        ~D2Q9() {
#ifdef _USE_AVX_DEFINES
            _mm_free(this->f0);
            _mm_free(this->f);
            _mm_free(this->fnext);
#else
            delete[] this->f0;
            delete[] this->f;
            delete[] this->fnext;
#endif
            delete[] this->fsend_xmin;
            delete[] this->fsend_xmax;
            delete[] this->fsend_ymin;
            delete[] this->fsend_ymax;
            delete[] this->frecv_xmin;
            delete[] this->frecv_xmax;
            delete[] this->frecv_ymin;
            delete[] this->frecv_ymax;
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
            return (D2Q9<T>::nc - 1)*_idx + (_c - 1);
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
        void iStream();

        void SmoothCornerAt(int _i, int _j, int _nx, int _ny);
        
        template<class Ff>
        void BoundaryCondition(Ff _bctype);
        template<class Ff>
        void iBoundaryCondition(Ff _bctype);
        void SmoothCorner() {
            //  Corner at xmin, ymin
            if (this->PEx == 0 && this->PEy == 0) {
                this->SmoothCornerAt(0, 0, -1, -1);
            }
            //  Corner at xmin, ymax
            if (this->PEx == 0 && this->PEy == this->my - 1) {
                this->SmoothCornerAt(0, this->ny - 1, -1, 1);
            }
            //  Corner at xmax, ymin
            if (this->PEx == this->mx - 1 && this->PEy == 0) {
                this->SmoothCornerAt(this->nx - 1, 0, 1, -1);
            }
            //  Corner at xmax, ymax
            if (this->PEx == this->mx - 1 && this->PEy == this->my - 1) {
                this->SmoothCornerAt(this->nx - 1, this->ny - 1, 1, 1);
            }
        }

        const int lx, ly, lz, PEid, mx, my, mz, PEx, PEy, PEz, nx, ny, nz, nxyz, offsetx, offsety, offsetz;
        static const int nc = 9, nd = 2, cx[nc], cy[nc], cz[nc];
        static const T ei[nc];
        T *f0, *f;

#ifdef _USE_AVX_DEFINES
        static const int packsize = 32/sizeof(T);
        static __m256d __cx[nc], __cy[nc], __cz[nc], __ei[nc];  //  If you use any type except double, cast these values.
        static void LoadCxCyCzEi(); 
        template<class mmT>
        void LoadF(int _idx, mmT *__f);
        template<class mmT>
        void StoreF(int _idx, const mmT *__f);
#endif

private:
        T *fnext;
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
                int idx = this->Index(i, j);
                for (int c = 1; c < D2Q9<T>::nc; ++c) {
                    int idxstream = this->Index(i - D2Q9<T>::cx[c], j - D2Q9<T>::cy[c]);
                    this->fnext[D2Q9<T>::IndexF(idx, c)] = this->f[D2Q9<T>::IndexF(idxstream, c)];
                }
            }
        }

        //  Swap
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;

#ifdef _USE_MPI_DEFINES
        int idx;

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
        int neib = 0;
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
    void D2Q9<T>::iStream() {
        //  Stream
#pragma omp parallel for
        for (int j = 0; j < this->ny; ++j) {
            for (int i = 0; i < this->nx; ++i) {
                int idx = this->Index(i, j);
                for (int c = 1; c < D2Q9<T>::nc; ++c) {
                    int idxstream = this->Index(i + D2Q9<T>::cx[c], j + D2Q9<T>::cy[c]);
                    this->fnext[D2Q9<T>::IndexF(idx, c)] = this->f[D2Q9<T>::IndexF(idxstream, c)];
                }
            }
        }

        //  Swap
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;

#ifdef _USE_MPI_DEFINES
        int idx;

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
        int neib = 0;
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

    template<class T>
    void D2Q9<T>::SmoothCornerAt(int _i, int _j, int _nx, int _ny) {
        int idx = this->Index(_i, _j), idxx = this->Index(_i - _nx, _j), idxy = this->Index(_i, _j - _ny);
        this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
        for (int c = 1; c < D2Q9<T>::nc; ++c) {
            this->f[D2Q9<T>::IndexF(idx, c)] = 0.5*(this->f[D2Q9<T>::IndexF(idxx, c)] + this->f[D2Q9<T>::IndexF(idxy, c)]);
        }
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

#ifdef _USE_AVX_DEFINES
    template<>__m256d D2Q9<double>::__cx[D2Q9<double>::nc] = { 0 };
    template<>__m256d D2Q9<double>::__cy[D2Q9<double>::nc] = { 0 };
    template<>__m256d D2Q9<double>::__cz[D2Q9<double>::nc] = { 0 };
    template<>__m256d D2Q9<double>::__ei[D2Q9<double>::nc] = { 0 };

    template<>
    void D2Q9<double>::LoadCxCyCzEi() {
        for (int c = 0; c < D2Q9<double>::nc; ++c) {
            D2Q9<double>::__cx[c] = _mm256_set1_pd((double)D2Q9<double>::cx[c]);
            D2Q9<double>::__cy[c] = _mm256_set1_pd((double)D2Q9<double>::cy[c]);
            D2Q9<double>::__cz[c] = _mm256_set1_pd((double)D2Q9<double>::cz[c]);
            D2Q9<double>::__ei[c] = _mm256_set1_pd((double)D2Q9<double>::ei[c]);
        }
    }

    template<>
    template<>
    void D2Q9<double>::LoadF<__m256d>(int _idx, __m256d *__f) {
        //  fijkl_m ijkl : c, m : idx
        const int offsetf = D2Q9<double>::IndexF(_idx, 1);
        __m256d f4321_0 = _mm256_load_pd(&this->f[offsetf +  0]);   //  f4(0) f3(0) f2(0) f1(0)
        __m256d f8765_0 = _mm256_load_pd(&this->f[offsetf +  4]);   //  f8(0) f7(0) f6(0) f5(0)  
        __m256d f4321_1 = _mm256_load_pd(&this->f[offsetf +  8]);   //  f4(1) f3(1) f2(1) f1(1)
        __m256d f8765_1 = _mm256_load_pd(&this->f[offsetf + 12]);   //  f8(1) f7(1) f6(1) f5(1)
        __m256d f4321_2 = _mm256_load_pd(&this->f[offsetf + 16]);   //  f4(2) f3(2) f2(2) f1(2)
        __m256d f8765_2 = _mm256_load_pd(&this->f[offsetf + 20]);   //  f8(2) f7(2) f6(2) f5(2)
        __m256d f4321_3 = _mm256_load_pd(&this->f[offsetf + 24]);   //  f4(3) f3(3) f2(3) f1(3)
        __m256d f8765_3 = _mm256_load_pd(&this->f[offsetf + 28]);   //  f8(3) f7(3) f6(3) f5(3)

        //  fij_kl ij : c, kl : idx
        __m256d f31_10 = _mm256_unpacklo_pd(f4321_0, f4321_1);      //  f3(1) f3(0) f1(1) f1(0)
        __m256d f42_10 = _mm256_unpackhi_pd(f4321_0, f4321_1);      //  f4(1) f4(0) f2(1) f2(0)
        __m256d f75_10 = _mm256_unpacklo_pd(f8765_0, f8765_1);      //  f7(1) f7(0) f5(1) f5(0)
        __m256d f86_10 = _mm256_unpackhi_pd(f8765_0, f8765_1);      //  f8(1) f8(0) f6(1) f6(0)
        __m256d f31_32 = _mm256_unpacklo_pd(f4321_2, f4321_3);      //  f3(3) f3(2) f1(3) f1(2)
        __m256d f42_32 = _mm256_unpackhi_pd(f4321_2, f4321_3);      //  f4(3) f4(2) f2(3) f2(2)
        __m256d f75_32 = _mm256_unpacklo_pd(f8765_2, f8765_3);      //  f7(3) f7(2) f5(3) f5(2)
        __m256d f86_32 = _mm256_unpackhi_pd(f8765_2, f8765_3);      //  f8(3) f6(2) f8(3) f6(2)

        //  fi i : c
        const int mm0 = 2*16 + 0*1, mm1 = 3*16 + 1*1;
        __f[0] = _mm256_load_pd(&this->f0[_idx]);                   //  f0(3) f0(2) f0(1) f0(0)
        __f[1] = _mm256_permute2f128_pd(f31_10, f31_32, mm0);       //  f1(3) f1(2) f1(1) f1(0)
        __f[2] = _mm256_permute2f128_pd(f42_10, f42_32, mm0);       //  f2(3) f2(2) f2(1) f2(0)
        __f[3] = _mm256_permute2f128_pd(f31_10, f31_32, mm1);       //  f3(3) f3(2) f3(1) f3(0)
        __f[4] = _mm256_permute2f128_pd(f42_10, f42_32, mm1);       //  f4(3) f4(2) f4(1) f4(0)
        __f[5] = _mm256_permute2f128_pd(f75_10, f75_32, mm0);       //  f5(3) f5(2) f5(1) f5(0)
        __f[6] = _mm256_permute2f128_pd(f86_10, f86_32, mm0);       //  f6(3) f6(2) f6(1) f6(0)
        __f[7] = _mm256_permute2f128_pd(f75_10, f75_32, mm1);       //  f7(3) f7(2) f7(1) f7(0)
        __f[8] = _mm256_permute2f128_pd(f86_10, f86_32, mm1);       //  f8(3) f8(2) f8(1) f8(0)
    }

    template<>
    template<>
    void D2Q9<double>::StoreF<__m256d>(int _idx, const __m256d *__f) {
        const int mm0 = 2*16 + 0*1, mm1 = 3*16 + 1*1;
        __m256d f31_10 = _mm256_permute2f128_pd(__f[1], __f[3], mm0);                   //  f3(1) f3(0) f1(1) f1(0)
        __m256d f42_10 = _mm256_permute2f128_pd(__f[2], __f[4], mm0);                   //  f4(1) f4(0) f2(1) f2(0)
        __m256d f75_10 = _mm256_permute2f128_pd(__f[5], __f[7], mm0);                   //  f7(1) f7(0) f5(1) f5(0)
        __m256d f86_10 = _mm256_permute2f128_pd(__f[6], __f[8], mm0);                   //  f8(1) f6(0) f8(1) f6(0)
        __m256d f31_32 = _mm256_permute2f128_pd(__f[1], __f[3], mm1);                   //  f3(3) f3(2) f1(3) f1(2)
        __m256d f42_32 = _mm256_permute2f128_pd(__f[2], __f[4], mm1);                   //  f4(3) f4(2) f2(3) f2(2)
        __m256d f75_32 = _mm256_permute2f128_pd(__f[5], __f[7], mm1);                   //  f7(3) f7(2) f5(3) f5(2)
        __m256d f86_32 = _mm256_permute2f128_pd(__f[6], __f[8], mm1);                   //  f8(3) f6(2) f8(3) f6(2)

        const int offsetf = D2Q9<double>::IndexF(_idx, 1);
        _mm256_store_pd(&this->f0[_idx], __f[0]);                                       //  f0(3) f0(2) f0(1) f0(0)
        _mm256_store_pd(&this->f[offsetf +  0], _mm256_unpacklo_pd(f31_10, f42_10));    //  f4(0) f3(0) f2(0) f1(0)
        _mm256_store_pd(&this->f[offsetf +  4], _mm256_unpacklo_pd(f75_10, f86_10));    //  f8(0) f7(0) f6(0) f5(0)  
        _mm256_store_pd(&this->f[offsetf +  8], _mm256_unpackhi_pd(f31_10, f42_10));    //  f4(1) f3(1) f2(1) f1(1)
        _mm256_store_pd(&this->f[offsetf + 12], _mm256_unpackhi_pd(f75_10, f86_10));    //  f8(1) f7(1) f6(1) f5(1)
        _mm256_store_pd(&this->f[offsetf + 16], _mm256_unpacklo_pd(f31_32, f42_32));    //  f4(2) f3(2) f2(2) f1(2)
        _mm256_store_pd(&this->f[offsetf + 20], _mm256_unpacklo_pd(f75_32, f86_32));    //  f8(2) f7(2) f6(2) f5(2)
        _mm256_store_pd(&this->f[offsetf + 24], _mm256_unpackhi_pd(f31_32, f42_32));    //  f4(3) f3(3) f2(3) f1(3)
        _mm256_store_pd(&this->f[offsetf + 28], _mm256_unpackhi_pd(f75_32, f86_32));    //  f8(3) f7(3) f6(3) f5(3)
    }
#endif
}