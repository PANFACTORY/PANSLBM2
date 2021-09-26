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
#ifdef _USE_AVX_DEFINES
    #include <immintrin.h>
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
            nxyz(this->nx*this->ny*this->nz),
            offsetx(this->mx - this->PEx > this->lx%this->mx ? this->PEx*this->nx : this->lx - (this->mx - this->PEx)*this->nx),
            offsety(this->my - this->PEy > this->ly%this->my ? this->PEy*this->ny : this->ly - (this->my - this->PEy)*this->ny),
            offsetz(this->mz - this->PEz > this->lz%this->mz ? this->PEz*this->nz : this->lz - (this->mz - this->PEz)*this->nz)
        {
            assert(0 < _lx && 0 < _ly && 0 < _lz && 0 <= _PEid && 0 < _mx && 0 < _my && 0 < _mz);
            this->f0 = new(std::align_val_t{32}) T[(this->nxyz/D3Q15<T>::packsize + 1)*D3Q15<T>::packsize];
            this->f = new(std::align_val_t{32}) T[(this->nxyz/D3Q15<T>::packsize + 1)*D3Q15<T>::packsize*(D3Q15<T>::nc - 1)];
            this->fnext = new(std::align_val_t{32}) T[(this->nxyz/D3Q15<T>::packsize + 1)*D3Q15<T>::packsize*(D3Q15<T>::nc - 1)];
            this->fsend_xmin = new T[this->ny*this->nz*5];
            this->fsend_xmax = new T[this->ny*this->nz*5];
            this->fsend_ymin = new T[this->nz*this->nx*5];
            this->fsend_ymax = new T[this->nz*this->nx*5];
            this->fsend_zmin = new T[this->nx*this->ny*5];
            this->fsend_zmax = new T[this->nx*this->ny*5];
            this->frecv_xmin = new T[this->ny*this->nz*5];
            this->frecv_xmax = new T[this->ny*this->nz*5];
            this->frecv_ymin = new T[this->nz*this->nx*5];
            this->frecv_ymax = new T[this->nz*this->nx*5];
            this->frecv_zmin = new T[this->nx*this->ny*5];
            this->frecv_zmax = new T[this->nx*this->ny*5];
            this->fsend_ymin_zmin = new T[this->nx*2];
            this->fsend_ymin_zmax = new T[this->nx*2];
            this->fsend_ymax_zmin = new T[this->nx*2];
            this->fsend_ymax_zmax = new T[this->nx*2];
            this->frecv_ymin_zmin = new T[this->nx*2];
            this->frecv_ymin_zmax = new T[this->nx*2];
            this->frecv_ymax_zmin = new T[this->nx*2];
            this->frecv_ymax_zmax = new T[this->nx*2]; 
            this->fsend_zmin_xmin = new T[this->ny*2];
            this->fsend_zmin_xmax = new T[this->ny*2];
            this->fsend_zmax_xmin = new T[this->ny*2];
            this->fsend_zmax_xmax = new T[this->ny*2];
            this->frecv_zmin_xmin = new T[this->ny*2];
            this->frecv_zmin_xmax = new T[this->ny*2];
            this->frecv_zmax_xmin = new T[this->ny*2];
            this->frecv_zmax_xmax = new T[this->ny*2];
            this->fsend_xmin_ymin = new T[this->nz*2];
            this->fsend_xmin_ymax = new T[this->nz*2];
            this->fsend_xmax_ymin = new T[this->nz*2];
            this->fsend_xmax_ymax = new T[this->nz*2];
            this->frecv_xmin_ymin = new T[this->nz*2];
            this->frecv_xmin_ymax = new T[this->nz*2];
            this->frecv_xmax_ymin = new T[this->nz*2];
            this->frecv_xmax_ymax = new T[this->nz*2];
#ifdef _USE_AVX_DEFINES
            for (int c = 0; c < D3Q15<T>::nc; ++c) {
                D3Q15<T>::__cx[c] = _mm256_set1_pd((double)D3Q15<T>::cx[c]);
                D3Q15<T>::__cy[c] = _mm256_set1_pd((double)D3Q15<T>::cy[c]);
                D3Q15<T>::__cz[c] = _mm256_set1_pd((double)D3Q15<T>::cz[c]);
                D3Q15<T>::__ei[c] = _mm256_set1_pd(D3Q15<T>::ei[c]);
            }
#endif
        }
        D3Q15(const D3Q15<T>& _p) = delete;
        ~D3Q15() {
            delete[] this->f0, this->f, this->fnext;
            delete[] this->fsend_xmin, this->fsend_xmax, this->fsend_ymin, this->fsend_ymax, this->fsend_zmin, this->fsend_zmax;
            delete[] this->frecv_xmin, this->frecv_xmax, this->frecv_ymin, this->frecv_ymax, this->frecv_zmin, this->frecv_zmax;
            delete[] this->fsend_ymin_zmin, this->fsend_ymin_zmax, this->fsend_ymax_zmin, this->fsend_ymax_zmax;
            delete[] this->frecv_ymin_zmin, this->frecv_ymin_zmax, this->frecv_ymax_zmin, this->frecv_ymax_zmax; 
            delete[] this->fsend_zmin_xmin, this->fsend_zmin_xmax, this->fsend_zmax_xmin, this->fsend_zmax_xmax;
            delete[] this->frecv_zmin_xmin, this->frecv_zmin_xmax, this->frecv_zmax_xmin, this->frecv_zmax_xmax;
            delete[] this->fsend_xmin_ymin, this->fsend_xmin_ymax, this->fsend_xmax_ymin, this->fsend_xmax_ymax;
            delete[] this->frecv_xmin_ymin, this->frecv_xmin_ymax, this->frecv_xmax_ymin, this->frecv_xmax_ymax;  
        }

        int Index(int _i, int _j, int _k) const {
            int i = _i == -1 ? this->nx - 1 : (_i == this->nx ? 0 : _i);
            int j = _j == -1 ? this->ny - 1 : (_j == this->ny ? 0 : _j);
            int k = _k == -1 ? this->nz - 1 : (_k == this->nz ? 0 : _k);
            return i + this->nx*j + this->nx*this->ny*k;
        }
        static int IndexF(int _idx, int _c) {
            return (_idx/D3Q15<T>::packsize)*D3Q15<T>::packsize*(D3Q15<T>::nc - 1) + D3Q15<T>::packsize*(_c - 1) + _idx%D3Q15<T>::packsize;
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

        void Stream();
        void iStream();

        template<class Ff>
        void BoundaryCondition(Ff _bctype);
        template<class Ff>
        void iBoundaryCondition(Ff _bctype);
        void SmoothCorner();

        const int lx, ly, lz, PEid, mx, my, mz, PEx, PEy, PEz, nx, ny, nz, nxyz, offsetx, offsety, offsetz;
        static const int nc = 15, nd = 3, cx[nc], cy[nc], cz[nc];
        static const T ei[nc];
        T *f0, *f;
#ifdef _USE_AVX_DEFINES
        static const int packsize = 32/sizeof(T);
        static __m256d __cx[nc], __cy[nc], __cz[nc], __ei[nc];
#else
        static const int packsize = 1;
#endif
        
private:
        T *fnext;
        T *fsend_xmin, *fsend_xmax, *fsend_ymin, *fsend_ymax, *fsend_zmin, *fsend_zmax, *frecv_xmin, *frecv_xmax, *frecv_ymin, *frecv_ymax, *frecv_zmin, *frecv_zmax;
        T *fsend_ymin_zmin, *fsend_ymin_zmax, *fsend_ymax_zmin, *fsend_ymax_zmax, *frecv_ymin_zmin, *frecv_ymin_zmax, *frecv_ymax_zmin, *frecv_ymax_zmax; 
        T *fsend_zmin_xmin, *fsend_zmin_xmax, *fsend_zmax_xmin, *fsend_zmax_xmax, *frecv_zmin_xmin, *frecv_zmin_xmax, *frecv_zmax_xmin, *frecv_zmax_xmax;
        T *fsend_xmin_ymin, *fsend_xmin_ymax, *fsend_xmax_ymin, *fsend_xmax_ymax, *frecv_xmin_ymin, *frecv_xmin_ymax, *frecv_xmax_ymin, *frecv_xmax_ymax;  
        T fsend_corner[8], frecv_corner[8];
#ifdef _USE_MPI_DEFINES
        MPI_Status status[52];
        MPI_Request request[52];
#endif
    };

    template<class T>const int D3Q15<T>::cx[D3Q15<T>::nc] = { 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, -1, -1 };
    template<class T>const int D3Q15<T>::cy[D3Q15<T>::nc] = { 0, 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, -1, -1, 1, -1 };
    template<class T>const int D3Q15<T>::cz[D3Q15<T>::nc] = { 0, 0, 0, 1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1 };
    template<class T>const T D3Q15<T>::ei[D3Q15<T>::nc] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0 };

#ifdef _USE_AVX_DEFINES
    template<class T>__m256d D2Q9<T>::__cx[D2Q9<T>::nc] = { 0 };
    template<class T>__m256d D2Q9<T>::__cy[D2Q9<T>::nc] = { 0 };
    template<class T>__m256d D2Q9<T>::__cz[D2Q9<T>::nc] = { 0 };
    template<class T>__m256d D2Q9<T>::__ei[D2Q9<T>::nc] = { 0 };
#endif

    template<class T>
    void D3Q15<T>::Stream() {
        //  Stream
#pragma omp parallel for
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idx = this->Index(i, j, k), idxstream = this->Index(i - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        }
        
        //  Swap
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;

#ifdef _USE_MPI_DEFINES
        int idx, idxface;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->fsend_xmin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->fsend_xmin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_xmin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_xmin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->fsend_xmin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on xmax
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->fsend_xmax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->fsend_xmax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_xmax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_xmax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_xmax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 12)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i);
                    this->fsend_ymin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->fsend_ymin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_ymin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_ymin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend_ymin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on ymax
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i);
                    this->fsend_ymax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->fsend_ymax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_ymax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_ymax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_ymax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j);
                    this->fsend_zmin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->fsend_zmin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_zmin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_zmin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend_zmin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];

                    //  Face on zmax
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j);
                    this->fsend_zmax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->fsend_zmax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_zmax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_zmax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_zmax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                this->fsend_ymin_zmin[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_ymin_zmin[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on ymin and zmax
                idx = this->Index(i, this->ny - 1, 0);
                this->fsend_ymin_zmax[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend_ymin_zmax[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on ymax and zmin
                idx = this->Index(i, 0, this->nz - 1);
                this->fsend_ymax_zmin[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend_ymax_zmin[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on ymax and zmax
                idx = this->Index(i, 0, 0);
                this->fsend_ymax_zmax[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_ymax_zmax[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                this->fsend_zmin_xmin[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_zmin_xmin[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on zmin and xmax
                idx = this->Index(0, j, this->nz - 1);
                this->fsend_zmin_xmax[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend_zmin_xmax[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on zmax and xmin
                idx = this->Index(this->nx - 1, j, 0);
                this->fsend_zmax_xmin[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend_zmax_xmin[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on zmax and xmax
                idx = this->Index(0, j, 0);
                this->fsend_zmax_xmax[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_zmax_xmax[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                this->fsend_xmin_ymin[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_xmin_ymin[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on xmin and ymax
                idx = this->Index(this->nx - 1, 0, k);
                this->fsend_xmin_ymax[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend_xmin_ymax[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on xmax and ymin
                idx = this->Index(0, this->ny - 1, k);
                this->fsend_xmax_ymin[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend_xmax_ymin[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on xmax and ymax
                idx = this->Index(0, 0, k);
                this->fsend_xmax_ymax[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_xmax_ymax[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend_corner[0] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)];   //  Corner at xmin, ymin and zmin
            this->fsend_corner[1] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)];              //  Corner at xmax, ymin and zmin
            this->fsend_corner[2] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)];              //  Corner at xmin, ymax and zmin 
            this->fsend_corner[3] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)];                         //  Corner at xmax, ymax and zmin
            this->fsend_corner[4] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)];              //  Corner at xmin, ymin and zmax
            this->fsend_corner[5] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)];                          //  Corner at xmax, ymin and zmax
            this->fsend_corner[6] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)];                          //  Corner at xmin, ymax and zmax
            this->fsend_corner[7] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)];                                     //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        int neib = 0;
        if (this->mx != 1) {
            //  To xmin
            MPI_Isend(this->fsend_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  To xmax
            MPI_Isend(this->fsend_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  To ymin
            MPI_Isend(this->fsend_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  To ymax
            MPI_Isend(this->fsend_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1) {
            //  To zmin
            MPI_Isend(this->fsend_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  To zmax
            MPI_Isend(this->fsend_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1 || this->mz != 1) {
            //  To ymin and zmin
            MPI_Isend(this->fsend_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymin and zmax
            MPI_Isend(this->fsend_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmin
            MPI_Isend(this->fsend_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmax
            MPI_Isend(this->fsend_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1 || this->mx != 1) {
            //  To zmin and xmin
            MPI_Isend(this->fsend_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To zmin and xmax
            MPI_Isend(this->fsend_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
                 
            //  To zmax and xmin
            MPI_Isend(this->fsend_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To zmax and xmax
            MPI_Isend(this->fsend_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  To xmin and ymin
            MPI_Isend(this->fsend_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To xmin and ymax
            MPI_Isend(this->fsend_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymin
            MPI_Isend(this->fsend_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymax
            MPI_Isend(this->fsend_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            MPI_Isend(&this->fsend_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmin
            MPI_Isend(&this->fsend_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmin
            MPI_Isend(&this->fsend_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmin
            MPI_Isend(&this->fsend_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmin
            MPI_Isend(&this->fsend_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmax
            MPI_Isend(&this->fsend_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmax
            MPI_Isend(&this->fsend_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmax
            MPI_Isend(&this->fsend_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmax
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
                    idxface = this->IndexBCx(j, k);
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->frecv_xmin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_xmin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_xmin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_xmin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_xmin[idxface*5 + 4];

                    //  Face on xmax
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->frecv_xmax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_xmax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_xmax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_xmax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_xmax[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i);
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->frecv_ymin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_ymin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_ymin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_ymin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_ymin[idxface*5 + 4];

                    //  Face on ymax
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i);
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->frecv_ymax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_ymax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_ymax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_ymax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_ymax[idxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j);
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->frecv_zmin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_zmin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_zmin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_zmin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_zmin[idxface*5 + 4];

                    //  Face on zmax
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j);
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->frecv_zmax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_zmax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_zmax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_zmax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_zmax[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, 0, 0);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_ymin_zmin[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_ymin_zmin[i*2 + 1];

                //  Edge on ymin and zmax
                idx = this->Index(i, 0, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_ymin_zmax[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_ymin_zmax[i*2 + 1];

                //  Edge on ymax and zmin
                idx = this->Index(i, this->ny - 1, 0);
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_ymax_zmin[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_ymax_zmin[i*2 + 1];

                //  Edge on ymax and zmax
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_ymax_zmax[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_ymax_zmax[i*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(0, j, 0);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_zmin_xmin[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_zmin_xmin[j*2 + 1];

                //  Edge on zmin and xmax
                idx = this->Index(this->nx - 1, j, 0);
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_zmin_xmax[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_zmin_xmax[j*2 + 1];

                //  Edge on zmax and xmin
                idx = this->Index(0, j, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_zmax_xmin[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_zmax_xmin[j*2 + 1];

                //  Edge on zmax and xmax
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_zmax_xmax[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_zmax_xmax[j*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(0, 0, k);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_xmin_ymin[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_xmin_ymin[k*2 + 1];

                //  Edge on xmin and ymax
                idx = this->Index(0, this->ny - 1, k);
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_xmin_ymax[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_xmin_ymax[k*2 + 1];

                //  Edge on xmax and ymin
                idx = this->Index(this->nx - 1, 0, k);
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_xmax_ymin[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_xmax_ymin[k*2 + 1];

                //  Edge on xmax and ymax
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_xmax_ymax[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_xmax_ymax[k*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)] = this->frecv_corner[0];                                     //  Corner at xmin, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)] = this->frecv_corner[1];                          //  Corner at xmax, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)] = this->frecv_corner[2];                          //  Corner at xmin, ymax and zmin 
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)] = this->frecv_corner[3];              //  Corner at xmax, ymax and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)] = this->frecv_corner[4];                         //  Corner at xmin, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)] = this->frecv_corner[5];              //  Corner at xmax, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)] = this->frecv_corner[6];              //  Corner at xmin, ymax and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)] = this->frecv_corner[7];   //  Corner at xmax, ymax and zmax
        }
#endif
    }

    template<class T>
    void D3Q15<T>::iStream() {
        //  Stream
#pragma omp parallel for
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idx = this->Index(i, j, k), idxstream = this->Index(i + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        }
        
        //  Swap
        T *tmp = this->f;
        this->f = this->fnext;
        this->fnext = tmp;

#ifdef _USE_MPI_DEFINES
        int idx, idxface;

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->fsend_xmin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 1)];
                    this->fsend_xmin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_xmin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_xmin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_xmin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 12)];

                    //  Face on xmax
                    idx = this->Index(0, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->fsend_xmax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 4)];
                    this->fsend_xmax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_xmax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_xmax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 13)];
                    this->fsend_xmax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i);
                    this->fsend_ymin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 2)];
                    this->fsend_ymin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_ymin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_ymin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_ymin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];

                    //  Face on ymax
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i);
                    this->fsend_ymax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 5)];
                    this->fsend_ymax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_ymax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_ymax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend_ymax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j);
                    this->fsend_zmin[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 3)];
                    this->fsend_zmin[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 7)];
                    this->fsend_zmin[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    this->fsend_zmin[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    this->fsend_zmin[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 14)];

                    //  Face on zmax
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j);
                    this->fsend_zmax[idxface*5 + 0] = this->f[D3Q15<T>::IndexF(idx, 6)];
                    this->fsend_zmax[idxface*5 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    this->fsend_zmax[idxface*5 + 2] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    this->fsend_zmax[idxface*5 + 3] = this->f[D3Q15<T>::IndexF(idx, 12)];
                    this->fsend_zmax[idxface*5 + 4] = this->f[D3Q15<T>::IndexF(idx, 13)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                this->fsend_ymin_zmin[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_ymin_zmin[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 8)];

                //  Edge on ymin and zmax
                idx = this->Index(i, this->ny - 1, 0);
                this->fsend_ymin_zmax[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend_ymin_zmax[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on ymax and zmin
                idx = this->Index(i, 0, this->nz - 1);
                this->fsend_ymax_zmin[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend_ymax_zmin[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on ymax and zmax
                idx = this->Index(i, 0, 0);
                this->fsend_ymax_zmax[i*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_ymax_zmax[i*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                this->fsend_zmin_xmin[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_zmin_xmin[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 9)];

                //  Edge on zmin and xmax
                idx = this->Index(0, j, this->nz - 1);
                this->fsend_zmin_xmax[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend_zmin_xmax[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];

                //  Edge on zmax and xmin
                idx = this->Index(this->nx - 1, j, 0);
                this->fsend_zmax_xmin[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 10)];
                this->fsend_zmax_xmin[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];

                //  Edge on zmax and xmax
                idx = this->Index(0, j, 0);
                this->fsend_zmax_xmax[j*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_zmax_xmax[j*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                this->fsend_xmin_ymin[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 7)];
                this->fsend_xmin_ymin[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 10)];

                //  Edge on xmin and ymax
                idx = this->Index(this->nx - 1, 0, k);
                this->fsend_xmin_ymax[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 9)];
                this->fsend_xmin_ymax[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 12)];
                
                //  Edge on xmax and ymin
                idx = this->Index(0, this->ny - 1, k);
                this->fsend_xmax_ymin[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 8)];
                this->fsend_xmax_ymin[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 13)];

                //  Edge on xmax and ymax
                idx = this->Index(0, 0, k);
                this->fsend_xmax_ymax[k*2 + 0] = this->f[D3Q15<T>::IndexF(idx, 11)];
                this->fsend_xmax_ymax[k*2 + 1] = this->f[D3Q15<T>::IndexF(idx, 14)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend_corner[0] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 7)];    //  Corner at xmin, ymin and zmin
            this->fsend_corner[1] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 8)];               //  Corner at xmax, ymin and zmin
            this->fsend_corner[2] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 9)];               //  Corner at xmin, ymax and zmin 
            this->fsend_corner[3] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 14)];                         //  Corner at xmax, ymax and zmin
            this->fsend_corner[4] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 10)];              //  Corner at xmin, ymin and zmax
            this->fsend_corner[5] = this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 13)];                         //  Corner at xmax, ymin and zmax
            this->fsend_corner[6] = this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 12)];                         //  Corner at xmin, ymax and zmax
            this->fsend_corner[7] = this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 11)];                                    //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        int neib = 0;
        if (this->mx != 1) {
            //  To xmin
            MPI_Isend(this->fsend_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &this->request[neib++]);

            //  To xmax
            MPI_Isend(this->fsend_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 1, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  To ymin
            MPI_Isend(this->fsend_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 2, MPI_COMM_WORLD, &this->request[neib++]);

            //  To ymax
            MPI_Isend(this->fsend_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 3, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1) {
            //  To zmin
            MPI_Isend(this->fsend_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 4, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 4, MPI_COMM_WORLD, &this->request[neib++]);

            //  To zmax
            MPI_Isend(this->fsend_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 5, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1 || this->mz != 1) {
            //  To ymin and zmin
            MPI_Isend(this->fsend_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 6, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymin and zmax
            MPI_Isend(this->fsend_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 7, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmin
            MPI_Isend(this->fsend_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 8, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmax
            MPI_Isend(this->fsend_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 9, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1 || this->mx != 1) {
            //  To zmin and xmin
            MPI_Isend(this->fsend_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 10, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To zmin and xmax
            MPI_Isend(this->fsend_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 11, MPI_COMM_WORLD, &this->request[neib++]);
                 
            //  To zmax and xmin
            MPI_Isend(this->fsend_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 12, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To zmax and xmax
            MPI_Isend(this->fsend_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 13, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  To xmin and ymin
            MPI_Isend(this->fsend_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 14, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To xmin and ymax
            MPI_Isend(this->fsend_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 15, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymin
            MPI_Isend(this->fsend_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 16, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymax
            MPI_Isend(this->fsend_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 17, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            MPI_Isend(&this->fsend_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmin
            MPI_Isend(&this->fsend_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmin
            MPI_Isend(&this->fsend_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmin
            MPI_Isend(&this->fsend_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmin
            MPI_Isend(&this->fsend_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmax
            MPI_Isend(&this->fsend_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmax
            MPI_Isend(&this->fsend_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmax
            MPI_Isend(&this->fsend_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmax
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
                    idxface = this->IndexBCx(j, k);
                    this->f[D3Q15<T>::IndexF(idx, 4)] = this->frecv_xmin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_xmin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_xmin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_xmin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_xmin[idxface*5 + 4];

                    //  Face on xmax
                    idx = this->Index(this->nx - 1, j, k);
                    idxface = this->IndexBCx(j, k);
                    this->f[D3Q15<T>::IndexF(idx, 1)] = this->frecv_xmax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_xmax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_xmax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_xmax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_xmax[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    idx = this->Index(i, 0, k);
                    idxface = this->IndexBCy(k, i);
                    this->f[D3Q15<T>::IndexF(idx, 5)] = this->frecv_ymin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_ymin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_ymin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_ymin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_ymin[idxface*5 + 4];

                    //  Face on ymax
                    idx = this->Index(i, this->ny - 1, k);
                    idxface = this->IndexBCy(k, i);
                    this->f[D3Q15<T>::IndexF(idx, 2)] = this->frecv_ymax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_ymax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_ymax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_ymax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_ymax[idxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    idx = this->Index(i, j, 0);
                    idxface = this->IndexBCz(i, j);
                    this->f[D3Q15<T>::IndexF(idx, 6)] = this->frecv_zmin[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_zmin[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_zmin[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_zmin[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_zmin[idxface*5 + 4];

                    //  Face on zmax
                    idx = this->Index(i, j, this->nz - 1);
                    idxface = this->IndexBCz(i, j);
                    this->f[D3Q15<T>::IndexF(idx, 3)] = this->frecv_zmax[idxface*5 + 0];
                    this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_zmax[idxface*5 + 1];
                    this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_zmax[idxface*5 + 2];
                    this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_zmax[idxface*5 + 3];
                    this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_zmax[idxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                idx = this->Index(i, 0, 0);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_ymin_zmin[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_ymin_zmin[i*2 + 1];

                //  Edge on ymin and zmax
                idx = this->Index(i, 0, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_ymin_zmax[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_ymin_zmax[i*2 + 1];

                //  Edge on ymax and zmin
                idx = this->Index(i, this->ny - 1, 0);
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_ymax_zmin[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_ymax_zmin[i*2 + 1];

                //  Edge on ymax and zmax
                idx = this->Index(i, this->ny - 1, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_ymax_zmax[i*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_ymax_zmax[i*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                idx = this->Index(0, j, 0);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_zmin_xmin[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_zmin_xmin[j*2 + 1];

                //  Edge on zmin and xmax
                idx = this->Index(this->nx - 1, j, 0);
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_zmin_xmax[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_zmin_xmax[j*2 + 1];

                //  Edge on zmax and xmin
                idx = this->Index(0, j, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_zmax_xmin[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_zmax_xmin[j*2 + 1];

                //  Edge on zmax and xmax
                idx = this->Index(this->nx - 1, j, this->nz - 1);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_zmax_xmax[j*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_zmax_xmax[j*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                idx = this->Index(0, 0, k);
                this->f[D3Q15<T>::IndexF(idx, 11)] = this->frecv_xmin_ymin[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 14)] = this->frecv_xmin_ymin[k*2 + 1];

                //  Edge on xmin and ymax
                idx = this->Index(0, this->ny - 1, k);
                this->f[D3Q15<T>::IndexF(idx, 8)] = this->frecv_xmin_ymax[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 13)] = this->frecv_xmin_ymax[k*2 + 1];

                //  Edge on xmax and ymin
                idx = this->Index(this->nx - 1, 0, k);
                this->f[D3Q15<T>::IndexF(idx, 9)] = this->frecv_xmax_ymin[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 12)] = this->frecv_xmax_ymin[k*2 + 1];

                //  Edge on xmax and ymax
                idx = this->Index(this->nx - 1, this->ny - 1, k);
                this->f[D3Q15<T>::IndexF(idx, 7)] = this->frecv_xmax_ymax[k*2 + 0];
                this->f[D3Q15<T>::IndexF(idx, 10)] = this->frecv_xmax_ymax[k*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, 0), 11)] = this->frecv_corner[0];                                    //  Corner at xmin, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 12)] = this->frecv_corner[1];                         //  Corner at xmax, ymin and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 13)] = this->frecv_corner[2];                         //  Corner at xmin, ymax and zmin 
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 10)] = this->frecv_corner[3];              //  Corner at xmax, ymax and zmin
            this->f[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 14)] = this->frecv_corner[4];                         //  Corner at xmin, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 9)] = this->frecv_corner[5];               //  Corner at xmax, ymin and zmax
            this->f[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 8)] = this->frecv_corner[6];               //  Corner at xmin, ymax and zmax
            this->f[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 7)] = this->frecv_corner[7];    //  Corner at xmax, ymax and zmax
        }
#endif
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::BoundaryCondition(Ff _bctype) {
        //  On xmin
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    if (_bctype(0 + this->offsetx, j + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(0, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    } else if (_bctype(0 + this->offsetx, j + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(0, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
        //  On xmax
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(this->nx - 1, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(this->nx - 1, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    }
                }
            }
        }
        //  On ymin
        if (this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    if (_bctype(i + this->offsetx, 0 + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, 0, k);
                        this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    } else if (_bctype(i + this->offsetx, 0 + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, 0, k);
                        this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
        //  On ymax
        if (this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, this->ny - 1, k);
                        this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, this->ny - 1, k);
                        this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    }
                }
            }
        }
        //  On zmin
        if (this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    if (_bctype(i + this->offsetx, j + this->offsety, 0 + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, j, 0);
                        this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype(i + this->offsetx, j + this->offsety, 0 + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, j, 0);
                        this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
        //  On zmax
        if (this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {                    
                    if (_bctype(i + this->offsetx, j + this->offsety, (this->nz - 1) + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, j, this->nz - 1);
                        this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    } else if (_bctype(i + this->offsetx, j + this->offsety, (this->nz - 1) + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, j, this->nz - 1);
                        this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    }
                }
            }
        }
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::iBoundaryCondition(Ff _bctype) {
        //  On xmin
        if (this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    if (_bctype(0 + this->offsetx, j + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(0, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype(0 + this->offsetx, j + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(0, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    }
                }
            }
        }
        //  On xmax
        if (this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(this->nx - 1, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    } else if (_bctype((this->nx - 1) + this->offsetx, j + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(this->nx - 1, j, k);
                        this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
        //  On ymin
        if (this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    if (_bctype(i + this->offsetx, 0 + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, 0, k);
                        this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype(i + this->offsetx, 0 + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, 0, k);
                        this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    }
                }
            }
        }
        //  On ymax
        if (this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety, k + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, this->ny - 1, k);
                        this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    } else if (_bctype(i + this->offsetx, (this->ny - 1) + this->offsety, k + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, this->ny - 1, k);
                        this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
        //  On zmin
        if (this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    if (_bctype(i + this->offsetx, j + this->offsety, 0 + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, j, 0);
                        this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                    } else if (_bctype(i + this->offsetx, j + this->offsety, 0 + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, j, 0);
                        this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                        this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                        this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                        this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                    }
                }
            }
        }
        //  On zmax
        if (this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    if (_bctype(i + this->offsetx, j + this->offsety, (this->nz - 1) + this->offsetz) == BARRIER) {
                        int idx = this->Index(i, j, this->nz - 1);
                        this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                    } else if (_bctype(i + this->offsetx, j + this->offsety, (this->nz - 1) + this->offsetz) == MIRROR) {
                        int idx = this->Index(i, j, this->nz - 1);
                        this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                        this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                        this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                        this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                    }
                }
            }
        }
    }

    template<class T>
    void D3Q15<T>::SmoothCorner() {
        //  Along line ymin and zmin
        if (this->PEy == 0 && this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                int idx = this->Index(i, 0, 0), idxy = this->Index(i, 1, 0), idxz = this->Index(i, 0, 1);
                this->f0[idx] = 0.5*(this->f0[idxy] + this->f0[idxz]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymax and zmin
        if (this->PEy == this->my - 1 && this->PEz == 0) {
            for (int i = 0; i < this->nx; ++i) {
                int idx = this->Index(i, this->ny - 1, 0), idxy = this->Index(i, this->ny - 2, 0), idxz = this->Index(i, this->ny - 1, 1);
                this->f0[idx] = 0.5*(this->f0[idxy] + this->f0[idxz]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymax and zmax
        if (this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                int idx = this->Index(i, this->ny - 1, this->nz - 1), idxy = this->Index(i, this->ny - 2, this->nz - 1), idxz = this->Index(i, this->ny - 1, this->nz - 2);
                this->f0[idx] = 0.5*(this->f0[idxy] + this->f0[idxz]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }
        //  Along line ymin and zmax
        if (this->PEy == 0 && this->PEz == this->mz - 1) {
            for (int i = 0; i < this->nx; ++i) {
                int idx = this->Index(i, 0, this->nz - 1), idxy = this->Index(i, 1, this->nz - 1), idxz = this->Index(i, 0, this->nz - 2);
                this->f0[idx] = 0.5*(this->f0[idxy] + this->f0[idxz]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                }
            }
        }

        //  Along line zmin and xmin
        if (this->PEz == 0 && this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {
                int idx = this->Index(0, j, 0), idxz = this->Index(1, j, 0), idxx = this->Index(0, j, 1);
                this->f0[idx] = 0.5*(this->f0[idxz] + this->f0[idxx]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmax and xmin
        if (this->PEz == this->mz - 1 && this->PEx == 0) {
            for (int j = 0; j < this->ny; ++j) {                
                int idx = this->Index(0, j, this->nz - 1), idxz = this->Index(0, j, this->nz - 2), idxx = this->Index(1, j, this->nz - 1);
                this->f0[idx] = 0.5*(this->f0[idxz] + this->f0[idxx]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmax and xmax
        if (this->PEz == this->mz - 1 && this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {                                        
                int idx = this->Index(this->nx - 1, j, this->nz - 1), idxz = this->Index(this->nx - 1, j, this->nz - 2), idxx = this->Index(this->nx - 2, j, this->nz - 1);
                this->f0[idx] = 0.5*(this->f0[idxz] + this->f0[idxx]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        //  Along line zmin and xmax
        if (this->PEz == 0 && this->PEx == this->mx - 1) {
            for (int j = 0; j < this->ny; ++j) {                                        
                int idx = this->Index(this->nx - 1, j, 0), idxz = this->Index(this->nx - 1, j, 1), idxx = this->Index(this->nx - 2, j, 0);
                this->f0[idx] = 0.5*(this->f0[idxz] + this->f0[idxx]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                }
            }
        }
        
        //  Along line xmin and ymin
        if (this->PEx == 0 && this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {
                int idx = this->Index(0, 0, k), idxx = this->Index(1, 0, k), idxy = this->Index(0, 1, k);
                this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmax and ymin
        if (this->PEx == this->mx - 1 && this->PEy == 0) {
            for (int k = 0; k < this->nz; ++k) {               
                int idx = this->Index(this->nx - 1, 0, k), idxx = this->Index(this->nx - 2, 0, k), idxy = this->Index(this->nx - 1, 1, k);
                this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmax and ymax
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {               
                int idx = this->Index(this->nx - 1, this->ny - 1, k), idxx = this->Index(this->nx - 2, this->ny - 1, k), idxy = this->Index(this->nx - 1, this->ny - 2, k);
                this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        //  Along line xmin and ymax
        if (this->PEx == 0 && this->PEy == this->my - 1) {
            for (int k = 0; k < this->nz; ++k) {               
                int idx = this->Index(0, this->ny - 1, k), idxx = this->Index(1, this->ny - 1, k), idxy = this->Index(0, this->ny - 2, k);
                this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                }
            }
        }
        
        //  Corner at xmin, ymin and zmin
        if (this->PEx == 0 && this->PEy == 0 && this->PEz == 0) {
            int idx = this->Index(0, 0, 0), idxx = this->Index(1, 0, 0), idxy = this->Index(0, 1, 0), idxz = this->Index(0, 0, 1);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymin and zmin
        if (this->PEx == this->mx - 1 && this->PEy == 0 && this->PEz == 0) {
            int idx = this->Index(this->nx - 1, 0, 0), idxx = this->Index(this->nx - 2, 0, 0), idxy = this->Index(this->nx - 1, 1, 0), idxz = this->Index(this->nx - 1, 0, 1);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymax and zmin
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1 && this->PEz == 0) {
            int idx = this->Index(this->nx - 1, this->ny - 1, 0), idxx = this->Index(this->nx - 2, this->ny - 1, 0), idxy = this->Index(this->nx - 1, this->ny - 2, 0), idxz = this->Index(this->nx - 1, this->ny - 1, 1);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymax and zmin
        if (this->PEx == 0 && this->PEy == this->my - 1 && this->PEz == 0) {
            int idx = this->Index(0, this->ny - 1, 0), idxx = this->Index(1, this->ny - 1, 0), idxy = this->Index(0, this->ny - 2, 0), idxz = this->Index(0, this->ny - 1, 1);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymin and zmax
        if (this->PEx == 0 && this->PEy == 0 && this->PEz == this->mz - 1) {
            int idx = this->Index(0, 0, this->nz - 1), idxx = this->Index(1, 0, this->nz - 1), idxy = this->Index(0, 1, this->nz - 1), idxz = this->Index(0, 0, this->nz - 2);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymin and zmax
        if (this->PEx == this->mx - 1 && this->PEy == 0 && this->PEz == this->mz - 1) {
            int idx = this->Index(this->nx - 1, 0, this->nz - 1), idxx = this->Index(this->nx - 2, 0, this->nz - 1), idxy = this->Index(this->nx - 1, 1, this->nz - 1), idxz = this->Index(this->nx - 1, 0, this->nz - 2);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmax, ymax and zmax
        if (this->PEx == this->mx - 1 && this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            int idx = this->Index(this->nx - 1, this->ny - 1, this->nz - 1), idxx = this->Index(this->nx - 2, this->ny - 1, this->nz - 1), idxy = this->Index(this->nx - 1, this->ny - 2, this->nz - 1), idxz = this->Index(this->nx - 1, this->ny - 1, this->nz - 2);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
        //  Corner at xmin, ymax and zmax
        if (this->PEx == 0 && this->PEy == this->my - 1 && this->PEz == this->mz - 1) {
            int idx = this->Index(0, this->ny - 1, this->nz - 1), idxx = this->Index(1, this->ny - 1, this->nz - 1), idxy = this->Index(0, this->ny - 2, this->nz - 1), idxz = this->Index(0, this->ny - 1, this->nz - 2);
            this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
            for (int c = 1; c < D3Q15<T>::nc; ++c) {
                this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
            }
        }
    }
}