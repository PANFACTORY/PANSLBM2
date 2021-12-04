//*****************************************************************************
//  Title       :   src/particle/d3q15.h
//  Author      :   Tanabe Yuta
//  Date        :   2020/10/29
//  Copyright   :   (C)2020 TanabeYuta
//*****************************************************************************

#pragma once
#include <utility>
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
#ifdef _USE_AVX_DEFINES
            this->f0 = (T*)_mm_malloc(sizeof(T)*this->nxyz, 32);
            this->f = (T*)_mm_malloc(sizeof(T)*this->nxyz*(D3Q15<T>::nc - 1), 32);
            this->fnext = (T*)_mm_malloc(sizeof(T)*this->nxyz*(D3Q15<T>::nc - 1), 32);
#else
            this->f0 = new T[this->nxyz];
            this->f = new T[this->nxyz*(D3Q15<T>::nc - 1)];
            this->fnext = new T[this->nxyz*(D3Q15<T>::nc - 1)];
#endif
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
            D3Q15<T>::LoadCxCyCzEi();
#endif
        }
        D3Q15(const D3Q15<T>& _p) = delete;
        ~D3Q15() {
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
            delete[] this->fsend_zmin;
            delete[] this->fsend_zmax;
            delete[] this->frecv_xmin;
            delete[] this->frecv_xmax;
            delete[] this->frecv_ymin;
            delete[] this->frecv_ymax;
            delete[] this->frecv_zmin;
            delete[] this->frecv_zmax;
            delete[] this->fsend_ymin_zmin;
            delete[] this->fsend_ymin_zmax;
            delete[] this->fsend_ymax_zmin;
            delete[] this->fsend_ymax_zmax;
            delete[] this->frecv_ymin_zmin;
            delete[] this->frecv_ymin_zmax;
            delete[] this->frecv_ymax_zmin;
            delete[] this->frecv_ymax_zmax; 
            delete[] this->fsend_zmin_xmin;
            delete[] this->fsend_zmin_xmax;
            delete[] this->fsend_zmax_xmin;
            delete[] this->fsend_zmax_xmax;
            delete[] this->frecv_zmin_xmin;
            delete[] this->frecv_zmin_xmax;
            delete[] this->frecv_zmax_xmin;
            delete[] this->frecv_zmax_xmax;
            delete[] this->fsend_xmin_ymin;
            delete[] this->fsend_xmin_ymax;
            delete[] this->fsend_xmax_ymin;
            delete[] this->fsend_xmax_ymax;
            delete[] this->frecv_xmin_ymin;
            delete[] this->frecv_xmin_ymax;
            delete[] this->frecv_xmax_ymin;
            delete[] this->frecv_xmax_ymax;  
        }

        int Index(int _i, int _j, int _k) const {
            int i = _i == -1 ? this->nx - 1 : (_i == this->nx ? 0 : _i);
            int j = _j == -1 ? this->ny - 1 : (_j == this->ny ? 0 : _j);
            int k = _k == -1 ? this->nz - 1 : (_k == this->nz ? 0 : _k);
            return i + this->nx*j + this->nx*this->ny*k;
        }
        static int IndexF(int _idx, int _c) {
            return (D3Q15<T>::nc - 1)*_idx + (_c - 1);
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

        void Stream(int _offset = 0);
        void iStream(int _offset = 0);

        template<class Ff>
        void BoundaryConditionAlongXFace(int _i, int _directionx, Ff _bctype);
        template<class Ff>
        void BoundaryConditionAlongYFace(int _j, int _directiony, Ff _bctype);
        template<class Ff>
        void BoundaryConditionAlongZFace(int _k, int _directionz, Ff _bctype);
        template<class Ff>
        void iBoundaryConditionAlongXFace(int _i, int _directionx, Ff _bctype);
        template<class Ff>
        void iBoundaryConditionAlongYFace(int _j, int _directiony, Ff _bctype);
        template<class Ff>
        void iBoundaryConditionAlongZFace(int _k, int _directionz, Ff _bctype);
        void SmoothCornerAlongYZ(int _j, int _k, int _directiony, int _directionz);
        void SmoothCornerAlongZX(int _k, int _i, int _directionz, int _directionx);
        void SmoothCornerAlongXY(int _i, int _j, int _directionx, int _directiony);
        void SmoothCornerAt(int _i, int _j, int _k, int _directionx, int _directiony, int _directionz);

        template<class Ff>
        void BoundaryCondition(Ff _bctype) {
            this->BoundaryConditionAlongXFace(0, -1, _bctype);              //  On xmin
            this->BoundaryConditionAlongXFace(this->lx - 1, 1, _bctype);    //  On xmax
            this->BoundaryConditionAlongYFace(0, -1, _bctype);              //  On ymin
            this->BoundaryConditionAlongYFace(this->ly - 1, 1, _bctype);    //  On ymax
            this->BoundaryConditionAlongZFace(0, -1, _bctype);              //  On zmin
            this->BoundaryConditionAlongZFace(this->lz - 1, 1, _bctype);    //  On zmax
        }
        template<class Ff>
        void iBoundaryCondition(Ff _bctype) {
            this->iBoundaryConditionAlongXFace(0, -1, _bctype);             //  On xmin
            this->iBoundaryConditionAlongXFace(this->lx - 1, 1, _bctype);   //  On xmax
            this->iBoundaryConditionAlongYFace(0, -1, _bctype);             //  On ymin
            this->iBoundaryConditionAlongYFace(this->ly - 1, 1, _bctype);   //  On ymax
            this->iBoundaryConditionAlongZFace(0, -1, _bctype);             //  On zmin
            this->iBoundaryConditionAlongZFace(this->lz - 1, 1, _bctype);   //  On zmax
        }
        void SmoothCorner() {
            this->SmoothCornerAlongYZ(0, 0, -1, -1);                        //  Along line ymin and zmin
            this->SmoothCornerAlongYZ(this->ly - 1, 0, 1, -1);              //  Along line ymax and zmin
            this->SmoothCornerAlongYZ(this->ly - 1, this->lz - 1, 1, 1);    //  Along line ymax and zmax
            this->SmoothCornerAlongYZ(0, this->lz - 1, -1, 1);              //  Along line ymin and zmax
            this->SmoothCornerAlongZX(0, 0, -1, -1);                        //  Along line zmin and xmin
            this->SmoothCornerAlongZX(this->lz - 1, 0, 1, -1);              //  Along line zmax and xmin
            this->SmoothCornerAlongZX(this->lz - 1, this->lx - 1, 1, 1);    //  Along line zmax and xmax
            this->SmoothCornerAlongZX(0, this->lx - 1, -1, 1);              //  Along line zmin and xmax
            this->SmoothCornerAlongXY(0, 0, -1, -1);                        //  Along line xmin and ymin
            this->SmoothCornerAlongXY(this->lx - 1, 0, 1, -1);              //  Along line xmax and ymin
            this->SmoothCornerAlongXY(this->lx - 1, this->ly - 1, 1, 1);    //  Along line xmax and ymax
            this->SmoothCornerAlongXY(0, this->ly - 1, -1, 1);              //  Along line xmin and ymax
            this->SmoothCornerAt(0, 0, 0, -1, -1, -1);                                  //  Corner at xmin, ymin and zmin
            this->SmoothCornerAt(this->lx - 1, 0, 0, 1, -1, -1);                        //  Corner at xmax, ymin and zmin
            this->SmoothCornerAt(this->lx - 1, this->ly - 1, 0, 1, 1, -1);              //  Corner at xmax, ymax and zmin
            this->SmoothCornerAt(0, this->ly - 1, 0, -1, 1, -1);                        //  Corner at xmin, ymax and zmin
            this->SmoothCornerAt(0, 0, this->lz - 1, -1, -1, 1);                        //  Corner at xmin, ymin and zmax
            this->SmoothCornerAt(this->lx - 1, 0, this->lz - 1, 1, -1, 1);              //  Corner at xmax, ymin and zmax
            this->SmoothCornerAt(this->lx - 1, this->ly - 1, this->lz - 1, 1, 1, 1);    //  Corner at xmax, ymax and zmax
            this->SmoothCornerAt(0, this->ly - 1, this->lz - 1, -1, 1, 1);              //  Corner at xmin, ymax and zmax
        }

        const int lx, ly, lz, PEid, mx, my, mz, PEx, PEy, PEz, nx, ny, nz, nxyz, offsetx, offsety, offsetz;
        static const int nc = 15, nd = 3, cx[nc], cy[nc], cz[nc];
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
        T *fsend_xmin, *fsend_xmax, *fsend_ymin, *fsend_ymax, *fsend_zmin, *fsend_zmax, *frecv_xmin, *frecv_xmax, *frecv_ymin, *frecv_ymax, *frecv_zmin, *frecv_zmax;
        T *fsend_ymin_zmin, *fsend_ymin_zmax, *fsend_ymax_zmin, *fsend_ymax_zmax, *frecv_ymin_zmin, *frecv_ymin_zmax, *frecv_ymax_zmin, *frecv_ymax_zmax; 
        T *fsend_zmin_xmin, *fsend_zmin_xmax, *fsend_zmax_xmin, *fsend_zmax_xmax, *frecv_zmin_xmin, *frecv_zmin_xmax, *frecv_zmax_xmin, *frecv_zmax_xmax;
        T *fsend_xmin_ymin, *fsend_xmin_ymax, *fsend_xmax_ymin, *fsend_xmax_ymax, *frecv_xmin_ymin, *frecv_xmin_ymax, *frecv_xmax_ymin, *frecv_xmax_ymax;  
        T fsend_corner[8], frecv_corner[8];
#ifdef _USE_MPI_DEFINES
        MPI_Status status[52];
        MPI_Request request[52];
#endif
        int Communicate(int _offset);
    };

    template<class T>const int D3Q15<T>::cx[D3Q15<T>::nc] = { 0, 1, 0, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, -1, -1 };
    template<class T>const int D3Q15<T>::cy[D3Q15<T>::nc] = { 0, 0, 1, 0, 0, -1, 0, 1, 1, -1, 1, -1, -1, 1, -1 };
    template<class T>const int D3Q15<T>::cz[D3Q15<T>::nc] = { 0, 0, 0, 1, 0, 0, -1, 1, 1, 1, -1, -1, -1, -1, 1 };
    template<class T>const T D3Q15<T>::ei[D3Q15<T>::nc] = { 2.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0, 1.0/72.0 };

    template<class T>
    void D3Q15<T>::Stream(int _offset) {
#ifdef _USE_MPI_DEFINES
        //  Stream of boundary points
        for (int j = 0; j < this->ny; ++j) {
            for (int k = 0; k < this->nz; ++k) {
                int idxmin = this->Index(0, j, k), idxmax = this->Index(this->nx - 1, j, k);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(0 - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index((this->nx - 1) - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }
        for (int k = 0; k < this->nz; ++k) {
            for (int i = 0; i < this->nx; ++i) {
                int idxmin = this->Index(i, 0, k), idxmax = this->Index(i, this->ny - 1, k);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(i - D3Q15<T>::cx[c], 0 - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index(i - D3Q15<T>::cx[c], (this->ny - 1) - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }
        for (int i = 0; i < this->nx; ++i) {
            for (int j = 0; j < this->ny; ++j) {
                int idxmin = this->Index(i, j, 0), idxmax = this->Index(i, j, this->nz - 1);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(i - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], 0 - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index(i - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], (this->nz - 1) - D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    int idxmin = this->Index(this->nx - 1, j, k), idxminface = this->IndexBCx(j, k);
                    this->fsend_xmin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 4)];
                    this->fsend_xmin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 8)];
                    this->fsend_xmin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 11)];
                    this->fsend_xmin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 13)];
                    this->fsend_xmin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 14)];

                    //  Face on xmax
                    int idxmax = this->Index(0, j, k), idxmaxface = this->IndexBCx(j, k);
                    this->fsend_xmax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 1)];
                    this->fsend_xmax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 7)];
                    this->fsend_xmax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 9)];
                    this->fsend_xmax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 10)];
                    this->fsend_xmax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 12)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    int idxmin = this->Index(i, this->ny - 1, k), idxminface = this->IndexBCy(k, i);
                    this->fsend_ymin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 5)];
                    this->fsend_ymin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 9)];
                    this->fsend_ymin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 11)];
                    this->fsend_ymin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 12)];
                    this->fsend_ymin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 14)];

                    //  Face on ymax
                    int idxmax = this->Index(i, 0, k), idxmaxface = this->IndexBCy(k, i);
                    this->fsend_ymax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 2)];
                    this->fsend_ymax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 7)];
                    this->fsend_ymax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 8)];
                    this->fsend_ymax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 10)];
                    this->fsend_ymax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 13)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    int idxmin = this->Index(i, j, this->nz - 1), idxminface = this->IndexBCz(i, j);
                    this->fsend_zmin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 6)];
                    this->fsend_zmin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 10)];
                    this->fsend_zmin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 11)];
                    this->fsend_zmin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 12)];
                    this->fsend_zmin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 13)];

                    //  Face on zmax
                    int idxmax = this->Index(i, j, 0), idxmaxface = this->IndexBCz(i, j);
                    this->fsend_zmax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 3)];
                    this->fsend_zmax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 7)];
                    this->fsend_zmax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 8)];
                    this->fsend_zmax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 9)];
                    this->fsend_zmax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 14)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                int idxminmin = this->Index(i, this->ny - 1, this->nz - 1);
                this->fsend_ymin_zmin[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 11)];
                this->fsend_ymin_zmin[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 12)];

                //  Edge on ymin and zmax
                int idxminmax = this->Index(i, this->ny - 1, 0);
                this->fsend_ymin_zmax[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 9)];
                this->fsend_ymin_zmax[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 14)];

                //  Edge on ymax and zmin
                int idxmaxmin = this->Index(i, 0, this->nz - 1);
                this->fsend_ymax_zmin[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 10)];
                this->fsend_ymax_zmin[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 13)];

                //  Edge on ymax and zmax
                int idxmaxmax = this->Index(i, 0, 0);
                this->fsend_ymax_zmax[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)];
                this->fsend_ymax_zmax[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 8)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                int idxminmin = this->Index(this->nx - 1, j, this->nz - 1);
                this->fsend_zmin_xmin[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 11)];
                this->fsend_zmin_xmin[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 13)];

                //  Edge on zmin and xmax
                int idxminmax = this->Index(0, j, this->nz - 1);
                this->fsend_zmin_xmax[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 10)];
                this->fsend_zmin_xmax[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 12)];

                //  Edge on zmax and xmin
                int idxmaxmin = this->Index(this->nx - 1, j, 0);
                this->fsend_zmax_xmin[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 8)];
                this->fsend_zmax_xmin[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 14)];

                //  Edge on zmax and xmax
                int idxmaxmax = this->Index(0, j, 0);
                this->fsend_zmax_xmax[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)];
                this->fsend_zmax_xmax[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 9)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                int idxminmin = this->Index(this->nx - 1, this->ny - 1, k);
                this->fsend_xmin_ymin[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 11)];
                this->fsend_xmin_ymin[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 14)];

                //  Edge on xmin and ymax
                int idxminmax = this->Index(this->nx - 1, 0, k);
                this->fsend_xmin_ymax[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 8)];
                this->fsend_xmin_ymax[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 13)];

                //  Edge on xmax and ymin
                int idxmaxmin = this->Index(0, this->ny - 1, k);
                this->fsend_xmax_ymin[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 9)];
                this->fsend_xmax_ymin[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 12)];

                //  Edge on xmax and ymax
                int idxmaxmax = this->Index(0, 0, k);
                this->fsend_xmax_ymax[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)];
                this->fsend_xmax_ymax[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 10)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend_corner[0] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)];   //  Corner at xmin, ymin and zmin
            this->fsend_corner[1] = this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)];              //  Corner at xmax, ymin and zmin
            this->fsend_corner[2] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)];              //  Corner at xmin, ymax and zmin 
            this->fsend_corner[3] = this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)];                         //  Corner at xmax, ymax and zmin
            this->fsend_corner[4] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)];              //  Corner at xmin, ymin and zmax
            this->fsend_corner[5] = this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)];                          //  Corner at xmax, ymin and zmax
            this->fsend_corner[6] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)];                          //  Corner at xmin, ymax and zmax
            this->fsend_corner[7] = this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)];                                     //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        int neib = this->Communicate(_offset);
        
        //  Stream of inner points
#pragma omp parallel for
        for (int k = 1; k < this->nz - 1; ++k) {
            for (int j = 1; j < this->ny - 1; ++j) {
                for (int i = 1; i < this->nx - 1; ++i) {
                    int idx = this->Index(i, j, k);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idxstream = this->Index(i - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        } 
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    int idxmin = this->Index(0, j, k), idxminface = this->IndexBCx(j, k);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 1)] = this->frecv_xmin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 7)] = this->frecv_xmin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 9)] = this->frecv_xmin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 10)] = this->frecv_xmin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 12)] = this->frecv_xmin[idxminface*5 + 4];

                    //  Face on xmax
                    int idxmax = this->Index(this->nx - 1, j, k), idxmaxface = this->IndexBCx(j, k);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 4)] = this->frecv_xmax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 8)] = this->frecv_xmax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 11)] = this->frecv_xmax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 13)] = this->frecv_xmax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 14)] = this->frecv_xmax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    int idxmin = this->Index(i, 0, k), idxminface = this->IndexBCy(k, i);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 2)] = this->frecv_ymin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 7)] = this->frecv_ymin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 8)] = this->frecv_ymin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 10)] = this->frecv_ymin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 13)] = this->frecv_ymin[idxminface*5 + 4];

                    //  Face on ymax
                    int idxmax = this->Index(i, this->ny - 1, k), idxmaxface = this->IndexBCy(k, i);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 5)] = this->frecv_ymax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 9)] = this->frecv_ymax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 11)] = this->frecv_ymax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 12)] = this->frecv_ymax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 14)] = this->frecv_ymax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    int idxmin = this->Index(i, j, 0), idxminface = this->IndexBCz(i, j);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 3)] = this->frecv_zmin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 7)] = this->frecv_zmin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 8)] = this->frecv_zmin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 9)] = this->frecv_zmin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 14)] = this->frecv_zmin[idxminface*5 + 4];

                    //  Face on zmax
                    int idxmax = this->Index(i, j, this->nz - 1), idxmaxface = this->IndexBCz(i, j);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 6)] = this->frecv_zmax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 10)] = this->frecv_zmax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 11)] = this->frecv_zmax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 12)] = this->frecv_zmax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 13)] = this->frecv_zmax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                int idxminmin = this->Index(i, 0, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 7)] = this->frecv_ymin_zmin[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 8)] = this->frecv_ymin_zmin[i*2 + 1];

                //  Edge on ymin and zmax
                int idxminmax = this->Index(i, 0, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 10)] = this->frecv_ymin_zmax[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 13)] = this->frecv_ymin_zmax[i*2 + 1];

                //  Edge on ymax and zmin
                int idxmaxmin = this->Index(i, this->ny - 1, 0);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 9)] = this->frecv_ymax_zmin[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 14)] = this->frecv_ymax_zmin[i*2 + 1];

                //  Edge on ymax and zmax
                int idxmaxmax = this->Index(i, this->ny - 1, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)] = this->frecv_ymax_zmax[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 12)] = this->frecv_ymax_zmax[i*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                int idxminmin = this->Index(0, j, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 7)] = this->frecv_zmin_xmin[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 9)] = this->frecv_zmin_xmin[j*2 + 1];

                //  Edge on zmin and xmax
                int idxminmax = this->Index(this->nx - 1, j, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 8)] = this->frecv_zmin_xmax[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 14)] = this->frecv_zmin_xmax[j*2 + 1];

                //  Edge on zmax and xmin
                int idxmaxmin = this->Index(0, j, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 10)] = this->frecv_zmax_xmin[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 12)] = this->frecv_zmax_xmin[j*2 + 1];

                //  Edge on zmax and xmax
                int idxmaxmax = this->Index(this->nx - 1, j, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)] = this->frecv_zmax_xmax[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 13)] = this->frecv_zmax_xmax[j*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                int idxminmin = this->Index(0, 0, k);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 7)] = this->frecv_xmin_ymin[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 10)] = this->frecv_xmin_ymin[k*2 + 1];

                //  Edge on xmin and ymax
                int idxminmax = this->Index(0, this->ny - 1, k);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 9)] = this->frecv_xmin_ymax[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 12)] = this->frecv_xmin_ymax[k*2 + 1];

                //  Edge on xmax and ymin
                int idxmaxmin = this->Index(this->nx - 1, 0, k);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 8)] = this->frecv_xmax_ymin[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 13)] = this->frecv_xmax_ymin[k*2 + 1];

                //  Edge on xmax and ymax
                int idxmaxmax = this->Index(this->nx - 1, this->ny - 1, k);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)] = this->frecv_xmax_ymax[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 14)] = this->frecv_xmax_ymax[k*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, 0), 7)] = this->frecv_corner[0];                                     //  Corner at xmin, ymin and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 8)] = this->frecv_corner[1];                          //  Corner at xmax, ymin and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 9)] = this->frecv_corner[2];                          //  Corner at xmin, ymax and zmin 
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 14)] = this->frecv_corner[3];              //  Corner at xmax, ymax and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 10)] = this->frecv_corner[4];                         //  Corner at xmin, ymin and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 13)] = this->frecv_corner[5];              //  Corner at xmax, ymin and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 12)] = this->frecv_corner[6];              //  Corner at xmin, ymax and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 11)] = this->frecv_corner[7];   //  Corner at xmax, ymax and zmax
        }
#else
#pragma omp parallel for
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    int idx = this->Index(i, j, k);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idxstream = this->Index(i - D3Q15<T>::cx[c], j - D3Q15<T>::cy[c], k - D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        }   
#endif
        //  Swap
        std::swap(this->f, this->fnext);
    }

    template<class T>
    void D3Q15<T>::iStream(int _offset) {        
#ifdef _USE_MPI_DEFINES
        //  Stream of boundary points
        for (int j = 0; j < this->ny; ++j) {
            for (int k = 0; k < this->nz; ++k) {
                int idxmin = this->Index(0, j, k), idxmax = this->Index(this->nx - 1, j, k);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(0 + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index((this->nx - 1) + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }
        for (int k = 0; k < this->nz; ++k) {
            for (int i = 0; i < this->nx; ++i) {
                int idxmin = this->Index(i, 0, k), idxmax = this->Index(i, this->ny - 1, k);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(i + D3Q15<T>::cx[c], 0 + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index(i + D3Q15<T>::cx[c], (this->ny - 1) + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }
        for (int i = 0; i < this->nx; ++i) {
            for (int j = 0; j < this->ny; ++j) {
                int idxmin = this->Index(i, j, 0), idxmax = this->Index(i, j, this->nz - 1);
                for (int c = 1; c < D3Q15<T>::nc; ++c) {
                    int idxminstream = this->Index(i + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], 0 + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmin, c)] = this->f[D3Q15<T>::IndexF(idxminstream, c)];
                    int idxmaxstream = this->Index(i + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], (this->nz - 1) + D3Q15<T>::cz[c]);
                    this->fnext[D3Q15<T>::IndexF(idxmax, c)] = this->f[D3Q15<T>::IndexF(idxmaxstream, c)];
                }
            }
        }

        //  Copy from f to fsend along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    int idxmin = this->Index(this->nx - 1, j, k), idxminface = this->IndexBCx(j, k);
                    this->fsend_xmin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 1)];
                    this->fsend_xmin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 7)];
                    this->fsend_xmin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 9)];
                    this->fsend_xmin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 10)];
                    this->fsend_xmin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 12)];

                    //  Face on xmax
                    int idxmax = this->Index(0, j, k), idxmaxface = this->IndexBCx(j, k);
                    this->fsend_xmax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 4)];
                    this->fsend_xmax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 8)];
                    this->fsend_xmax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 11)];
                    this->fsend_xmax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 13)];
                    this->fsend_xmax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 14)];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    int idxmin = this->Index(i, this->ny - 1, k), idxminface = this->IndexBCy(k, i);
                    this->fsend_ymin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 2)];
                    this->fsend_ymin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 7)];
                    this->fsend_ymin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 8)];
                    this->fsend_ymin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 10)];
                    this->fsend_ymin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 13)];

                    //  Face on ymax
                    int idxmax = this->Index(i, 0, k), idxmaxface = this->IndexBCy(k, i);
                    this->fsend_ymax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 5)];
                    this->fsend_ymax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 9)];
                    this->fsend_ymax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 11)];
                    this->fsend_ymax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 12)];
                    this->fsend_ymax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 14)];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    int idxmin = this->Index(i, j, this->nz - 1), idxminface = this->IndexBCz(i, j);
                    this->fsend_zmin[idxminface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmin, 3)];
                    this->fsend_zmin[idxminface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmin, 7)];
                    this->fsend_zmin[idxminface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmin, 8)];
                    this->fsend_zmin[idxminface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmin, 9)];
                    this->fsend_zmin[idxminface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmin, 14)];

                    //  Face on zmax
                    int idxmax = this->Index(i, j, 0), idxmaxface = this->IndexBCz(i, j);
                    this->fsend_zmax[idxmaxface*5 + 0] = this->fnext[D3Q15<T>::IndexF(idxmax, 6)];
                    this->fsend_zmax[idxmaxface*5 + 1] = this->fnext[D3Q15<T>::IndexF(idxmax, 10)];
                    this->fsend_zmax[idxmaxface*5 + 2] = this->fnext[D3Q15<T>::IndexF(idxmax, 11)];
                    this->fsend_zmax[idxmaxface*5 + 3] = this->fnext[D3Q15<T>::IndexF(idxmax, 12)];
                    this->fsend_zmax[idxmaxface*5 + 4] = this->fnext[D3Q15<T>::IndexF(idxmax, 13)];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                int idxminmin = this->Index(i, this->ny - 1, this->nz - 1);
                this->fsend_ymin_zmin[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 7)];
                this->fsend_ymin_zmin[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 8)];

                //  Edge on ymin and zmax
                int idxminmax = this->Index(i, this->ny - 1, 0);
                this->fsend_ymin_zmax[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 10)];
                this->fsend_ymin_zmax[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 13)];

                //  Edge on ymax and zmin
                int idxmaxmin = this->Index(i, 0, this->nz - 1);
                this->fsend_ymax_zmin[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 9)];
                this->fsend_ymax_zmin[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 14)];

                //  Edge on ymax and zmax
                int idxmaxmax = this->Index(i, 0, 0);
                this->fsend_ymax_zmax[i*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)];
                this->fsend_ymax_zmax[i*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 12)];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                int idxminmin = this->Index(this->nx - 1, j, this->nz - 1);
                this->fsend_zmin_xmin[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 7)];
                this->fsend_zmin_xmin[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 9)];

                //  Edge on zmin and xmax
                int idxminmax = this->Index(0, j, this->nz - 1);
                this->fsend_zmin_xmax[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 8)];
                this->fsend_zmin_xmax[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 14)];

                //  Edge on zmax and xmin
                int idxmaxmin = this->Index(this->nx - 1, j, 0);
                this->fsend_zmax_xmin[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 10)];
                this->fsend_zmax_xmin[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 12)];

                //  Edge on zmax and xmax
                int idxmaxmax = this->Index(0, j, 0);
                this->fsend_zmax_xmax[j*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)];
                this->fsend_zmax_xmax[j*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 13)];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                int idxminmin = this->Index(this->nx - 1, this->ny - 1, k);
                this->fsend_xmin_ymin[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmin, 7)];
                this->fsend_xmin_ymin[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmin, 10)];

                //  Edge on xmin and ymax
                int idxminmax = this->Index(this->nx - 1, 0, k);
                this->fsend_xmin_ymax[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxminmax, 9)];
                this->fsend_xmin_ymax[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxminmax, 12)];
                
                //  Edge on xmax and ymin
                int idxmaxmin = this->Index(0, this->ny - 1, k);
                this->fsend_xmax_ymin[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 8)];
                this->fsend_xmax_ymin[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmin, 13)];

                //  Edge on xmax and ymax
                int idxmaxmax = this->Index(0, 0, k);
                this->fsend_xmax_ymax[k*2 + 0] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 11)];
                this->fsend_xmax_ymax[k*2 + 1] = this->fnext[D3Q15<T>::IndexF(idxmaxmax, 14)];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fsend_corner[0] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 7)];    //  Corner at xmin, ymin and zmin
            this->fsend_corner[1] = this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 8)];               //  Corner at xmax, ymin and zmin
            this->fsend_corner[2] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 9)];               //  Corner at xmin, ymax and zmin 
            this->fsend_corner[3] = this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 14)];                         //  Corner at xmax, ymax and zmin
            this->fsend_corner[4] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 10)];              //  Corner at xmin, ymin and zmax
            this->fsend_corner[5] = this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 13)];                         //  Corner at xmax, ymin and zmax
            this->fsend_corner[6] = this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 12)];                         //  Corner at xmin, ymax and zmax
            this->fsend_corner[7] = this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, 0), 11)];                                    //  Corner at xmax, ymax and zmax
        }

        //  Communicate with other PE
        int neib = this->Communicate(_offset);

        //  Stream of inner points
#pragma omp parallel for
        for (int k = 1; k < this->nz - 1; ++k) {
            for (int j = 1; j < this->ny - 1; ++j) {
                for (int i = 1; i < this->nx - 1; ++i) {
                    int idx = this->Index(i, j, k);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idxstream = this->Index(i + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        }
        if (neib > 0) {
            MPI_Waitall(neib, this->request, this->status);
        }

        //  Copy to f from frecv along edge or at corner
        if (this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    //  Face on xmin
                    int idxmin = this->Index(0, j, k), idxminface = this->IndexBCx(j, k);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 4)] = this->frecv_xmin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 8)] = this->frecv_xmin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 11)] = this->frecv_xmin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 13)] = this->frecv_xmin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 14)] = this->frecv_xmin[idxminface*5 + 4];

                    //  Face on xmax
                    int idxmax = this->Index(this->nx - 1, j, k), idxmaxface = this->IndexBCx(j, k);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 1)] = this->frecv_xmax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 7)] = this->frecv_xmax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 9)] = this->frecv_xmax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 10)] = this->frecv_xmax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 12)] = this->frecv_xmax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    //  Face on ymin
                    int idxmin = this->Index(i, 0, k), idxminface = this->IndexBCy(k, i);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 5)] = this->frecv_ymin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 9)] = this->frecv_ymin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 11)] = this->frecv_ymin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 12)] = this->frecv_ymin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 14)] = this->frecv_ymin[idxminface*5 + 4];

                    //  Face on ymax
                    int idxmax = this->Index(i, this->ny - 1, k), idxmaxface = this->IndexBCy(k, i);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 2)] = this->frecv_ymax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 7)] = this->frecv_ymax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 8)] = this->frecv_ymax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 10)] = this->frecv_ymax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 13)] = this->frecv_ymax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    //  Face on zmin
                    int idxmin = this->Index(i, j, 0), idxminface = this->IndexBCz(i, j);
                    this->fnext[D3Q15<T>::IndexF(idxmin, 6)] = this->frecv_zmin[idxminface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 10)] = this->frecv_zmin[idxminface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 11)] = this->frecv_zmin[idxminface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 12)] = this->frecv_zmin[idxminface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmin, 13)] = this->frecv_zmin[idxminface*5 + 4];

                    //  Face on zmax
                    int idxmax = this->Index(i, j, this->nz - 1), idxmaxface = this->IndexBCz(i, j);
                    this->fnext[D3Q15<T>::IndexF(idxmax, 3)] = this->frecv_zmax[idxmaxface*5 + 0];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 7)] = this->frecv_zmax[idxmaxface*5 + 1];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 8)] = this->frecv_zmax[idxmaxface*5 + 2];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 9)] = this->frecv_zmax[idxmaxface*5 + 3];
                    this->fnext[D3Q15<T>::IndexF(idxmax, 14)] = this->frecv_zmax[idxmaxface*5 + 4];
                }
            }
        }
        if (this->my != 1 || this->mz != 1) {
            for (int i = 0; i < this->nx; ++i) {
                //  Edge on ymin and zmin
                int idxminmin = this->Index(i, 0, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 11)] = this->frecv_ymin_zmin[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 12)] = this->frecv_ymin_zmin[i*2 + 1];

                //  Edge on ymin and zmax
                int idxminmax = this->Index(i, 0, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 9)] = this->frecv_ymin_zmax[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 14)] = this->frecv_ymin_zmax[i*2 + 1];

                //  Edge on ymax and zmin
                int idxmaxmin = this->Index(i, this->ny - 1, 0);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 10)] = this->frecv_ymax_zmin[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 13)] = this->frecv_ymax_zmin[i*2 + 1];

                //  Edge on ymax and zmax
                int idxmaxmax = this->Index(i, this->ny - 1, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)] = this->frecv_ymax_zmax[i*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 8)] = this->frecv_ymax_zmax[i*2 + 1];
            }
        }
        if (this->mz != 1 || this->mx != 1) {
            for (int j = 0; j < this->ny; ++j) {
                //  Edge on zmin and xmin
                int idxminmin = this->Index(0, j, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 11)] = this->frecv_zmin_xmin[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 13)] = this->frecv_zmin_xmin[j*2 + 1];

                //  Edge on zmin and xmax
                int idxminmax = this->Index(this->nx - 1, j, 0);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 10)] = this->frecv_zmin_xmax[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 12)] = this->frecv_zmin_xmax[j*2 + 1];

                //  Edge on zmax and xmin
                int idxmaxmin = this->Index(0, j, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 8)] = this->frecv_zmax_xmin[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 14)] = this->frecv_zmax_xmin[j*2 + 1];

                //  Edge on zmax and xmax
                int idxmaxmax = this->Index(this->nx - 1, j, this->nz - 1);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)] = this->frecv_zmax_xmax[j*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 9)] = this->frecv_zmax_xmax[j*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1) {
            for (int k = 0; k < this->nz; ++k) {
                //  Edge on xmin and ymin
                int idxminmin = this->Index(0, 0, k);
                this->fnext[D3Q15<T>::IndexF(idxminmin, 11)] = this->frecv_xmin_ymin[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmin, 14)] = this->frecv_xmin_ymin[k*2 + 1];

                //  Edge on xmin and ymax
                int idxminmax = this->Index(0, this->ny - 1, k);
                this->fnext[D3Q15<T>::IndexF(idxminmax, 8)] = this->frecv_xmin_ymax[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxminmax, 13)] = this->frecv_xmin_ymax[k*2 + 1];

                //  Edge on xmax and ymin
                int idxmaxmin = this->Index(this->nx - 1, 0, k);
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 9)] = this->frecv_xmax_ymin[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmin, 12)] = this->frecv_xmax_ymin[k*2 + 1];

                //  Edge on xmax and ymax
                int idxmaxmax = this->Index(this->nx - 1, this->ny - 1, k);
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 7)] = this->frecv_xmax_ymax[k*2 + 0];
                this->fnext[D3Q15<T>::IndexF(idxmaxmax, 10)] = this->frecv_xmax_ymax[k*2 + 1];
            }
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {
            this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, 0), 11)] = this->frecv_corner[0];                                    //  Corner at xmin, ymin and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, 0), 12)] = this->frecv_corner[1];                         //  Corner at xmax, ymin and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, 0), 13)] = this->frecv_corner[2];                         //  Corner at xmin, ymax and zmin 
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, 0), 10)] = this->frecv_corner[3];              //  Corner at xmax, ymax and zmin
            this->fnext[D3Q15<T>::IndexF(this->Index(0, 0, this->nz - 1), 14)] = this->frecv_corner[4];                         //  Corner at xmin, ymin and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, 0, this->nz - 1), 9)] = this->frecv_corner[5];               //  Corner at xmax, ymin and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(0, this->ny - 1, this->nz - 1), 8)] = this->frecv_corner[6];               //  Corner at xmin, ymax and zmax
            this->fnext[D3Q15<T>::IndexF(this->Index(this->nx - 1, this->ny - 1, this->nz - 1), 7)] = this->frecv_corner[7];    //  Corner at xmax, ymax and zmax
        }
#else
#pragma omp parallel for
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    int idx = this->Index(i, j, k);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        int idxstream = this->Index(i + D3Q15<T>::cx[c], j + D3Q15<T>::cy[c], k + D3Q15<T>::cz[c]);
                        this->fnext[D3Q15<T>::IndexF(idx, c)] = this->f[D3Q15<T>::IndexF(idxstream, c)];
                    }
                }
            }
        }
#endif
        //  Swap
        std::swap(this->f, this->fnext);
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::BoundaryConditionAlongXFace(int _i, int _directionx, Ff _bctype) {
        int i = _i - this->offsetx;
        if (0 <= i && i < this->nx) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directionx == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        }
                    } else if (_directionx == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
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
        }
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::BoundaryConditionAlongYFace(int _j, int _directiony, Ff _bctype) {
        int j = _j - this->offsety;
        if (0 <= j && j < this->ny) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directiony == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        }
                    } else if (_directiony == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        }
                    }
                }
            }
        }
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::BoundaryConditionAlongZFace(int _k, int _directionz, Ff _bctype) {
        int k = _k - this->offsetz;
        if (0 <= k && k < this->nz) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directionz == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        }
                    } else if (_directionz == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        } else if (bctype == MIRROR) {
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
    }
        
    template<class T>
    template<class Ff>
    void D3Q15<T>::iBoundaryConditionAlongXFace(int _i, int _directionx, Ff _bctype) {
        int i = _i - this->offsetx;
        if (0 <= i && i < this->nx) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directionx == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 4)] = this->f[D3Q15<T>::IndexF(idx, 1)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        }
                    } else if (_directionx == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 1)] = this->f[D3Q15<T>::IndexF(idx, 4)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        }
                    }
                }
            }
        }
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::iBoundaryConditionAlongYFace(int _j, int _directiony, Ff _bctype) {
        int j = _j - this->offsety;
        if (0 <= j && j < this->ny) {
            for (int k = 0; k < this->nz; ++k) {
                for (int i = 0; i < this->nx; ++i) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directiony == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 5)] = this->f[D3Q15<T>::IndexF(idx, 2)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        }
                    } else if (_directiony == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 2)] = this->f[D3Q15<T>::IndexF(idx, 5)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                        }
                    }
                }
            }
        }
    }

    template<class T>
    template<class Ff>
    void D3Q15<T>::iBoundaryConditionAlongZFace(int _k, int _directionz, Ff _bctype) {
        int k = _k - this->offsetz;
        if (0 <= k && k < this->nz) {
            for (int i = 0; i < this->nx; ++i) {
                for (int j = 0; j < this->ny; ++j) {
                    int idx = this->Index(i, j, k);
                    int bctype = _bctype(i + this->offsetx, j + this->offsety, k + this->offsetz);
                    if (_directionz == -1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                        } else if (bctype == MIRROR) {
                            this->f[D3Q15<T>::IndexF(idx, 6)] = this->f[D3Q15<T>::IndexF(idx, 3)];
                            this->f[D3Q15<T>::IndexF(idx, 10)] = this->f[D3Q15<T>::IndexF(idx, 7)];
                            this->f[D3Q15<T>::IndexF(idx, 11)] = this->f[D3Q15<T>::IndexF(idx, 14)];
                            this->f[D3Q15<T>::IndexF(idx, 12)] = this->f[D3Q15<T>::IndexF(idx, 9)];
                            this->f[D3Q15<T>::IndexF(idx, 13)] = this->f[D3Q15<T>::IndexF(idx, 8)];
                        }
                    } else if (_directionz == 1) {
                        if (bctype == BARRIER) {
                            this->f[D3Q15<T>::IndexF(idx, 3)] = this->f[D3Q15<T>::IndexF(idx, 6)];
                            this->f[D3Q15<T>::IndexF(idx, 7)] = this->f[D3Q15<T>::IndexF(idx, 11)];
                            this->f[D3Q15<T>::IndexF(idx, 8)] = this->f[D3Q15<T>::IndexF(idx, 12)];
                            this->f[D3Q15<T>::IndexF(idx, 9)] = this->f[D3Q15<T>::IndexF(idx, 13)];
                            this->f[D3Q15<T>::IndexF(idx, 14)] = this->f[D3Q15<T>::IndexF(idx, 10)];
                        } else if (bctype == MIRROR) {
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
    }

    template<class T>
    void D3Q15<T>::SmoothCornerAlongYZ(int _j, int _k, int _directiony, int _directionz) {
        int j = _j - this->offsety, k = _k - this->offsetz;
        if (0 <= j && j < this->ny) {
            if (0 <= k && k < this->nz) {
                for (int i = 0; i < this->nx; ++i) {
                    int idx = this->Index(i, j, k), idxy = this->Index(i, j - _directiony, k), idxz = this->Index(i, j, k - _directionz);
                    this->f0[idx] = 0.5*(this->f0[idxy] + this->f0[idxz]);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)]);
                    }
                }
            }    
        }
    }

    template<class T>
    void D3Q15<T>::SmoothCornerAlongZX(int _k, int _i, int _directionz, int _directionx) {
        int k = _k - this->offsetz, i = _i - this->offsetx;
        if (0 <= k && k < this->nz) {
            if (0 <= i && i < this->nx) {
                for (int j = 0; j < this->ny; ++j) {
                    int idx = this->Index(i, j, k), idxz = this->Index(i, j, k - _directionz), idxx = this->Index(i - _directionx, j, k);
                    this->f0[idx] = 0.5*(this->f0[idxz] + this->f0[idxx]);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxz, c)] + this->f[D3Q15<T>::IndexF(idxx, c)]);
                    }
                }
            }
        }
    }

    template<class T>
    void D3Q15<T>::SmoothCornerAlongXY(int _i, int _j, int _directionx, int _directiony) {
        int i = _i - this->offsetx, j = _j - this->offsety;
        if (0 <= i && i < this->nx) {
            if (0 <= j && j < this->ny) {
                for (int k = 0; k < this->nz; ++k) {
                    int idx = this->Index(i, j, k), idxx = this->Index(i - _directionx, j, k), idxy = this->Index(i, j - _directiony, k);
                    this->f0[idx] = 0.5*(this->f0[idxx] + this->f0[idxy]);
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        this->f[D3Q15<T>::IndexF(idx, c)] = 0.5*(this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)]);
                    }
                }
            }
        }
    }

    template<class T>
    void D3Q15<T>::SmoothCornerAt(int _i, int _j, int _k, int _directionx, int _directiony, int _directionz) {
        int i = _i - this->offsetx, j = _j - this->offsety, k = _k - this->offsetz;
        if (0 <= i && i < this->nx) {
            if (0 <= j && j < this->ny) {
                if (0 <= k && k < this->nz) {
                    int idx = this->Index(i, j, k), idxx = this->Index(i - _directionx, j, k), idxy = this->Index(i, j - _directiony, k), idxz = this->Index(i, j, k - _directionz);
                    this->f0[idx] = (this->f0[idxx] + this->f0[idxy] + this->f0[idxz])/3.0;
                    for (int c = 1; c < D3Q15<T>::nc; ++c) {
                        this->f[D3Q15<T>::IndexF(idx, c)] = (this->f[D3Q15<T>::IndexF(idxx, c)] + this->f[D3Q15<T>::IndexF(idxy, c)] + this->f[D3Q15<T>::IndexF(idxz, c)])/3.0;
                    }
                }
            }
        }
    }

    template<class T>
    int D3Q15<T>::Communicate(int _offset) {
        int neib = 0;
#ifdef _USE_MPI_DEFINES
        if (this->mx != 1) {
            //  To xmin
            MPI_Isend(this->fsend_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0 + _offset, MPI_COMM_WORLD, &this->request[neib++]);

            //  To xmax
            MPI_Isend(this->fsend_xmax, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 1 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin, this->ny*this->nz*5, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 1 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1) {
            //  To ymin
            MPI_Isend(this->fsend_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 2 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 2 + _offset, MPI_COMM_WORLD, &this->request[neib++]);

            //  To ymax
            MPI_Isend(this->fsend_ymax, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz), 3 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin, this->nz*this->nx*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz), 3 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1) {
            //  To zmin
            MPI_Isend(this->fsend_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 4 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 4 + _offset, MPI_COMM_WORLD, &this->request[neib++]);

            //  To zmax
            MPI_Isend(this->fsend_zmax, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz + 1), 5 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin, this->nx*this->ny*5, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy, this->PEz - 1), 5 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->my != 1 || this->mz != 1) {
            //  To ymin and zmin
            MPI_Isend(this->fsend_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 6 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 6 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymin and zmax
            MPI_Isend(this->fsend_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 7 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 7 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmin
            MPI_Isend(this->fsend_ymax_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz - 1), 8 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz + 1), 8 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To ymax and zmax
            MPI_Isend(this->fsend_ymax_zmax, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy + 1, this->PEz + 1), 9 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_ymin_zmin, this->nx*2, MPI_DOUBLE, this->IndexPE(this->PEx, this->PEy - 1, this->PEz - 1), 9 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mz != 1 || this->mx != 1) {
            //  To zmin and xmin
            MPI_Isend(this->fsend_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 10 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 10 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                
            //  To zmin and xmax
            MPI_Isend(this->fsend_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 11 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 11 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                 
            //  To zmax and xmin
            MPI_Isend(this->fsend_zmax_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz + 1), 12 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz - 1), 12 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To zmax and xmax
            MPI_Isend(this->fsend_zmax_xmax, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz + 1), 13 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_zmin_xmin, this->ny*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz - 1), 13 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1) {
            //  To xmin and ymin
            MPI_Isend(this->fsend_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 14 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 14 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                             
            //  To xmin and ymax
            MPI_Isend(this->fsend_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 15 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 15 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymin
            MPI_Isend(this->fsend_xmax_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz), 16 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz), 16 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
                           
            //  To xmax and ymax
            MPI_Isend(this->fsend_xmax_ymax, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz), 17 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(this->frecv_xmin_ymin, this->nz*2, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz), 17 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
        }
        if (this->mx != 1 || this->my != 1 || this->mz != 1) {  //ここにミス（D2Q9と比較せよ）
            MPI_Isend(&this->fsend_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 18 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 18 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmin
            MPI_Isend(&this->fsend_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 19 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 19 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmin
            MPI_Isend(&this->fsend_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 20 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 20 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmin
            MPI_Isend(&this->fsend_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 21 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 21 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmin
            MPI_Isend(&this->fsend_corner[4], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz + 1), 22 + _offset, MPI_COMM_WORLD, &this->request[neib++]);
            MPI_Irecv(&this->frecv_corner[3], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz - 1), 22 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymin and zmax
            MPI_Isend(&this->fsend_corner[5], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz + 1), 23 + _offset, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[2], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz - 1), 23 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymin and zmax
            MPI_Isend(&this->fsend_corner[6], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy + 1, this->PEz + 1), 24 + _offset, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[1], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy - 1, this->PEz - 1), 24 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmin, ymax and zmax
            MPI_Isend(&this->fsend_corner[7], 1, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy + 1, this->PEz + 1), 25 + _offset, MPI_COMM_WORLD, &this->request[neib++]); 
            MPI_Irecv(&this->frecv_corner[0], 1, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy - 1, this->PEz - 1), 25 + _offset, MPI_COMM_WORLD, &this->request[neib++]);   //  Corner at xmax, ymax and zmax
        }
#endif
        return neib;
    }

#ifdef _USE_AVX_DEFINES
    template<>__m256d D3Q15<double>::__cx[D3Q15<double>::nc] = { 0 };
    template<>__m256d D3Q15<double>::__cy[D3Q15<double>::nc] = { 0 };
    template<>__m256d D3Q15<double>::__cz[D3Q15<double>::nc] = { 0 };
    template<>__m256d D3Q15<double>::__ei[D3Q15<double>::nc] = { 0 };

    template<>
    void D3Q15<double>::LoadCxCyCzEi() {
        for (int c = 0; c < D3Q15<double>::nc; ++c) {
            D3Q15<double>::__cx[c] = _mm256_set1_pd((double)D3Q15<double>::cx[c]);
            D3Q15<double>::__cy[c] = _mm256_set1_pd((double)D3Q15<double>::cy[c]);
            D3Q15<double>::__cz[c] = _mm256_set1_pd((double)D3Q15<double>::cz[c]);
            D3Q15<double>::__ei[c] = _mm256_set1_pd((double)D3Q15<double>::ei[c]);
        }
    }

    template<>
    template<>
    void D3Q15<double>::LoadF<__m256d>(int _idx, __m256d *__f) {
        const int offsetf = D3Q15<double>::IndexF(_idx, 1);
        __m256d load0 = _mm256_load_pd(&this->f[offsetf +  0]);   //  f 4(0) f 3(0) f 2(0) f 1(0)
        __m256d load1 = _mm256_load_pd(&this->f[offsetf +  4]);   //  f 8(0) f 7(0) f 6(0) f 5(0)
        __m256d load2 = _mm256_load_pd(&this->f[offsetf +  8]);   //  f12(0) f11(0) f10(0) f 9(0)
        __m256d load3 = _mm256_load_pd(&this->f[offsetf + 12]);   //  f 2(1) f 1(1) f14(0) f13(0)
        __m256d load4 = _mm256_load_pd(&this->f[offsetf + 16]);   //  f 6(1) f 5(1) f 4(1) f 3(1)
        __m256d load5 = _mm256_load_pd(&this->f[offsetf + 20]);   //  f10(1) f 9(1) f 8(1) f 7(1)
        __m256d load6 = _mm256_load_pd(&this->f[offsetf + 24]);   //  f14(1) f13(1) f12(1) f11(1)
        __m256d load7 = _mm256_load_pd(&this->f[offsetf + 28]);   //  f 4(2) f 3(2) f 2(2) f 1(2)
        __m256d load8 = _mm256_load_pd(&this->f[offsetf + 32]);   //  f 8(2) f 7(2) f 6(2) f 5(2)
        __m256d load9 = _mm256_load_pd(&this->f[offsetf + 36]);   //  f12(2) f11(2) f10(2) f 9(2)
        __m256d load10 = _mm256_load_pd(&this->f[offsetf + 40]);  //  f 2(3) f 1(3) f14(2) f13(2)
        __m256d load11 = _mm256_load_pd(&this->f[offsetf + 44]);  //  f 6(3) f 5(3) f 4(3) f 3(3)
        __m256d load12 = _mm256_load_pd(&this->f[offsetf + 48]);  //  f10(3) f 9(3) f 8(3) f 7(3)
        __m256d load13 = _mm256_load_pd(&this->f[offsetf + 52]);  //  f14(3) f13(3) f12(3) f11(3)
        
        const int mm0 = 2*16 + 0*1, mm1 = 3*16 + 1*1; 
        __m256d permute0 = _mm256_permute2f128_pd(load0, load7, mm0);    //  f 2(2) f 1(2) f 2(0) f 1(0)
        __m256d permute1 = _mm256_permute2f128_pd(load0, load7, mm1);    //  f 4(2) f 3(2) f 4(0) f 3(0)
        __m256d permute2 = _mm256_permute2f128_pd(load1, load8, mm0);    //  f 6(2) f 5(2) f 6(0) f 5(0)
        __m256d permute3 = _mm256_permute2f128_pd(load1, load8, mm1);    //  f 8(2) f 7(2) f 8(0) f 7(0)
        __m256d permute4 = _mm256_permute2f128_pd(load2, load9, mm0);    //  f10(2) f 9(2) f10(0) f 9(0)
        __m256d permute5 = _mm256_permute2f128_pd(load2, load9, mm1);    //  f12(2) f11(2) f12(0) f11(0)
        __m256d permute6 = _mm256_permute2f128_pd(load3, load10, mm0);   //  f14(2) f13(2) f14(0) f13(0)
        __m256d permute7 = _mm256_permute2f128_pd(load3, load10, mm1);   //  f 2(3) f 1(3) f 2(1) f 1(1)
        __m256d permute8 = _mm256_permute2f128_pd(load4, load11, mm0);   //  f 4(3) f 3(3) f 4(1) f 3(1)
        __m256d permute9 = _mm256_permute2f128_pd(load4, load11, mm1);   //  f 6(3) f 5(3) f 6(1) f 5(1)
        __m256d permute10 = _mm256_permute2f128_pd(load5, load12, mm0);  //  f 8(3) f 7(3) f 8(1) f 7(1)
        __m256d permute11 = _mm256_permute2f128_pd(load5, load12, mm1);  //  f10(3) f 9(3) f10(1) f 9(1)
        __m256d permute12 = _mm256_permute2f128_pd(load6, load13, mm0);  //  f12(3) f11(3) f12(1) f11(1)
        __m256d permute13 = _mm256_permute2f128_pd(load6, load13, mm1);  //  f14(3) f13(3) f14(1) f13(1)

        __f[ 0] = _mm256_load_pd(&this->f0[_idx]);          //  f 0(3) f 0(2) f 0(1) f 0(0)
        __f[ 1] = _mm256_unpacklo_pd(permute0, permute7);   //  f 1(3) f 1(2) f 1(1) f 1(0)
        __f[ 2] = _mm256_unpackhi_pd(permute0, permute7);   //  f 2(3) f 2(2) f 2(1) f 2(0)
        __f[ 3] = _mm256_unpacklo_pd(permute1, permute8);   //  f 3(3) f 3(2) f 3(1) f 3(0)
        __f[ 4] = _mm256_unpackhi_pd(permute1, permute8);   //  f 4(3) f 4(2) f 4(1) f 4(0)
        __f[ 5] = _mm256_unpacklo_pd(permute2, permute9);   //  f 5(3) f 5(2) f 5(1) f 5(0)
        __f[ 6] = _mm256_unpackhi_pd(permute2, permute9);   //  f 6(3) f 6(2) f 6(1) f 6(0)
        __f[ 7] = _mm256_unpacklo_pd(permute3, permute10);  //  f 7(3) f 7(2) f 7(1) f 7(0)
        __f[ 8] = _mm256_unpackhi_pd(permute3, permute10);  //  f 8(3) f 8(2) f 8(1) f 8(0)
        __f[ 9] = _mm256_unpacklo_pd(permute4, permute11);  //  f 9(3) f 9(2) f 9(1) f 9(0)
        __f[10] = _mm256_unpackhi_pd(permute4, permute11);  //  f10(3) f10(2) f10(1) f10(0)
        __f[11] = _mm256_unpacklo_pd(permute5, permute12);  //  f11(3) f11(2) f11(1) f11(0)
        __f[12] = _mm256_unpackhi_pd(permute5, permute12);  //  f12(3) f12(2) f12(1) f12(0)
        __f[13] = _mm256_unpacklo_pd(permute6, permute13);  //  f13(3) f13(2) f13(1) f13(0)
        __f[14] = _mm256_unpackhi_pd(permute6, permute13);  //  f14(3) f14(2) f14(1) f14(0)
    }

    template<>
    template<>
    void D3Q15<double>::StoreF<__m256d>(int _idx, const __m256d *__f) {
        __m256d unpack0 = _mm256_unpacklo_pd(__f[1], __f[2]);       //  f 2(2) f 1(2) f 2(0) f 1(0)
        __m256d unpack1 = _mm256_unpackhi_pd(__f[1], __f[2]);       //  f 2(3) f 1(3) f 2(1) f 1(1)
        __m256d unpack2 = _mm256_unpacklo_pd(__f[3], __f[4]);       //  f 4(2) f 3(2) f 4(0) f 3(0)
        __m256d unpack3 = _mm256_unpackhi_pd(__f[3], __f[4]);       //  f 4(3) f 3(3) f 4(1) f 3(1)
        __m256d unpack4 = _mm256_unpacklo_pd(__f[5], __f[6]);       //  f 6(2) f 5(2) f 6(0) f 5(0)
        __m256d unpack5 = _mm256_unpackhi_pd(__f[5], __f[6]);       //  f 6(3) f 5(3) f 6(1) f 5(1)
        __m256d unpack6 = _mm256_unpacklo_pd(__f[7], __f[8]);       //  f 8(2) f 7(2) f 8(0) f 7(0)
        __m256d unpack7 = _mm256_unpackhi_pd(__f[7], __f[8]);       //  f 8(3) f 7(3) f 8(1) f 7(1)
        __m256d unpack8 = _mm256_unpacklo_pd(__f[9], __f[10]);      //  f10(2) f 9(2) f10(0) f 9(0)
        __m256d unpack9 = _mm256_unpackhi_pd(__f[9], __f[10]);      //  f10(3) f 9(3) f10(1) f 9(1)
        __m256d unpack10 = _mm256_unpacklo_pd(__f[11], __f[12]);    //  f12(2) f11(2) f12(0) f11(0)
        __m256d unpack11 = _mm256_unpackhi_pd(__f[11], __f[12]);    //  f12(3) f11(3) f12(1) f11(1)
        __m256d unpack12 = _mm256_unpacklo_pd(__f[13], __f[14]);    //  f14(2) f13(2) f14(0) f13(0)
        __m256d unpack13 = _mm256_unpackhi_pd(__f[13], __f[14]);    //  f14(3) f13(3) f14(1) f13(1)

        const int mm0 = 2*16 + 0*1, mm1 = 3*16 + 1*1, offsetf = D3Q15<double>::IndexF(_idx, 1);
        _mm256_store_pd(&this->f0[_idx], __f[0]);                                                   //  f 0(3) f 0(2) f 0(1) f 0(0)
        _mm256_store_pd(&this->f[offsetf +  0], _mm256_permute2f128_pd(unpack0, unpack2, mm0));     //  f 4(0) f 3(0) f 2(0) f 1(0)
        _mm256_store_pd(&this->f[offsetf +  4], _mm256_permute2f128_pd(unpack4, unpack6, mm0));     //  f 8(0) f 7(0) f 6(0) f 5(0)
        _mm256_store_pd(&this->f[offsetf +  8], _mm256_permute2f128_pd(unpack8, unpack10, mm0));    //  f12(0) f11(0) f10(0) f 9(0)
        _mm256_store_pd(&this->f[offsetf + 12], _mm256_permute2f128_pd(unpack12, unpack1, mm0));    //  f 2(1) f 1(1) f14(0) f13(0)
        _mm256_store_pd(&this->f[offsetf + 16], _mm256_permute2f128_pd(unpack3, unpack5, mm0));     //  f 6(1) f 5(1) f 4(1) f 3(1)
        _mm256_store_pd(&this->f[offsetf + 20], _mm256_permute2f128_pd(unpack7, unpack9, mm0));     //  f10(1) f 9(1) f 8(1) f 7(1)
        _mm256_store_pd(&this->f[offsetf + 24], _mm256_permute2f128_pd(unpack11, unpack13, mm0));   //  f14(1) f13(1) f12(1) f11(1)
        _mm256_store_pd(&this->f[offsetf + 28], _mm256_permute2f128_pd(unpack0, unpack2, mm1));     //  f 4(2) f 3(2) f 2(2) f 1(2)
        _mm256_store_pd(&this->f[offsetf + 32], _mm256_permute2f128_pd(unpack4, unpack6, mm1));     //  f 8(2) f 7(2) f 6(2) f 5(2)
        _mm256_store_pd(&this->f[offsetf + 36], _mm256_permute2f128_pd(unpack8, unpack10, mm1));    //  f12(2) f11(2) f10(2) f 9(2)
        _mm256_store_pd(&this->f[offsetf + 40], _mm256_permute2f128_pd(unpack12, unpack1, mm1));    //  f 2(3) f 1(3) f14(2) f13(2)
        _mm256_store_pd(&this->f[offsetf + 44], _mm256_permute2f128_pd(unpack3, unpack5, mm1));     //  f 6(3) f 5(3) f 4(3) f 3(3)
        _mm256_store_pd(&this->f[offsetf + 48], _mm256_permute2f128_pd(unpack7, unpack9, mm1));     //  f10(3) f 9(3) f 8(3) f 7(3)
        _mm256_store_pd(&this->f[offsetf + 52], _mm256_permute2f128_pd(unpack11, unpack13, mm1));   //  f14(3) f13(3) f12(3) f11(3)
    }
#endif
}