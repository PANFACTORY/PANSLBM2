#pragma once
#include <vector>
#include <cassert>
#include "mpi.h"

namespace PANSLBM2 {
    namespace DensityFilter {
        template<class T, template<class>class P>
        std::vector<T> GetFilteredValue(P<T>& _p, T _R, const std::vector<T> &_v) {
            assert(_R > T());
            int nR = (int)_R;

            auto IndexFX = [&](int _i, int _j) { return _j + _p.ny*_i; };
            auto IndexFY = [&](int _i, int _j) { return _i + _p.nx*_j; };
            auto IndexC = [&](int _i, int _j) { return _i + nR*_j; };

            T *send_xmin = new T[nR*_p.ny], *send_xmax = new T[nR*_p.ny];
            T *send_ymin = new T[_p.nx*nR], *send_ymax = new T[_p.nx*nR];
            T *send_xmin_ymin = new T[nR*nR], *send_xmax_ymin = new T[nR*nR], *send_xmin_ymax = new T[nR*nR], *send_xmax_ymax = new T[nR*nR];
            T *recv_xmin = new T[nR*_p.ny], *recv_xmax = new T[nR*_p.ny];
            T *recv_ymin = new T[_p.nx*nR], *recv_ymax = new T[_p.nx*nR];
            T *recv_xmin_ymin = new T[nR*nR], *recv_xmax_ymin = new T[nR*nR], *recv_xmin_ymax = new T[nR*nR], *recv_xmax_ymax = new T[nR*nR];

            MPI_Status status[52];
            MPI_Request request[52];    

            //  Copy from _v to send buffer along edge or at corner
            if (_p.PEx != 0) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int j = 0; j < _p.ny; ++j) {
                        send_xmin[IndexFX(ii, j)] = _v[_p.Index(0 + ii, j)];                                     //  Edge along xmin 
                    }
                }
            }
            if (_p.PEx != _p.mx - 1) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int j = 0; j < _p.ny; ++j) {
                        send_xmax[IndexFX(ii, j)] = _v[_p.Index((_p.nx - nR) + ii, j)];                          //  Edge along xmax
                    }
                }
            }
            if (_p.PEy != 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_ymin[IndexFY(i, jj)] = _v[_p.Index(i, 0 + jj)];                                     //  Edge along ymin 
                    }
                }
            }
            if (_p.PEy != _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_ymax[IndexFY(i, jj)] = _v[_p.Index(i, (_p.ny - nR) + jj)];                          //  Edge along ymax
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != 0) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_xmin_ymin[IndexC(ii, jj)] = _v[_p.Index(0 + ii, 0 + jj)];                          //  Corner at xmin and ymin
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_xmax_ymin[IndexC(ii, jj)] = _v[_p.Index((_p.nx - nR) + ii, 0 + jj)];               //  Corner at xmax and ymin
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_xmin_ymax[IndexC(ii, jj)] = _v[_p.Index(0 + ii, (_p.ny - nR) + jj)];               //  Corner at xmin and ymax
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1) {
                for (int ii = 0; ii < nR; ++ii) {
                    for (int jj = 0; jj < nR; ++jj) {
                        send_xmax_ymax[IndexC(ii, jj)] = _v[_p.Index((_p.nx - nR) + ii, (_p.ny - nR) + jj)];    //  Corner at xmax and ymax
                    }
                }
            }
            
            //  Communicate with other PE
            int neib = 0;
            if (_p.PEx != 0) {
                MPI_Isend(send_xmin, nR*_p.ny, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy), 0, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmin, nR*_p.ny, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy), 1, MPI_COMM_WORLD, &_p.request[neib++]);    //  Edge along xmin 
            }
            if (_p.PEx != _p.mx - 1) {
                MPI_Isend(send_xmax, nR*_p.ny, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy), 1, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmax, nR*_p.ny, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy), 0, MPI_COMM_WORLD, &_p.request[neib++]);    //  Edge along xmax
            }
            if (_p.PEy != 0) {
                MPI_Isend(send_ymin, _p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1), 2, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_ymin, _p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1), 3, MPI_COMM_WORLD, &_p.request[neib++]);    //  Edge along ymin
            }
            if (_p.PEy != _p.my - 1) {
                MPI_Isend(send_ymax, _p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1), 3, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_ymax, _p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1), 2, MPI_COMM_WORLD, &_p.request[neib++]);    //  Edge along ymax
            }
            if (_p.PEx != 0 && _p.PEy != 0) {
                MPI_Isend(send_xmin_ymin, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1), 4, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmin_ymin, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1), 7, MPI_COMM_WORLD, &_p.request[neib++]);   //  Corner at xmin and ymin
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0) {
                MPI_Isend(send_xmax_ymin, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1), 6, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmax_ymin, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1), 5, MPI_COMM_WORLD, &_p.request[neib++]);   //  Corner at xmax and ymin
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1) {
                MPI_Isend(send_xmin_ymax, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1), 5, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmin_ymax, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1), 6, MPI_COMM_WORLD, &_p.request[neib++]);   //  Corner at xmin and ymax
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1) {
                MPI_Isend(send_xmax_ymax, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1), 7, MPI_COMM_WORLD, &_p.request[neib++]);
                MPI_Irecv(recv_xmax_ymax, nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1), 4, MPI_COMM_WORLD, &_p.request[neib++]);   //  Corner at xmax and ymax
            }
            if (neib > 0) {
                MPI_Waitall(neib, request, status);
            }

            //  Filter value
            std::vector<T> fv(_p.nxy, T());
            for(int i1 = 0; i1 < _p.nx; ++i1){
                for(int j1 = 0; j1 < _p.ny; ++j1){
                    int idx = _p.Index(i1, j1);
                    T wsum = T();
                    for (int i2 = i1 - nR, i2max = i1 + nR; i2 <= i2max; ++i2) {
                        for (int j2 = j1 - nR, j2max = j1 + nR; j2 < j2max; ++j2) {                     
                            double distance = sqrt(pow(i1 - i2, 2.0) + pow(j1 - j2, 2.0));
                            if (distance <= _R) {
                                T weight = (_R - distance)/_R;
                                wsum += weight;
                                if (0 <= i2 && i2 <= _p.nx - 1) {
                                    if (0 <= j2 && j2 <= _p.ny - 1) {
                                        fv[idx] += weight*_v[_p.Index(i2, j2)];                             //  In self PE
                                    } else if (j2 < 0 && _p.PEy != 0) {
                                        fv[idx] += weight*recv_ymin[IndexFY(0 + i2, nR + j2)];               //  In ymin PE 
                                    } else if (_p.ny - 1 < j2 && _p.PEy != _p.my - 1) {
                                        fv[idx] += weight*recv_ymax[IndexFY(0 + i2, j2 - _p.ny)];            //  In ymax PE
                                    }
                                } else if (i2 < 0 && _p.PEx != 0) {
                                    if (0 <= j2 && j2 <= _p.ny - 1) {
                                        fv[idx] += weight*recv_xmin[IndexFX(nR + i2, 0 + j2)];               //  In xmin PE
                                    } else if (j2 < 0 && _p.PEy != 0) {
                                        fv[idx] += weight*recv_xmin_ymin[IndexC(nR + i2, nR + j2)];         //  In xmin ymin PE
                                    } else if (_p.ny - 1 < j2 && _p.PEy != _p.my - 1) {
                                        fv[idx] += weight*recv_xmin_ymax[IndexC(nR + i2, j2 - _p.ny)];      //  In xmin ymax PE
                                    }
                                } else if (_p.nx - 1 < i2 && this-PEx != _p.mx - 1) {
                                    if (0 <= j2 && j2 <= _p.ny - 1) {
                                        fv[idx] += weight*recv_xmax[IndexFX(i2 - _p.nx, 0 + j2)];            //  In xmax PE
                                    } else if (j2 < 0 && _p.PEy != 0) {
                                        fv[idx] += weight*recv_xmax_ymin[IndexC(i2 - _p.nx, nR + j2)];      //  In xmax ymin PE
                                    } else if (_p.ny - 1 < j2 && _p.PEy != _p.my - 1) {
                                        fv[idx] += weight*recv_xmax_ymax[IndexC(i2 - _p.nx, j2 - _p.nx)];   //  In xmax ymax PE
                                    }
                                }
                            }
                        }
                    }
                    fv[idx] /= wsum;
                }
            }

            delete[] send_xmin, send_xmax, send_ymin, send_ymax, send_xmin_ymin, send_xmax_ymin, send_xmin_ymax, send_xmax_ymax;
            delete[] recv_xmin, recv_xmax, recv_ymin, recv_ymax, recv_xmin_ymin, recv_xmax_ymin, recv_xmin_ymax, recv_xmax_ymax;

            return fv;
        }
    }
}