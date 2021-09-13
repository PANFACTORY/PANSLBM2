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

            auto IndexFX = [&](int _i, int _j, int _k) { return _j + _p.ny*_k + _p.ny*_p.nz*_i; };
            auto IndexFY = [&](int _i, int _j, int _k) { return _k + _p.nz*_i + _p.nz*_p.nx*_j; };
            auto IndexFZ = [&](int _i, int _j, int _k) { return _i + _p.nx*_j + _p.nx*_p.ny*_k; };
            auto IndexEX = [&](int _i, int _j, int _k) { return _j + nR*_k + nR*nR*_i; };
            auto IndexEY = [&](int _i, int _j, int _k) { return _k + nR*_i + nR*nR*_j; };
            auto IndexEZ = [&](int _i, int _j, int _k) { return _i + nR*_j + nR*nR*_k; };
            auto IndexCC = [&](int _i, int _j, int _k) { return _i + nR*_j + nR*nR*_k; };

            T *send_xmin = new T[_p.ny*_p.nz*nR], *send_xmax = new T[_p.ny*_p.nz*nR];
            T *send_ymin = new T[_p.nz*_p.nx*nR], *send_ymax = new T[_p.nz*_p.nx*nR];
            T *send_zmin = new T[_p.nx*_p.ny*nR], *send_zmax = new T[_p.nx*_p.ny*nR];
            T *send_ymin_zmin = new T[_p.nx*nR*nR], *send_ymax_zmin = new T[_p.nx*nR*nR], *send_ymin_zmax = new T[_p.nx*nR*nR], *send_ymax_zmax = new T[_p.nx*nR*nR];
            T *send_zmin_xmin = new T[_p.ny*nR*nR], *send_zmax_xmin = new T[_p.ny*nR*nR], *send_zmin_xmax = new T[_p.ny*nR*nR], *send_zmax_xmax = new T[_p.ny*nR*nR];
            T *send_xmin_ymin = new T[_p.nz*nR*nR], *send_xmax_ymin = new T[_p.nz*nR*nR], *send_xmin_ymax = new T[_p.nz*nR*nR], *send_xmax_ymax = new T[_p.nz*nR*nR];
            T *send_xmin_ymin_zmin = new T[nR*nR*nR], *send_xmax_ymin_zmin = new T[nR*nR*nR], *send_xmin_ymax_zmin = new T[nR*nR*nR], *send_xmax_ymax_zmin = new T[nR*nR*nR];
            T *send_xmin_ymin_zmax = new T[nR*nR*nR], *send_xmax_ymin_zmax = new T[nR*nR*nR], *send_xmin_ymax_zmax = new T[nR*nR*nR], *send_xmax_ymax_zmax = new T[nR*nR*nR]; 
            T *recv_xmin = new T[_p.ny*_p.nz*nR], *recv_xmax = new T[_p.ny*_p.nz*nR];
            T *recv_ymin = new T[_p.nz*_p.nx*nR], *recv_ymax = new T[_p.nz*_p.nx*nR];
            T *recv_zmin = new T[_p.nx*_p.ny*nR], *recv_zmax = new T[_p.nx*_p.ny*nR];
            T *recv_ymin_zmin = new T[_p.nx*nR*nR], *recv_ymax_zmin = new T[_p.nx*nR*nR], *recv_ymin_zmax = new T[_p.nx*nR*nR], *recv_ymax_zmax = new T[_p.nx*nR*nR];
            T *recv_zmin_xmin = new T[_p.ny*nR*nR], *recv_zmax_xmin = new T[_p.ny*nR*nR], *recv_zmin_xmax = new T[_p.ny*nR*nR], *recv_zmax_xmax = new T[_p.ny*nR*nR];
            T *recv_xmin_ymin = new T[_p.nz*nR*nR], *recv_xmax_ymin = new T[_p.nz*nR*nR], *recv_xmin_ymax = new T[_p.nz*nR*nR], *recv_xmax_ymax = new T[_p.nz*nR*nR];
            T *recv_xmin_ymin_zmin = new T[nR*nR*nR], *recv_xmax_ymin_zmin = new T[nR*nR*nR], *recv_xmin_ymax_zmin = new T[nR*nR*nR], *recv_xmax_ymax_zmin = new T[nR*nR*nR];
            T *recv_xmin_ymin_zmax = new T[nR*nR*nR], *recv_xmax_ymin_zmax = new T[nR*nR*nR], *recv_xmin_ymax_zmax = new T[nR*nR*nR], *recv_xmax_ymax_zmax = new T[nR*nR*nR]; 

            MPI_Status status[52];
            MPI_Request request[52];    

            //  Copy from _v to send buffer along edge or at corner
            if (_p.PEx != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmin[IndexFX(ri, j, k)] = _v[_p.Index(0 + ri, j, k)];                                     //  Face on xmin
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1) {
               for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmax[IndexFX(ri, j, k)] = _v[_p.Index((_p.nx - nR) + ri, j, k)];                          //  Face on xmax
                        }
                    }
                }
            }
            if (_p.PEy != 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_ymin[IndexFY(i, rj, k)] = _v[_p.Index(i, 0 + rj, k)];                                     //  Face on ymin 
                        }
                    }
                }
            }
            if (_p.PEy != _p.my - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_ymax[IndexFY(i, rj, k)] = _v[_p.Index(i, (_p.ny - nR) + rj, k)];                          //  Face on ymax
                        }
                    }
                }
            }
            if (_p.PEz != 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmin[IndexFZ(i, j, rk)] = _v[_p.Index(i, j, 0 + rk)];                                      //  Face on zmin
                        }
                    }
                }
            }
            if (_p.PEz != _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmax[IndexFZ(i, j, rk)] = _v[_p.Index(i, j, (_p.nz - nR) + rk)];                            //  Face on zmax
                        }
                    }
                }
            }
            if (_p.PEy != 0 && _p.PEz != 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_ymin_zmin[IndexEX(i, rj, rk)] = _v[_p.Index(i, 0 + rj, 0 + rk)];                           //  Edge along ymin zmin
                        }
                    }
                }
            }
            if (_p.PEy != _p.my - 1 && _p.PEz != 0) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_ymax_zmin[IndexEX(i, rj, rk)] = _v[_p.Index(i, (_p.ny - nR) + rj, 0 + rk)];                           //  Edge along ymin zmin
                        }
                    }
                }
            }
            if (_p.PEy != 0 && _p.PEz != _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_ymin_zmax[IndexEX(i, rj, rk)] = _v[_p.Index(i, 0 + rj, (_p.nz - nR) + rk)];                           //  Edge along ymin zmin
                        }
                    }
                }
            }
            if (_p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
                for (int i = 0; i < _p.nx; ++i) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_ymax_zmax[IndexEX(i, rj, rk)] = _v[_p.Index(i, (_p.ny - nR) + rj, (_p.nz - nR) + rk)];                           //  Edge along ymin zmin
                        }
                    }
                }
            }
            if (_p.PEz != 0 && _p.PEx != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmin_xmin[IndexEY(ri, j, rk)] = _v[_p.Index(0 + ri, j, 0 + rk)];   //  Edge along zmin xmin
                        }
                    }
                }
            }
            if (_p.PEz != _p.mz - 1 && _p.PEx != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmax_xmin[IndexEY(ri, j, rk)] = _v[_p.Index(0 + ri, j, (_p.nz - nR) + rk)];   //  Edge along zmax xmin
                        }
                    }
                }
            }
            if (_p.PEz != 0 && _p.PEx != _p.mx - 1) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmin_xmax[IndexEY(ri, j, rk)] = _v[_p.Index((_p.nx - nR) + ri, j, 0 + rk)];   //  Edge along zmin xmax
                        }
                    }
                }
            }
            if (_p.PEz != _p.mz - 1 && _p.PEx != _p.mx - 1) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int j = 0; j < _p.ny; ++j) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_zmax_xmax[IndexEY(ri, j, rk)] = _v[_p.Index((_p.nx - nR) + ri, j, (_p.nz - nR) + rk)];   //  Edge along zmax xmax
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmin_ymin[IndexEZ(ri, rj, k)] = _v[_p.Index(0 + ri, 0 + rj, k)];   //  Edge along xmin ymin
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmax_ymin[IndexEZ(ri, rj, k)] = _v[_p.Index((_p.nx - nR) + ri, 0 + rj, k)];   //  Edge along xmax ymin
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1) {
               for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmin_ymax[IndexEZ(ri, rj, k)] = _v[_p.Index(0 + ri, (_p.ny - nR) + rj, k)];   //  Edge along xmin ymax
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1) {
               for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int k = 0; k < _p.nz; ++k) {
                            send_xmax_ymax[IndexEZ(ri, rj, k)] = _v[_p.Index((_p.nx - nR) + ri, (_p.ny - nR) + rj, k)];   //  Edge along xmax ymax
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmin_ymin_zmin[IndexCC(ri, rj, rk)] = _v[_p.Index(0 + ri, 0 + rj, 0 + rk)];   //  Edge along xmin ymin zmin
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0 && _p.PEz != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmax_ymin_zmin[IndexCC(ri, rj, rk)] = _v[_p.Index((_p.nx - nR) + ri, 0 + rj, 0 + rk)];   //  Edge along xmax ymin zmin
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1 && _p.PEz != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmin_ymax_zmin[IndexCC(ri, rj, rk)] = _v[_p.Index(0 + ri, (_p.ny - nR) + rj, 0 + rk)];   //  Edge along xmin ymax zmin
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1 && _p.PEz != 0) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmax_ymax_zmin[IndexCC(ri, rj, rk)] = _v[_p.Index((_p.nx - nR) + ri, (_p.ny - nR) + rj, 0 + rk)];   //  Edge along xmax ymax zmin
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != _p.mz - 1) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmin_ymin_zmax[IndexCC(ri, rj, rk)] = _v[_p.Index(0 + ri, 0 + rj, (_p.nz - nR) + rk)];   //  Edge along xmin ymin zmax
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0 && _p.PEz != _p.mz - 1) {
               for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmax_ymin_zmax[IndexCC(ri, rj, rk)] = _v[_p.Index((_p.nx - nR) + ri, 0 + rj, (_p.nz - nR) + rk)];   //  Edge along xmax ymin zmax
                        }
                    }
                }
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
               for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmin_ymax_zmax[IndexCC(ri, rj, rk)] = _v[_p.Index(0 + ri, (_p.ny - nR) + rj, (_p.nz - nR) + rk)];   //  Edge along xmin ymax zmax
                        }
                    }
                }
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
                for (int ri = 0; ri < nR; ++ri) {
                    for (int rj = 0; rj < nR; ++rj) {
                        for (int rk = 0; rk < nR; ++rk) {
                            send_xmax_ymax_zmax[IndexCC(ri, rj, rk)] = _v[_p.Index((_p.nx - nR) + ri, (_p.ny - nR) + rj, (_p.nz - nR) + rk)];   //  Edge along xmax ymax zmax
                        }
                    }
                }
            }

            //  Communicate with other PE
            int neib = 0;
            if (_p.PEx != 0) {
                MPI_Isend(send_xmin, _p.ny*_p.nz*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz), 0, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin, _p.ny*_p.nz*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz), 1, MPI_COMM_WORLD, &request[neib++]);  //  Face on xmin
            }
            if (_p.PEx != _p.mx - 1) {
                MPI_Isend(send_xmax, _p.ny*_p.nz*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz), 1, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax, _p.ny*_p.nz*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz), 0, MPI_COMM_WORLD, &request[neib++]);  //  Face on xmax
            }
            if (_p.PEy != 0) {
                MPI_Isend(send_ymin, _p.nz*_p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz), 2, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymin, _p.nz*_p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz), 3, MPI_COMM_WORLD, &request[neib++]);  //  Face on ymin
            }
            if (_p.PEy != _p.my - 1) {
                MPI_Isend(send_ymax, _p.nz*_p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz), 3, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymax, _p.nz*_p.nx*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz), 2, MPI_COMM_WORLD, &request[neib++]);  //  Face on ymax
            }
            if (_p.PEz != 0) {
                MPI_Isend(send_zmin, _p.nx*_p.ny*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz - 1), 4, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmin, _p.nx*_p.ny*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz - 1), 5, MPI_COMM_WORLD, &request[neib++]);  //  Face on zmin
            }
            if (_p.PEz != _p.mz - 1) {
                MPI_Isend(send_zmax, _p.nx*_p.ny*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz + 1), 5, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmax, _p.nx*_p.ny*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz + 1), 4, MPI_COMM_WORLD, &request[neib++]);  //  Face on zmax
            }
            if (_p.PEy != 0 && _p.PEz != 0) {
                MPI_Isend(send_ymin_zmin, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz - 1), 6, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymin_zmin, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz - 1), 9, MPI_COMM_WORLD, &request[neib++]);    //  Edge along ymin zmin
            }
            if (_p.PEy != _p.my - 1 && _p.PEz != 0) {
                MPI_Isend(send_ymax_zmin, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz - 1), 8, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymax_zmin, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz - 1), 7, MPI_COMM_WORLD, &request[neib++]);    //  Edge along ymax zmin
            }
            if (_p.PEy != 0 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_ymin_zmax, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz + 1), 7, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymin_zmax, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz + 1), 8, MPI_COMM_WORLD, &request[neib++]);    //  Edge along ymin zmax
            }
            if (_p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_ymax_zmax, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz + 1), 9, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_ymax_zmax, _p.nx*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz + 1), 6, MPI_COMM_WORLD, &request[neib++]);    //  Edge along ymax zmax
            }
            if (_p.PEz != 0 && _p.PEx != 0) {
                MPI_Isend(send_zmin_xmin, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz - 1), 10, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmin_xmin, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz - 1), 13, MPI_COMM_WORLD, &request[neib++]);   //  Edge along zmin xmin
            }
            if (_p.PEz != _p.mz - 1 && _p.PEx != 0) {
                MPI_Isend(send_zmax_xmin, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz + 1), 12, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmax_xmin, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz + 1), 11, MPI_COMM_WORLD, &request[neib++]);   //  Edge along zmax xmin
            }
            if (_p.PEz != 0 && _p.PEx != _p.mx - 1) {
                MPI_Isend(send_zmin_xmax, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz - 1), 11, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmin_xmax, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz - 1), 12, MPI_COMM_WORLD, &request[neib++]);   //  Edge along zmin xmax
            }
            if (_p.PEz != _p.mz - 1 && _p.PEx != _p.mx - 1) {
                MPI_Isend(send_zmax_xmax, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz + 1), 13, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_zmax_xmax, _p.ny*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz + 1), 10, MPI_COMM_WORLD, &request[neib++]);   //  Edge along zmax xmax
            }
            if (_p.PEx != 0 && _p.PEy != 0) {
                MPI_Isend(send_xmin_ymin, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz), 14, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin_ymin, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz), 17, MPI_COMM_WORLD, &request[neib++]);   //  Edge along xmin ymin
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0) {
                MPI_Isend(send_xmax_ymin, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz), 16, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax_ymin, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz), 15, MPI_COMM_WORLD, &request[neib++]);   //  Edge along xmax ymin
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1) {
                MPI_Isend(send_xmin_ymax, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz), 15, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin_ymax, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz), 16, MPI_COMM_WORLD, &request[neib++]);   //  Edge along xmin ymax
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1) {
                MPI_Isend(send_xmax_ymax, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz), 17, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax_ymax, _p.nz*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz), 14, MPI_COMM_WORLD, &request[neib++]);   //  Edge along xmax ymax
            }
            if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != 0) {
                MPI_Isend(send_xmin_ymin_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz - 1), 18, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin_ymin_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz - 1), 18, MPI_COMM_WORLD, &request[neib++]); 
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0 && _p.PEz != 0) {
                MPI_Isend(send_xmax_ymin_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz - 1), 19, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax_ymin_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz - 1), 19, MPI_COMM_WORLD, &request[neib++]); 
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1 && _p.PEz != 0) {
                MPI_Isend(send_xmin_ymax_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz - 1), 20, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin_ymax_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz - 1), 20, MPI_COMM_WORLD, &request[neib++]);
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1 && _p.PEz != 0) {
                MPI_Isend(send_xmax_ymax_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz - 1), 21, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax_ymax_zmin, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz - 1), 21, MPI_COMM_WORLD, &request[neib++]);  
            }
            if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_xmin_ymin_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz + 1), 22, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmin_ymin_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz + 1), 22, MPI_COMM_WORLD, &request[neib++]); 
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != 0 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_xmax_ymin_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz + 1), 23, MPI_COMM_WORLD, &request[neib++]);
                MPI_Irecv(recv_xmax_ymin_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy - 1, _p.PEz + 1), 23, MPI_COMM_WORLD, &request[neib++]);
            }
            if (_p.PEx != 0 && _p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_xmin_ymax_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz + 1), 24, MPI_COMM_WORLD, &request[neib++]); 
                MPI_Irecv(recv_xmin_ymax_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy + 1, _p.PEz + 1), 24, MPI_COMM_WORLD, &request[neib++]); 
            }
            if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
                MPI_Isend(send_xmax_ymax_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz + 1), 25, MPI_COMM_WORLD, &request[neib++]); 
                MPI_Irecv(recv_xmax_ymax_zmax, nR*nR*nR, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz + 1), 25, MPI_COMM_WORLD, &request[neib++]);
            }
            if (neib > 0) {
                MPI_Waitall(neib, request, status);
            }

            //  Filter value
            std::vector<T> fv(_p.nxyz, T());
            for(int i1 = 0; i1 < _p.nx; ++i1){
                for(int j1 = 0; j1 < _p.ny; ++j1){
                    for(int k1 = 0; k1 < _p.nz; ++k1){
                        int idx = _p.Index(i1, j1, k1);
                        T wsum = T();
                        for (int i2 = i1 - nR, i2max = i1 + nR; i2 <= i2max; ++i2) {
                            for (int j2 = j1 - nR, j2max = j1 + nR; j2 < j2max; ++j2) {
                                for (int k2 = k1 - nR, k2max = k1 + nR; k2 < k2max; ++k2) {                     
                                    double distance = sqrt(pow(i1 - i2, 2.0) + pow(j1 - j2, 2.0) + pow(k1 - k2, 2.0));
                                    if (distance <= _R) {
                                        T weight = (_R - distance)/_R;
                                        wsum += weight;
                                        if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                                                //  In self PE
                                            fv[idx] += weight*_v[_p.Index(i2, j2, k2)];                             
                                        } else if ((i2 < 0 && _p.PEx != 0) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                                         //  In xmin PE
                                            fv[idx] += weight*recv_xmin[IndexFX(nR + i2, j2, k2)];  
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (0 <= j2 && j2 < _p.ny) && (0 <= k2 && k2 < _p.nz)) {                            //  In xmax PE
                                            fv[idx] += weight*recv_xmax[IndexFX(i2 - _p.nx, j2, k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (j2 < 0 && _p.PEy != 0) && (0 <= k2 && k2 < _p.nz)) {                                         //  In ymin PE
                                            fv[idx] += weight*recv_ymin[IndexFY(i2, nR + j2, k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (0 <= k2 && k2 < _p.nz)) {                            //  In ymax PE
                                            fv[idx] += weight*recv_ymax[IndexFY(i2, j2 - _p.ny, k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (k2 < 0 && _p.PEz != 0)) {                                         //  In zmin PE
                                            fv[idx] += weight*recv_zmin[IndexFZ(i2, j2, nR + k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (0 <= j2 && j2 < _p.ny) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {                            //  In zmax PE
                                            fv[idx] += weight*recv_zmax[IndexFZ(i2, j2, k2 - _p.nz)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (j2 < 0 && _p.PEy != 0) && (k2 < 0 && _p.PEz != 0)) {                                         //  In ymin zmin PE
                                            fv[idx] += weight*recv_ymin_zmin[IndexEX(i2, nR + j2, nR + k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (k2 < 0 && _p.PEz != 0)) {                            //  In ymax zmin PE
                                            fv[idx] += weight*recv_ymax_zmin[IndexEX(i2, j2 - _p.ny, nR + k2)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (j2 < 0 && _p.PEy != 0) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {                            //  In ymin zmax PE
                                            fv[idx] += weight*recv_ymin_zmax[IndexEX(i2, nR + j2, k2 - _p.nz)];
                                        } else if ((0 <= i2 && i2 < _p.nx) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {               //  In ymax zmax PE
                                            fv[idx] += weight*recv_ymax_zmax[IndexEX(i2, j2 - _p.ny, k2 - _p.nz)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (0 <= j2 && j2 < _p.ny) && (k2 < 0 && _p.PEz != 0)) {                                         //  In zmin xmin PE
                                            fv[idx] += weight*recv_zmin_xmin[IndexEY(nR + i2, j2, nR + k2)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (0 <= j2 && j2 < _p.ny) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {                            //  In zmax xmin PE
                                            fv[idx] += weight*recv_zmax_xmin[IndexEY(nR + i2, j2, k2 - _p.nz)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (0 <= j2 && j2 < _p.ny) && (k2 < 0 && _p.PEz != 0)) {                            //  In zmin xmax PE
                                            fv[idx] += weight*recv_zmin_xmax[IndexEY(i2 - _p.nx, j2, nR + k2)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (0 <= j2 && j2 < _p.ny) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {               //  In zmax xmax PE
                                            fv[idx] += weight*recv_zmax_xmax[IndexEY(i2 - _p.nx, j2, k2 - _p.nz)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (j2 < 0 && _p.PEy != 0) && (0 <= k2 && k2 < _p.nz)) {                                         //  In xmin ymin PE
                                            fv[idx] += weight*recv_xmin_ymin[IndexEZ(nR + i2, nR + j2, k2)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (j2 < 0 && _p.PEy != 0) && (0 <= k2 && k2 < _p.nz)) {                            //  In xmax ymin PE
                                            fv[idx] += weight*recv_xmax_ymin[IndexEZ(i2 - _p.nx, nR + j2, k2)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (0 <= k2 && k2 < _p.nz)) {                            //  In xmin ymax PE
                                            fv[idx] += weight*recv_xmin_ymax[IndexEZ(nR + i2, j2 - _p.ny, k2)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (0 <= k2 && k2 < _p.nz)) {               //  In xmax ymax PE
                                            fv[idx] += weight*recv_xmax_ymax[IndexEZ(i2 - _p.nx, j2 - _p.ny, k2)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (j2 < 0 && _p.PEy != 0) && (k2 < 0 && _p.PEz != 0)) {                                         //  In xmin ymin zmin PE
                                            fv[idx] += weight*recv_xmin_ymin_zmin[IndexCC(nR + i2, nR + j2, nR + k2)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (j2 < 0 && _p.PEy != 0) && (k2 < 0 && _p.PEz != 0)) {                            //  In xmax ymin zmin PE
                                            fv[idx] += weight*recv_xmax_ymin_zmin[IndexCC(i2 - _p.nx, nR + j2, nR + k2)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (k2 < 0 && _p.PEz != 0)) {                            //  In xmin ymax zmin PE
                                            fv[idx] += weight*recv_xmin_ymax_zmin[IndexCC(nR + i2, j2 - _p.ny, nR + k2)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (k2 < 0 && _p.PEz != 0)) {               //  In xmax ymax zmin PE
                                            fv[idx] += weight*recv_xmax_ymax_zmin[IndexCC(i2 - _p.nx, j2 - _p.ny, nR + k2)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (j2 < 0 && _p.PEy != 0) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {                            //  In xmin ymin zmax PE
                                            fv[idx] += weight*recv_xmin_ymin_zmax[IndexCC(nR + i2, nR + j2, k2 - _p.nz)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (j2 < 0 && _p.PEy != 0) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {               //  In xmax ymin zmax PE
                                            fv[idx] += weight*recv_xmax_ymin_zmax[IndexCC(i2 - _p.nx, nR + j2, k2 - _p.nz)];
                                        } else if ((i2 < 0 && _p.PEx != 0) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {               //  In xmin ymax zmax PE
                                            fv[idx] += weight*recv_xmin_ymax_zmax[IndexCC(nR + i2, j2 - _p.ny, k2 - _p.nz)];
                                        } else if ((_p.nx <= i2 && _p.PEx != _p.mx - 1) && (_p.ny <= j2 && _p.PEy != _p.my - 1) && (_p.nz <= k2 && _p.PEz != _p.mz - 1)) {  //  In xmax ymax zmax PE
                                            fv[idx] += weight*recv_xmax_ymax_zmax[IndexCC(i2 - _p.nx, j2 - _p.ny, k2 - _p.nz)];
                                        }
                                    }
                                }
                            }
                        }
                        fv[idx] /= wsum;
                    }
                }
            }

            delete[] send_xmin, send_xmax, send_ymin, send_ymax, send_zmin, send_zmax;
            delete[] send_ymin_zmin, send_ymax_zmin, send_ymin_zmax, send_ymax_zmax; 
            delete[] send_zmin_xmin, send_zmax_xmin, send_zmin_xmax, send_zmax_xmax;
            delete[] send_xmin_ymin, send_xmax_ymin, send_xmin_ymax, send_xmax_ymax;
            delete[] send_xmin_ymin_zmin, send_xmax_ymin_zmin, send_xmin_ymax_zmin, send_xmax_ymax_zmin, send_xmin_ymin_zmax, send_xmax_ymin_zmax, send_xmin_ymax_zmax, send_xmax_ymax_zmax; 
            delete[] recv_xmin, recv_xmax, recv_ymin, recv_ymax, recv_zmin, recv_zmax;
            delete[] recv_ymin_zmin, recv_ymax_zmin, recv_ymin_zmax, recv_ymax_zmax;
            delete[] recv_zmin_xmin, recv_zmax_xmin, recv_zmin_xmax, recv_zmax_xmax;
            delete[] recv_xmin_ymin, recv_xmax_ymin, recv_xmin_ymax, recv_xmax_ymax;
            delete[] recv_xmin_ymin_zmin, recv_xmax_ymin_zmin, recv_xmin_ymax_zmin, recv_xmax_ymax_zmin, recv_xmin_ymin_zmax, recv_xmax_ymin_zmax, recv_xmin_ymax_zmax, recv_xmax_ymax_zmax; 

            return fv;
        }
    }
}