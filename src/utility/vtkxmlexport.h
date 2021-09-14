#pragma once
#include <string>
#include <fstream>
#include <cassert>
#ifdef _USE_MPI_DEFINES
    #include "mpi.h"
#endif

namespace PANSLBM2 {
    class VTKXMLExport {
public:
        VTKXMLExport() = delete;
        template<class T, template<class>class P>
        VTKXMLExport(P<T>& _p, std::string _fname) {
            this->isaddpointdata = false;

            //  Export .pvts file (Only host process)
            if (_p.PEid == 0) {
                int pos = std::max((int)_fname.rfind('\\'), (int)_fname.rfind('/'));
                std::string fname = _fname.substr(pos != std::string::npos ? pos + 1 : 0);
            
                this->fout0.open(_fname + ".pvts");
                this->fout0 << "<VTKFile type=\"PStructuredGrid\">" << std::endl;
                this->fout0 << "\t<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 " << _p.lx - 1 << " 0 " << _p.ly - 1 << " 0 " << _p.lz - 1 << "\">" << std::endl;
                this->fout0 << "\t\t<PPoints>" << std::endl;
                this->fout0 << "\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
                this->fout0 << "\t\t</PPoints>" << std::endl;
                for (int k = 0; k < _p.mz; ++k) {
                    for (int j = 0; j < _p.my; ++j) {
                        for (int i = 0; i < _p.mx; ++i) {
                            int peid = i + _p.mx*j + _p.mx*_p.my*k;
                            int nnx = (_p.lx + i)/_p.mx;
                            int nny = (_p.ly + j)/_p.my;
                            int nnz = (_p.lz + k)/_p.mz;
                            int osx = _p.mx - i > _p.lx%_p.mx ? i*nnx : _p.lx - (_p.mx - i)*nnx;
                            int osy = _p.my - j > _p.ly%_p.my ? j*nny : _p.ly - (_p.my - j)*nny;
                            int osz = _p.mz - k > _p.lz%_p.mz ? k*nnz : _p.lz - (_p.mz - k)*nnz;
                            nnx += i != _p.mx - 1 ? 1 : 0;
                            nny += j != _p.my - 1 ? 1 : 0;
                            nnz += k != _p.mz - 1 ? 1 : 0;
                            this->fout0 << "\t\t<Piece Extent=\"" 
                                << 0 + osx << " " << (nnx - 1) + osx << " " 
                                << 0 + osy << " " << (nny - 1) + osy << " " 
                                << 0 + osz << " " << (nnz - 1) + osz 
                                << "\" Source=\"" << fname + "_" + std::to_string(peid) + ".vts" << "\"/>" << std::endl;
                        }
                    }
                }
                this->addpos0 = this->fout0.tellp();
                this->fout0 << "\t</PStructuredGrid>" << std::endl;
                this->fout0 << "</VTKFile>";
            }

            //  Export .vts file
            this->fout.open(_fname + "_" + std::to_string(_p.PEid) + ".vts");
            this->fout << "<VTKFile type=\"StructuredGrid\">" << std::endl;
            this->fout << "\t<StructuredGrid WholeExtent=\"0 "<< _p.lx - 1 << " 0 " << _p.ly - 1 << " 0 " << _p.lz - 1 << "\">" << std::endl;
            int nnx = _p.nx + (_p.PEx != _p.mx - 1 ? 1 : 0);
            int nny = _p.ny + (_p.PEy != _p.my - 1 ? 1 : 0);
            int nnz = _p.nz + (_p.PEz != _p.mz - 1 ? 1 : 0);
            this->fout << "\t\t<Piece Extent=\"" 
                << 0 + _p.offsetx << " " << (nnx - 1) + _p.offsetx << " " 
                << 0 + _p.offsety << " " << (nny - 1) + _p.offsety << " " 
                << 0 + _p.offsetz << " " << (nnz - 1) + _p.offsetz << "\">" << std::endl;
            this->fout << "\t\t\t<Points>" << std::endl; 
            this->fout << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
            for (int k = 0; k < nnz; ++k) {
                for (int j = 0; j < nny; ++j) {
                    for (int i = 0; i < nnx; ++i) {
                        this->fout << "\t\t\t\t\t" << i + _p.offsetx << " " << j + _p.offsety << " " << k + _p.offsetz << std::endl;
                    }
                }
            }
            this->fout << "\t\t\t\t</DataArray>" << std::endl;
            this->fout << "\t\t\t</Points>" << std::endl;
            this->addpos = this->fout.tellp();
            this->fout << "\t\t</Piece>" << std::endl;
            this->fout << "\t</StructuredGrid>" << std::endl;
            this->fout << "</VTKFile>";
        }
        VTKXMLExport(const VTKXMLExport&) = delete;
        ~VTKXMLExport() {}

        template<class T, template<class>class P, class... Fs>
        void AddPointData(P<T>& _p, std::string _label, Fs ..._values);
        
private:
        std::ofstream fout0, fout;
        std::ofstream::pos_type addpos0, addpos;
        bool isaddpointdata;

        template<class T>
        void Values(T *_v, int _i, int _j, int _k) {}
        template<class T, class F, class... Fs>
        void Values(T *_v, int _i, int _j, int _k, F _value, Fs... _values) {
            *_v = _value(_i, _j, _k);
            this->Values(++_v, _i, _j, _k, _values...);
        }
    };

    template<class T, template<class>class P, class... Fs>
    void VTKXMLExport::AddPointData(P<T>& _p, std::string _label, Fs ..._values) {
        int cmp = sizeof...(Fs);
        assert(cmp == 1 || cmp == 3 || cmp == 9);
        T values[9] = { T() };

#ifdef _USE_MPI_DEFINES
        auto IndexFX = [&](int _j, int _k, int _c) { return _c + cmp*(_j + _p.ny*_k); };
        auto IndexFY = [&](int _k, int _i, int _c) { return _c + cmp*(_k + _p.nz*_i); };
        auto IndexFZ = [&](int _i, int _j, int _c) { return _c + cmp*(_i + _p.nx*_j); };
        auto IndexEX = [&](int _i, int _c) { return _c + cmp*_i; };
        auto IndexEY = [&](int _j, int _c) { return _c + cmp*_j; };
        auto IndexEZ = [&](int _k, int _c) { return _c + cmp*_k; };

        T *send_xmin = new T[_p.ny*_p.nz*cmp], *send_ymin = new T[_p.nz*_p.nx*cmp], *send_zmin = new T[_p.nx*_p.ny*cmp];
        T *send_ymin_zmin = new T[_p.nx*cmp], *send_zmin_xmin = new T[_p.ny*cmp], *send_xmin_ymin = new T[_p.nx*cmp];
        T *send_xmin_ymin_zmin = new T[cmp];
        T *recv_xmax = new T[_p.ny*_p.nz*cmp], *recv_ymax = new T[_p.nz*_p.nx*cmp], *recv_zmax = new T[_p.nx*_p.ny*cmp];
        T *recv_ymax_zmax = new T[_p.nx*cmp], *recv_zmax_xmax = new T[_p.ny*cmp], *recv_xmax_ymax = new T[_p.nx*cmp];
        T *recv_xmax_ymax_zmax = new T[cmp];

        MPI_Status status[14];
        MPI_Request request[14];

        //  Copy from _value to send buffer along edge or at corner
        if (_p.PEx != 0) {
            for (int j = 0; j < _p.ny; ++j) {
                for (int k = 0; k < _p.nz; ++k) {
                    this->Values(&send_xmin[IndexFX(j, k, 0)], 0, j, k, _values...);
                }
            }
        }
        if (_p.PEy != 0) {
            for (int k = 0; k < _p.nz; ++k) {
                for (int i = 0; i < _p.nx; ++i) {
                    this->Values(&send_ymin[IndexFY(k, i, 0)], i, 0, k, _values...);
                }
            }
        }
        if (_p.PEz != 0) {
            for (int i = 0; i < _p.nx; ++i) {
                for (int j = 0; j < _p.ny; ++j) {
                    this->Values(&send_zmin[IndexFZ(i, j, 0)], i, j, 0, _values...);
                }
            }
        }
        if (_p.PEy != 0 && _p.PEz != 0) {
            for (int i = 0; i < _p.nx; ++i) {
                this->Values(&send_ymin_zmin[IndexEX(i, c)], i, 0, 0, _values...);
            }
        }
        if (_p.PEz != 0 && _p.PEx != 0) {
            for (int j = 0; j < _p.ny; ++j) {
                this->Values(&send_zmin_xmin[IndexEY(j, c)], 0, j, 0, _values...);
            }
        }
        if (_p.PEx != 0 && _p.PEy != 0) {
            for (int k = 0; k < _p.nz; ++k) {
                this->Values(&send_xmin_ymin[IndexEZ(k, c)], 0, 0, k, _values...);
            }
        }
        if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != 0) {
            this->Values(send_xmin_ymin_zmin, 0, 0, 0, _values...);
        }

        //  Communicate with other PE
        int neib = 0;
        if (_p.PEx != 0) {
            MPI_Isend(send_xmin, _p.ny*_p.nz*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz), 0, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEx != _p.mx - 1) {
            MPI_Irecv(recv_xmax, _p.ny*_p.nz*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz), 0, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEy != 0) {
            MPI_Isend(send_ymin, _p.nz*_p.nx*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz), 1, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEy != _p.my - 1) {
            MPI_Irecv(recv_ymax, _p.nz*_p.nx*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz), 1, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEz != 0) {
            MPI_Isend(send_zmin, _p.nx*_p.ny*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz - 1), 2, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEz != _p.mz - 1) {
            MPI_Irecv(recv_zmax, _p.nx*_p.ny*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy, _p.PEz + 1), 2, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEy != 0 && _p.PEz != 0) {
            MPI_Isend(send_ymin_zmin, _p.nx*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy - 1, _p.PEz - 1), 3, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
            MPI_Irecv(recv_ymax_zmax, _p.nx*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx, _p.PEy + 1, _p.PEz + 1), 3, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEz != 0 && _p.PEx != 0) {
            MPI_Isend(send_zmin_xmin, _p.ny*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy, _p.PEz - 1), 4, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEz != _p.mz - 1 && _p.PEx != _p.mx - 1) {
            MPI_Irecv(recv_zmax_xmax, _p.ny*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy, _p.PEz + 1), 4, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEx != 0 && _p.PEy != 0) {
            MPI_Isend(send_xmin_ymin, _p.nz*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz), 5, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1) {
            MPI_Irecv(recv_xmax_ymax, _p.nz*cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz), 5, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEx != 0 && _p.PEy != 0 && _p.PEz != 0) {
            MPI_Isend(send_xmin_ymin_zmin, cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx - 1, _p.PEy - 1, _p.PEz - 1), 6, MPI_COMM_WORLD, &request[neib++]);
        }
        if (_p.PEx != _p.mx - 1 && _p.PEy != _p.my - 1 && _p.PEz != _p.mz - 1) {
            MPI_Irecv(recv_xmax_ymax_zmax, cmp, MPI_DOUBLE, _p.IndexPE(_p.PEx + 1, _p.PEy + 1, _p.PEz + 1), 6, MPI_COMM_WORLD, &request[neib++]);
        }
        if (neib > 0) {
            MPI_Waitall(neib, request, status);
        }
#endif
        //  Add header
        if (!this->isaddpointdata) {
            this->isaddpointdata = true;

            //  Export header of .pvts file
            if (_p.PEid == 0) {
                this->fout0.seekp(this->addpos0);
                this->fout0 << "\t\t<PPointData>" << std::endl;
                this->addpos0 = this->fout0.tellp();
            }

            //  Export header of .vts file
            this->fout.seekp(this->addpos);
            this->fout << "\t\t\t<PointData>" << std::endl;
            this->addpos = this->fout.tellp();
        }

        //  Export .pvts file
        if (_p.PEid == 0) {
            this->fout0.seekp(this->addpos0);
            this->fout0 << "\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"" << cmp << "\" format=\"ascii\"/>" << std::endl;
            this->addpos0 = this->fout0.tellp();
        }

        //  Export .vts file
        int nnx = _p.nx + (_p.PEx != _p.mx - 1 ? 1 : 0);
        int nny = _p.ny + (_p.PEy != _p.my - 1 ? 1 : 0);
        int nnz = _p.nz + (_p.PEz != _p.mz - 1 ? 1 : 0);
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"" << cmp << "\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < nnz; ++k) {
            for (int j = 0; j < nny; ++j) {
                for (int i = 0; i < nnx; ++i) {
                    this->fout << "\t\t\t\t\t";
#ifdef _USE_MPI_DEFINES                     
                    if ((_p.PEx != _p.mx - 1 && i == nnx - 1) && (_p.PEy != _p.my - 1 && j == nny - 1) && (_p.PEz != _p.mz - 1 && k == nnz - 1)) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_xmax_ymax_zmax[c] << " ";
                        }
                    } else if ((_p.PEy != _p.my - 1 && j == nny - 1) && (_p.PEz != _p.mz - 1 && k == nnz - 1)) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_ymax_zmax[IndexEX(i, c)] << " ";
                        }
                    } else if ((_p.PEz != _p.mz - 1 && k == nnz - 1) && (_p.PEx != _p.mx - 1 && i == nnx - 1)) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_zmax_xmax[IndexEY(j, c)] << " ";
                        }
                    } else if ((_p.PEx != _p.mx - 1 && i == nnx - 1) && (_p.PEy != _p.my - 1 && j == nny - 1)) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_xmax_ymax[IndexEZ(k, c)] << " ";
                        }
                    } else if (_p.PEx != _p.mx - 1 && i == nnx - 1) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_xmax[IndexFX(j, k, c)] << " ";
                        }
                    } else if (_p.PEy != _p.my - 1 && j == nny - 1) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_ymax[IndexFY(k, i, c)] << " ";
                        }
                    } else if (_p.PEz != _p.mz - 1 && k == nnz - 1) {
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << recv_zmax[IndexFZ(i, j, c)] << " ";
                        }
                    } else {
#endif
                        this->Values(values, i, j, k, _values...);
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << values[c] << " ";
                        }
#ifdef _USE_MPI_DEFINES
                    }
#endif
                    this->fout << std::endl;
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        //  Export Footer of .pvts file
        if (_p.PEid == 0) {
            this->fout0 << "\t\t</PPointData>" << std::endl;
            this->fout0 << "\t</PStructuredGrid>" << std::endl;
            this->fout0 << "</VTKFile>";
        }

        //  Export Footer of .vts file
        this->fout << "\t\t\t</PointData>" << std::endl;
        this->fout << "\t\t</Piece>" << std::endl;
        this->fout << "\t</StructuredGrid>" << std::endl;
        this->fout << "</VTKFile>";
#ifdef _USE_MPI_DEFINES
        delete[] send_xmin, send_ymin, send_zmin, send_ymin_zmin, send_zmin_xmin, send_xmin_ymin, send_xmin_ymin_zmin;
        delete[] recv_xmax, recv_ymax, recv_zmax, recv_ymax_zmax, recv_zmax_xmax, recv_xmax_ymax, recv_xmax_ymax_zmax;
#endif
    }
}