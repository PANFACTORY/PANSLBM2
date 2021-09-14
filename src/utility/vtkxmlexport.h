#pragma once
#include <string>
#include <fstream>
#include <cassert>

namespace PANSLBM2 {
    class VTKXMLExport {
public:
        VTKXMLExport() = delete;
        VTKXMLExport(std::string _fname, int _PEid, int _lx, int _ly, int _lz, int _mx, int _my, int _mz) :
            lx(_lx), ly(_ly), lz(_lz), mx(_mx), my(_my), mz(_mz), PEid(_PEid),
            PEx(this->PEid%this->mx), PEy((this->PEid/this->mx)%this->my), PEz(this->PEid/(this->mx*this->my)),
            nx((this->lx + this->PEx)/this->mx), ny((this->ly + this->PEy)/this->my), nz((this->lz + this->PEz)/this->mz),
            offsetx(this->mx - this->PEx > this->lx%this->mx ? this->PEx*this->nx : this->lx - (this->mx - this->PEx)*this->nx),
            offsety(this->my - this->PEy > this->ly%this->my ? this->PEy*this->ny : this->ly - (this->my - this->PEy)*this->ny),
            offsetz(this->mz - this->PEz > this->lz%this->mz ? this->PEz*this->nz : this->lz - (this->mz - this->PEz)*this->nz)
        {
            this->isaddpointdata = false;

            //  Export .pvts file (Only host process)
            if (this->PEid == 0) {
                int pos = std::max((int)_fname.rfind('\\'), (int)_fname.rfind('/'));
                std::string fname = _fname.substr(pos != std::string::npos ? pos + 1 : 0);
            
                this->fout0.open(_fname + ".pvts");
                this->fout0 << "<VTKFile type=\"PStructuredGrid\">" << std::endl;
                this->fout0 << "\t<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 " << this->lx - 1 << " 0 " << this->ly - 1 << " 0 " << this->lz - 1 << "\">" << std::endl;
                this->fout0 << "\t\t<PPoints>" << std::endl;
                this->fout0 << "\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
                this->fout0 << "\t\t</PPoints>" << std::endl;
                for (int k = 0; k < this->mz; ++k) {
                    for (int j = 0; j < this->my; ++j) {
                        for (int i = 0; i < this->mx; ++i) {
                            int peid = i + this->mx*j + this->mx*this->my*k;
                            int nnx = (this->lx + i)/this->mx;
                            int nny = (this->ly + j)/this->my;
                            int nnz = (this->lz + k)/this->mz;
                            int osx = this->mx - i > this->lx%this->mx ? i*nnx : this->lx - (this->mx - i)*nnx;
                            int osy = this->my - j > this->ly%this->my ? j*nny : this->ly - (this->my - j)*nny;
                            int osz = this->mz - k > this->lz%this->mz ? k*nnz : this->lz - (this->mz - k)*nnz;
                            nnx += i != this->mx - 1 ? 1 : 0;
                            nny += j != this->my - 1 ? 1 : 0;
                            nnz += k != this->mz - 1 ? 1 : 0;
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
            this->fout.open(_fname + "_" + std::to_string(this->PEid) + ".vts");
            this->fout << "<VTKFile type=\"StructuredGrid\">" << std::endl;
            this->fout << "\t<StructuredGrid WholeExtent=\"0 "<< this->lx - 1 << " 0 " << this->ly - 1 << " 0 " << this->lz - 1 << "\">" << std::endl;
            int nnx = this->nx + (this->PEx != this->mx - 1 ? 1 : 0);
            int nny = this->ny + (this->PEy != this->my - 1 ? 1 : 0);
            int nnz = this->nz + (this->PEz != this->mz - 1 ? 1 : 0);
            this->fout << "\t\t<Piece Extent=\"" 
                << 0 + this->offsetx << " " << (nnx - 1) + this->offsetx << " " 
                << 0 + this->offsety << " " << (nny - 1) + this->offsety << " " 
                << 0 + this->offsetz << " " << (nnz - 1) + this->offsetz << "\">" << std::endl;
            this->fout << "\t\t\t<Points>" << std::endl; 
            this->fout << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
            for (int k = 0; k < nnz; ++k) {
                for (int j = 0; j < nny; ++j) {
                    for (int i = 0; i < nnx; ++i) {
                        this->fout << "\t\t\t\t\t" << i + this->offsetx << " " << j + this->offsety << " " << k + this->offsetz << std::endl;
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

        template<class... Fs>
        void AddPointData(std::string _label, Fs ..._values);
        
private:
        std::ofstream fout0, fout;
        const int lx, ly, lz, mx, my, mz, PEid, PEx, PEy, PEz, nx, ny, nz, offsetx, offsety, offsetz;
        std::ofstream::pos_type addpos0, addpos;
        bool isaddpointdata;

        void AddHeader() {
            if (!this->isaddpointdata) {
                this->isaddpointdata = true;

                //  Export header of .pvts file
                if (this->PEid == 0) {
                    this->fout0.seekp(this->addpos0);
                    this->fout0 << "\t\t<PPointData>" << std::endl;
                    this->addpos0 = this->fout0.tellp();
                }

                //  Export header of .vts file
                this->fout.seekp(this->addpos);
                this->fout << "\t\t\t<PointData>" << std::endl;
                this->addpos = this->fout.tellp();
            }
        }
        void AddFooter() {
            //  Export Footer of .pvts file
            if (this->PEid == 0) {
                this->fout0 << "\t\t</PPointData>" << std::endl;
                this->fout0 << "\t</PStructuredGrid>" << std::endl;
                this->fout0 << "</VTKFile>";
            }

            //  Export Footer of .vts file
            this->fout << "\t\t\t</PointData>" << std::endl;
            this->fout << "\t\t</Piece>" << std::endl;
            this->fout << "\t</StructuredGrid>" << std::endl;
            this->fout << "</VTKFile>";
        }
    
        void Values(int _i, int _j, int _k) {}
        template <class F, class... Fs>
        void Values(int _i, int _j, int _k, F _value, Fs... _values) {
            this->fout << _value(_i, _j, _k) << " ";
            this->Values(_i, _j, _k, _values...);
        }

        //void Synchronize();
    };

    /*void VTKXMLExport::Synchronize() {
        int cmp;
        T *send_xmin = new T[this->ny*this->nz*cmp];
        T *send_ymin = new T[this->nz*this->nx*cmp];
        T *send_zmin = new T[this->nx*this->ny*cmp];
        T *send_ymin_zmin = new T[this->nx*cmp];
        T *send_zmin_xmin = new T[this->ny*cmp];
        T *send_xmin_ymin = new T[this->nx*cmp];
        T *send_xmin_ymin_zmin = new T[cmp];
        T *recv_xmax = new T[this->ny*this->nz*cmp];
        T *recv_ymax = new T[this->nz*this->nx*cmp];
        T *recv_zmax = new T[this->nx*this->ny*cmp];
        T *recv_ymax_zmax = new T[this->nx*cmp];
        T *recv_zmax_xmax = new T[this->ny*cmp];
        T *recv_xmax_ymax = new T[this->nx*cmp];
        T *recv_xmax_ymax_zmax = new T[cmp];

        MPI_Status status[14];
        MPI_Request request[14];

        //  Copy from _value to send buffer along edge or at corner
        if (this->PEx != 0) {
            for (int j = 0; j < this->ny; ++j) {
                for (int k = 0; k < this->nz; ++k) {
                    for (int c = 0; c < cmp; ++c) {
                        send_xmin[] = F(0, j, k);
                    }
                }
            }
        }

        //  Communicate with other PE
        int neib = 0;
        if (this->PEx != 0) {
            MPI_Isend(send_xmin, this->ny*this->nz*cmp, MPI_DOUBLE, this->IndexPE(this->PEx - 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &request[neib++]);
        }
        if (this->PEx != this->mx - 1) {
            MPI_Irecv(recv_xmax, this->ny*this->nz*cmp, MPI_DOUBLE, this->IndexPE(this->PEx + 1, this->PEy, this->PEz), 0, MPI_COMM_WORLD, &request[neib++]);
        }

        //  

    }*/

    

    template<class... Fs>
    void VTKXMLExport::AddPointData(std::string _label, Fs ..._values) {
        int cmp = sizeof...(Fs);
        assert(cmp == 1 || cmp == 3 || cmp == 9);
        
        this->AddHeader();

        //  Export .pvts file
        if (this->PEid == 0) {
            this->fout0.seekp(this->addpos0);
            this->fout0 << "\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"" << cmp << "\" format=\"ascii\"/>" << std::endl;
            this->addpos0 = this->fout0.tellp();
        }

        //  Export .vts file
        int nnx = this->nx + (this->PEx != this->mx - 1 ? 1 : 0);
        int nny = this->ny + (this->PEy != this->my - 1 ? 1 : 0);
        int nnz = this->nz + (this->PEz != this->mz - 1 ? 1 : 0);
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"" << cmp << "\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < nnz; ++k) {
            for (int j = 0; j < nny; ++j) {
                for (int i = 0; i < nnx; ++i) {
                    if (
                        (this->PEx != this->mx - 1 && i == nnx - 1) || 
                        (this->PEy != this->my - 1 && j == nny - 1) || 
                        (this->PEz != this->mz - 1 && k == nnz - 1)
                    ) {
                        this->fout << "\t\t\t\t\t";
                        for (int c = 0; c < cmp; ++c) {
                            this->fout << 0 << " ";
                        }
                        this->fout << std::endl;
                    } else {
                        this->fout << "\t\t\t\t\t";
                        this->Values(i, j, k, _values...);
                        this->fout << std::endl;
                    }
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        this->AddFooter();
    }
}