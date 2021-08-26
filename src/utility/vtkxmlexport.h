#pragma once
#include <string>
#include <fstream>

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

        template<class F>
        void AddPointScaler(std::string _label, F _value);
        template<class F0, class F1, class F2>
        void AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2);
        template<class F00, class F01, class F02, class F10, class F11, class F12, class F20, class F21, class F22>
        void AddPointTensor(std::string _label, F00 _v00, F01 _v01, F02 _v02, F10 _v10, F11 _v11, F12 _v12, F20 _v20, F21 _v21, F22 _v22);

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
    };

    template<class F>
    void VTKXMLExport::AddPointScaler(std::string _label, F _value) {
        this->AddHeader();

        //  Export .pvts file
        if (this->PEid == 0) {
            this->fout0.seekp(this->addpos0);
            this->fout0 << "\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" format=\"ascii\"/>" << std::endl;
            this->addpos0 = this->fout0.tellp();
        }

        //  Export .vts file
        int nnx = this->nx + (this->PEx != this->mx - 1 ? 1 : 0);
        int nny = this->ny + (this->PEy != this->my - 1 ? 1 : 0);
        int nnz = this->nz + (this->PEz != this->mz - 1 ? 1 : 0);
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < nnz; ++k) {
            for (int j = 0; j < nny; ++j) {
                for (int i = 0; i < nnx; ++i) {
                    if (
                        (this->PEx != this->mx - 1 && i == nnx - 1) || 
                        (this->PEy != this->my - 1 && j == nny - 1) || 
                        (this->PEz != this->mz - 1 && k == nnz - 1)
                    ) {
                        this->fout << "\t\t\t\t\t" << 0 << std::endl;
                    } else {
                        this->fout << "\t\t\t\t\t" << _value(i, j, k) << std::endl;
                    }
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        this->AddFooter();
    }

    template<class F0, class F1, class F2>
    void VTKXMLExport::AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2) {
        this->AddHeader();

        //  Export .pvts file
        if (this->PEid == 0) {
            this->fout0.seekp(this->addpos0);
            this->fout0 << "\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
            this->addpos0 = this->fout0.tellp();
        }

        //  Export .vts file
        int nnx = this->nx + (this->PEx != this->mx - 1 ? 1 : 0);
        int nny = this->ny + (this->PEy != this->my - 1 ? 1 : 0);
        int nnz = this->nz + (this->PEz != this->mz - 1 ? 1 : 0);
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < nnz; ++k) {
            for (int j = 0; j < nny; ++j) {
                for (int i = 0; i < nnx; ++i) {
                    if (
                        (this->PEx != this->mx - 1 && i == nnx - 1) || 
                        (this->PEy != this->my - 1 && j == nny - 1) || 
                        (this->PEz != this->mz - 1 && k == nnz - 1)
                    ) {
                        this->fout << "\t\t\t\t\t" << 0 << " " << 0 << " " << 0 << std::endl;
                    } else {
                        this->fout << "\t\t\t\t\t" << _value0(i, j, k) << " " << _value1(i, j, k) << " " << _value2(i, j, k) << std::endl;
                    }
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        this->AddFooter();
    }

    template<class F00, class F01, class F02, class F10, class F11, class F12, class F20, class F21, class F22>
    void VTKXMLExport::AddPointTensor(std::string _label, F00 _v00, F01 _v01, F02 _v02, F10 _v10, F11 _v11, F12 _v12, F20 _v20, F21 _v21, F22 _v22) {
        this->AddHeader();

        //  Export .pvts file
        if (this->PEid == 0) {
            this->fout0.seekp(this->addpos0);
            this->fout0 << "\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"9\" format=\"ascii\"/>" << std::endl;
            this->addpos0 = this->fout0.tellp();
        }

        //  Export .vts file
        int nnx = this->nx + (this->PEx != this->mx - 1 ? 1 : 0);
        int nny = this->ny + (this->PEy != this->my - 1 ? 1 : 0);
        int nnz = this->nz + (this->PEz != this->mz - 1 ? 1 : 0);
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"9\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < nnz; ++k) {
            for (int j = 0; j < nny; ++j) {
                for (int i = 0; i < nnx; ++i) {
                    if (
                        (this->PEx != this->mx - 1 && i == nnx - 1) || 
                        (this->PEy != this->my - 1 && j == nny - 1) || 
                        (this->PEz != this->mz - 1 && k == nnz - 1)
                    ) {
                        this->fout << "\t\t\t\t\t" << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                    } else {
                        this->fout << "\t\t\t\t\t" 
                            << _v00(i, j, k) << " " << _v01(i, j, k) << " " << _v02(i, j, k) << " "
                            << _v10(i, j, k) << " " << _v11(i, j, k) << " " << _v12(i, j, k) << " "
                            << _v20(i, j, k) << " " << _v21(i, j, k) << " " << _v22(i, j, k) << std::endl;
                    }
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        this->AddFooter();
    }
}