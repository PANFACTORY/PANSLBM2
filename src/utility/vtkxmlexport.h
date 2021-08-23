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
            nx(((this->lx - 1) + this->PEx)/this->mx + 1), ny(((this->ly - 1) + this->PEy)/this->my + 1), nz(((this->lz - 1) + this->PEz)/this->mz + 1),
            offsetx(this->mx - this->PEx > (this->lx - 1)%this->mx ? this->PEx*(this->nx - 1) : (this->lx - 1) - (this->mx - this->PEx)*(this->nx - 1)),
            offsety(this->my - this->PEy > (this->ly - 1)%this->my ? this->PEy*(this->ny - 1) : (this->ly - 1) - (this->my - this->PEy)*(this->ny - 1)),
            offsetz(this->mz - this->PEz > (this->lz - 1)%this->mz ? this->PEz*(this->nz - 1) : (this->lz - 1) - (this->mz - this->PEz)*(this->nz - 1))
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
                            int nnx = ((this->lx - 1) + i)/this->mx + 1, nny = ((this->ly - 1) + j)/this->my + 1, nnz = ((this->lz - 1) + k)/this->mz + 1;
                            int osx = this->mx - i > (this->lx - 1)%this->mx ? i*(nnx - 1) : (this->lx - 1) - (this->mx - i)*(nnx - 1);
                            int osy = this->my - j > (this->ly - 1)%this->my ? j*(nny - 1) : (this->ly - 1) - (this->my - j)*(nny - 1);
                            int osz = this->mz - k > (this->lz - 1)%this->mz ? k*(nnz - 1) : (this->lz - 1) - (this->mz - k)*(nnz - 1);
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
            this->fout << "\t\t<Piece Extent=\"" 
                << 0 + this->offsetx << " " << (this->nx - 1) + this->offsetx << " " 
                << 0 + this->offsety << " " << (this->ny - 1) + this->offsety << " " 
                << 0 + this->offsetz << " " << (this->nz - 1) + this->offsetz << "\">" << std::endl;
            this->fout << "\t\t\t<Points>" << std::endl; 
            this->fout << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
            for (int k = 0; k < this->nz; ++k) {
                for (int j = 0; j < this->ny; ++j) {
                    for (int i = 0; i < this->nx; ++i) {
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
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    this->fout << "\t\t\t\t\t" << _value(i, j, k) << std::endl;
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
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    this->fout << "\t\t\t\t\t" << _value0(i, j, k) << " " << _value1(i, j, k) << " " << _value2(i, j, k) << std::endl;
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
        this->fout.seekp(this->addpos);
        this->fout << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << _label << "\" NumberOfComponents=\"9\" format=\"ascii\">" << std::endl;
        for (int k = 0; k < this->nz; ++k) {
            for (int j = 0; j < this->ny; ++j) {
                for (int i = 0; i < this->nx; ++i) {
                    this->fout << "\t\t\t\t\t" 
                        << _value00(i, j, k) << " " << _value01(i, j, k) << " " << _value02(i, j, k) << " "
                        << _value10(i, j, k) << " " << _value11(i, j, k) << " " << _value12(i, j, k) << " "
                        << _value20(i, j, k) << " " << _value21(i, j, k) << " " << _value22(i, j, k) << std::endl;
                }
            }
        }
        this->fout << "\t\t\t\t</DataArray>" << std::endl;
        this->addpos = this->fout.tellp();

        this->AddFooter();
    }
}