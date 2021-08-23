#pragma once
#include <string>
#include <fstream>

namespace PANSLBM2 {
    class VTKXMLExport {
public:
        VTKXMLExport() = delete;
        VTKXMLExport(std::string _fname, int _PEid, int _lx, int _ly, int _lz, int _mx, int _my, int _mz) :
            lx(_lx), ly(_ly), lz(_lz), mx(_mx), my(_my), mz(_mz), PEid(_PEid)
            PEx(this->PEid%this->mx), PEy((this->PEid/this->mx)%this->my), PEz(this->PEid/(this->mx*this->my)),
            nx(((this->lx - 1) + this->PEx)/this->mx + 1), ny(((this->ly - 1) + this->PEy)/this->my + 1), nz(((this->lz - 1) + this->PEz)/this->mz + 1),
            offsetx(this->mx - this->PEx > (this->lx - 1)%this->mx ? this->PEx*(this->nx - 1) : (this->lx - 1) - (this->mx - this->PEx)*(this->nx - 1)),
            offsety(this->my - this->PEy > (this->ly - 1)%this->my ? this->PEy*(this->ny - 1) : (this->ly - 1) - (this->my - this->PEy)*(this->ny - 1)),
            offsetz(this->mz - this->PEz > (this->lz - 1)%this->mz ? this->PEz*(this->nz - 1) : (this->lz - 1) - (this->mz - this->PEz)*(this->nz - 1))
        {
            //  Export .pvts file (Only host process)
            if (this->PEid == 0) {
                std::ofstream fout0(_fname + ".pvts");
                fout0 << "<VTKFile type=\"PStructuredGrid\">" << std::endl;
                fout0 << "\t<PStructuredGrid GhostLevel=\"0\" WholeExtent=\"0 " << this->lx - 1 << " 0 " << this->ly - 1 << " 0 " << this->lz - 1 << "\">" << std::endl;
                fout0 << "\t\t<PPoints>" << std::endl;
                fout0 << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
                fout0 << "\t\t</PPoints>" << std::endl;
                for (int k = 0; k < this->mz; ++k) {
                    for (int j = 0; j < this->my; ++j) {
                        for (int i = 0; i < this->mx; ++i) {
                            int peid = i + this->mx*j + this->mx*this->my*k;
                            fout0 << "\t\t<Piece Extent=\"" << 0 << " " << 3 << " " << 0 << " " << 2 << " " << 0 << " " << 0 
                                << "\" Source=\"" << _fname + "_" + std::to_string(peid) + ".vts" << "\"/>" << std::endl;
                        }
                    }
                }
                fout0 << "\t</PStructuredGrid>" << std::endl;
                fout0 << "</VTKFile>";
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
            this->fout << "\t\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
            for (int k = 0; k < this->nz; ++k) {
                for (int j = 0; j < this->ny; ++j) {
                    for (int i = 0; i < this->nx; ++i) {
                        this->fout << "\t\t\t\t\t" << i + this->offsetx << " " << j + this->offsety << " " << k + this->offsetz << std::endl;
                    }
                }
            }
            this->fout << "\t\t\t\t</DataArray>" << std::endl;
            this->fout << "\t\t\t</Points>" << std::endl;
            this->fout << "\t\t</Piece>" << std::endl;
            this->fout << "\t</StructuredGrid>" << std::endl;
            this->fout << "</VTKFile>";
        }
        VTKXMLExport(const VTKXMLExport&) = delete;
        ~VTKXMLExport() {}

private:
        std::ofstream fout;
        const int lx, ly, lz, mx, my, mz, PEid, PEx, PEy, PEz, nx, ny, nz, offsetx, offsety, offsetz;
    }
}