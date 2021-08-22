#pragma once
#include <string>
#include <fstream>

namespace PANSLBM2 {
    //  VTK format file export class (single process version)
    class VTKExport {
public:
        VTKExport() = delete;
        VTKExport(std::string _fname, int _nx, int _ny = 1, int _nz = 1, std::string _title = "notitle") : nx(_nx), ny(_ny), nz(_nz) {
            this->fout.open(_fname);
            this->isaddpointdata = false;

            //  Add vtk header  
            this->fout << "# vtk DataFile Version 3.0" << std::endl;
            this->fout << _title << std::endl;
            this->fout << "ASCII" << std::endl;
            this->fout << "DATASET\tSTRUCTURED_GRID" << std::endl;
            this->fout << "DIMENSIONS\t" << this->nx << "\t" << this->ny << "\t" << this->nz << std::endl;
        
            //  Add point coordinates
            this->fout << "POINTS\t" << this->nx*this->ny*this->nz << "\t" << "float" << std::endl;
            for (int k = 0; k < this->nz; k++) {
                for (int j = 0; j < this->ny; j++) {
                    for (int i = 0; i < this->nx; i++) {
                        fout << i << "\t" << j << "\t" << k << std::endl;
                    }
                }    
            }
        }
        VTKExport(const VTKExport&) = delete;
        ~VTKExport() {}

        template<class F>
        void AddPointScaler(std::string _label, F _value) {
            //  Add point data header
            if (!this->isaddpointdata) {
                this->fout << "POINT_DATA\t" << this->nx*this->ny*this->nz << std::endl;
                this->isaddpointdata = true;
            }

            //  Add point data value
            this->fout << "SCALARS\t" << _label << "\tfloat" << std::endl;
            this->fout << "LOOKUP_TABLE\tdefault" << std::endl;
            for (int k = 0; k < this->nz; k++) {
                for (int j = 0; j < this->ny; j++) {
                    for (int i = 0; i < this->nx; i++) {
                        this->fout << _value(i, j, k) << std::endl;
                    }
                }
            }
        }
        template<class F0, class F1, class F2>
        void AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2) {
            //  Add point data header
            if (!this->isaddpointdata) {
                this->fout << "POINT_DATA\t" << this->nx*this->ny*this->nz << std::endl;
                this->isaddpointdata = true;
            }

            //  Add point data value
            this->fout << "VECTORS\t" << _label << "\tfloat" << std::endl;
            for (int k = 0; k < this->nz; k++) {
                for (int j = 0; j < this->ny; j++) {
                    for (int i = 0; i < this->nx; i++) {
                        this->fout << _value0(i, j, k) << "\t" << _value1(i, j, k) << "\t" << _value2(i, j, k) << std::endl;
                    }
                }
            }
        }
        template<class F00, class F01, class F02, class F10, class F11, class F12, class F20, class F21, class F22>
        void AddPointTensor(std::string _label, F00 _v00, F01 _v01, F02 _v02, F10 _v10, F11 _v11, F12 _v12, F20 _v20, F21 _v21, F22 _v22) {
            //  Add point data header
            if (!this->isaddpointdata) {
                this->fout << "POINT_DATA\t" << this->nx*this->ny*this->nz << std::endl;
                this->isaddpointdata = true;
            }

            //  Add point data value
            this->fout << "TENSORS\t" << _label << "\tfloat" << std::endl;
            for (int k = 0; k < this->nz; k++) {
                for (int j = 0; j < this->ny; j++) {
                    for (int i = 0; i < this->nx; i++) {
                        this->fout << _v00(i, j, k) << "\t" << _v01(i, j, k) << "\t" << _v02(i, j, k) << std::endl;
                        this->fout << _v10(i, j, k) << "\t" << _v11(i, j, k) << "\t" << _v12(i, j, k) << std::endl;
                        this->fout << _v20(i, j, k) << "\t" << _v21(i, j, k) << "\t" << _v22(i, j, k) << std::endl << std::endl;
                    }
                }
            }
        }

private:
        std::ofstream fout;
        const int nx, ny, nz;
        bool isaddpointdata;
    };

    //  VTK(XML) format file export class (multi process version)
    class VTKXMLExport {
public:
        VTKXMLExport() = delete;
        VTKXMLExport(std::string _fname, int _PEid, int _PETot, int _lx, int _ly = 1, int _lz = 1) : lx(_nx), ly(_ny), lz(_nz) {
            //  Make .pvts file
            if (_PEid == 0) {
                std::ofstream fout0(_fname);
                fout0 << "<VTKFile type=\"PStructuredGrid\">" << std::endl;
                fout0 << "\t<PStructuredGrid WholeExtent=\"0 " << this->lx << " 0 " << this->ly << " 0 " << this->lz << "\">" << std::endl;
                fout0 << "\t\t<PPoints>" << std::endl;
                fout0 << "\t\t\t<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\"/>" << std::endl;
                fout0 << "\t\t</PPoints>" << std::endl;
                for (int peid = 0; peid < _PETot; ++peid) {
                    fout << "\t\t<Piece Source=\"" << _fname + "_" + std::to_string(peid) + ".vts" << "\"/>" << std::endl;
                }
                fout0 << "\t</PStructuredGrid>" << std::endl;
                fout0 << "</VTKFile>";
            }
        }
        VTKXMLExport(const VTKXMLExport&) = delete;
        ~VTKXMLExport() {}

private:
        const int lx, ly, lz;
    };
}