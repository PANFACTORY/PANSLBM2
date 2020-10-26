#pragma once
#include <string>
#include <fstream>

namespace PANSLBM2 {
    class VTKExport {
public:
        VTKExport() = delete;
        VTKExport(std::string _fname, int _nx, int _ny = 1, int _nz = 1, std::string _title = "notitle");
        VTKExport(const VTKExport&) = delete;
        ~VTKExport() {}

        template<class F>
        void AddPointScaler(std::string _label, F _value);
        template<class F0, class F1, class F2>
        void AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2);

private:
        std::ofstream fout;
        const int nx, ny, nz;
        bool isaddpointdata;
    };

    VTKExport::VTKExport(std::string _fname, int _nx, int _ny, int _nz, std::string _title) : nx(_nx), ny(_ny), nz(_nz) {
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

    template<class F>
    void VTKExport::AddPointScaler(std::string _label, F _value) {
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
    void VTKExport::AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2) {
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
}