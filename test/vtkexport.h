#pragma once
#include <string>
#include <fstream>

namespace PANSLBM2 {
    //*************************************************************************
    //  This class exports datas to a file with VTK format.
    //  The dataset format is Structured Grid.
    //  This class is used to storage LBM datas.
    //*************************************************************************
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
        void AddPointScaler(std::string _label, F _value);
        template<class F0, class F1, class F2>
        void AddPointVector(std::string _label, F0 _value0, F1 _value1, F2 _value2);
        template<class F00, class F01, class F02, class F10, class F11, class F12, class F20, class F21, class F22>
        void AddPointTensor(std::string _label, F00 _v00, F01 _v01, F02 _v02, F10 _v10, F11 _v11, F12 _v12, F20 _v20, F21 _v21, F22 _v22);

private:
        std::ofstream fout;
        const int nx, ny, nz;
        bool isaddpointdata;
    };

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

    template<class F00, class F01, class F02, class F10, class F11, class F12, class F20, class F21, class F22>
    void VTKExport::AddPointTensor(std::string _label, F00 _v00, F01 _v01, F02 _v02, F10 _v10, F11 _v11, F12 _v12, F20 _v20, F21 _v21, F22 _v22) {
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

    //*************************************************************************
    //  This class exports datas to a file with VTK format.
    //  The dataset format is Polygonal Data.
    //  This class is used to storage Immersed-Boundary datas.
    //*************************************************************************
    class VTKExportIB {
public:
        VTKExportIB() = delete;
        VTKExportIB(std::string _fname, std::string _title = "notitle") {
            this->fout.open(_fname);

            //  Add vtk header  
            this->fout << "# vtk DataFile Version 3.0" << std::endl;
            this->fout << _title << std::endl;
            this->fout << "ASCII" << std::endl;
            this->fout << "DATASET POLYDATA" << std::endl;
        }
        VTKExportIB(const VTKExportIB&) = delete;
        ~VTKExportIB() {}

        template<class T>
        void AddPoint(int _n, T *_x, T *_y = nullptr, T *_z = nullptr);

private:
        std::ofstream fout;
    };

    template<class T>
    void VTKExportIB::AddPoint(int _n, T *_x, T *_y, T *_z) {
        this->fout << "POINTS\t" << _n << "\tfloat" << std::endl;
        for (int i = 0; i < _n; i++) {
            this->fout << _x[i] << "\t" << (_y ? _y[i] : T()) << "\t" << (_z ? _z[i] : T()) << std::endl;
        }
        this->fout << "VERTICES\t" << _n << "\t" << 2*_n << std::endl;
        for (int i = 0; i < _n; i++) {
            this->fout << "1\t" << i << std::endl;
        }
    }
}