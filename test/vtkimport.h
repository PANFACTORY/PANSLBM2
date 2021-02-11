#pragma once
#include <iostream>
#include <string>
#include <fstream>

namespace PANSLBM2 {
    class VTKImport {
public:
        VTKImport() = delete;
        VTKImport(std::string _fname) {
            this->fin.open(_fname);
            std::string tmp;
            while (this->fin >> tmp) {
                if (tmp == "DIMENSIONS") {
                    this->fin >> this->nx >> this->ny >> this->nz;
                } else if (tmp == "POINT_DATA") {
                    this->fin >> this->np;
                } else if (tmp == "CELL_DATA") {
                    this->fin >> this->nc;
                }
            }
        }
        VTKImport(const VTKImport&) = delete;
        ~VTKImport() {}

        template<class T>
        void GetPointScaler(std::string _label, T *_p);
        template<class T>
        void GetPointVector(std::string _label, T *_px, T *_py);
        template<class T>
        void GetPointVector(std::string _label, T *_px, T *_py, T *_pz);
        template<class T>
        void GetPointTensor(std::string _label, T *_pxx, T *_pxy, T *_pyx, T *_pyy);
        template<class T>
        void GetPointTensor(std::string _label, T *_pxx, T *_pxy, T *_pxz, T *_pyx, T *_pyy, T *_pyz, T *_pzx, T *_pzy, T *_pzz);
        
private:
        std::ifstream fin;
        int nx, ny, nz, np, nc;

        template<class F>
        void GetPointBase(std::string _label, std::string _type, F _f);
    };

    template<class F>
    void VTKImport::GetPointBase(std::string _label, std::string _type, F _f) {
        this->fin.clear();
        this->fin.seekg(0, std::ios_base::beg); 
        std::string tmp;
        while (this->fin >> tmp) {
            if (tmp == "POINT_DATA") {
                while (this->fin >> tmp) {
                    if (tmp == _type) {
                        this->fin >> tmp;
                        if (tmp == _label) {
                            _f();
                        }
                    } else if (tmp == "CELL_DATA") {
                        break;
                    }
                }
            } else if (tmp == "CELL_DATA") {
                break;
            }
        }
    }

    template<class T>
    void VTKImport::GetPointScaler(std::string _label, T *_p) {
        std::string tmp;
        this->GetPointBase(_label, "SCALERS", [&]() {
            this->fin >> tmp >> tmp >> tmp;
            for (int i = 0; i < this->np; i++) {
                this->fin >> _p[i];
            }
        });
    }

    template<class T>
    void VTKImport::GetPointVector(std::string _label, T *_px, T *_py) {
        std::string tmp;
        this->GetPointBase(_label, "VECTORS", [&]() {
            this->fin >> tmp;
            for (int i = 0; i < this->np; i++) {
                this->fin >> _px[i] >> _py[i] >> tmp;
            }
        });
    }

    template<class T>
    void VTKImport::GetPointVector(std::string _label, T *_px, T *_py, T *_pz) {
        std::string tmp;
        this->GetPointBase(_label, "VECTORS", [&]() {
            this->fin >> tmp;
            for (int i = 0; i < this->np; i++) {
                this->fin >> _px[i] >> _py[i] >> _pz[i];
            }
        });
    }

    template<class T>
    void VTKImport::GetPointTensor(std::string _label, T *_pxx, T *_pxy, T *_pyx, T *_pyy) {
        std::string tmp;
        this->GetPointBase(_label, "TENSORS", [&]() {
            this->fin >> tmp;
            for (int i = 0; i < this->np; i++) {
                this->fin >> _pxx[i] >> _pxy[i] >> tmp >> _pyx[i] >> _pyy[i] >> tmp >> tmp >> tmp >> tmp;
            }
        });
    }

    template<class T>
    void VTKImport::GetPointTensor(std::string _label, T *_pxx, T *_pxy, T *_pxz, T *_pyx, T *_pyy, T *_pyz, T *_pzx, T *_pzy, T *_pzz) {
        std::string tmp;
        this->GetPointBase(_label, "TENSORS", [&]() {
            this->fin >> tmp;
            for (int i = 0; i < this->np; i++) {
                this->fin >> _pxx[i] >> _pxy[i] >> _pxz[i] >> _pyx[i] >> _pyy[i] >> _pyz[i] >> _pzx[i] >> _pzy[i] >> _pzz[i];
            }
        });
    }
}