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
        
private:
        std::ifstream fin;
        int nx, ny, nz, np, nc;
    };

    template<class T>
    void VTKImport::GetPointScaler(std::string _label, T *_p) {
        this->fin.clear();
        this->fin.seekg(0, std::ios_base::beg); 
        std::string tmp;
        while (this->fin >> tmp) {
            if (tmp == "POINT_DATA") {
                while (this->fin >> tmp) {
                    if (tmp == "SCALARS") {
                        this->fin >> tmp;
                        if (tmp == _label) {
                            this->fin >> tmp >> tmp >> tmp;
                            for (int i = 0; i < this->np; i++) {
                                this->fin >> _p[i];
                            }
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
    void VTKImport::GetPointVector(std::string _label, T *_px, T *_py) {
        this->fin.clear();
        this->fin.seekg(0, std::ios_base::beg); 
        std::string tmp;
        while (this->fin >> tmp) {
            if (tmp == "POINT_DATA") {
                while (this->fin >> tmp) {
                    if (tmp == "VECTORS") {
                        this->fin >> tmp;
                        if (tmp == _label) {
                            this->fin >> tmp;
                            for (int i = 0; i < this->np; i++) {
                                this->fin >> _px[i] >> _py[i] >> tmp;
                            }
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
}