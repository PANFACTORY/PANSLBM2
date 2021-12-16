#pragma once
#include <string>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <vector>

namespace PANSLBM2 {
    class VTKXMLImport {
public:
        VTKXMLImport() = delete;
        VTKXMLImport(std::string _fname) : fname_vts(0) {
            int pos_lastsep = std::max((int)_fname.rfind('\\'), (int)_fname.rfind('/'));
            this->fpath = pos_lastsep != std::string::npos ? _fname.substr(0, pos_lastsep + 1) : "";
            
            this->fin_pvts.open(_fname);
            std::string tmp;
            std::string::size_type pos_keyword;
            while (std::getline(this->fin_pvts, tmp)) {
                //  全体の格子点数を読み込み
                pos_keyword = tmp.find("WholeExtent");
                if (pos_keyword != std::string::npos) {
                    std::string::size_type pos_dquot1 = tmp.find("\"", pos_keyword), pos_dquot2 = tmp.find("\"", pos_dquot1 + 1);
                    tmp = tmp.substr(pos_dquot1 + 1, pos_dquot2 - pos_dquot1 - 1);
                    std::istringstream tmpstream(tmp);
                    int nxmin, nxmax, nymin, nymax, nzmin, nzmax;
                    tmpstream >> nxmin >> nxmax >> nymin >> nymax >> nzmin >> nzmax;
                    this->nx = nxmax - nxmin;
                    this->ny = nymax - nymin;
                    this->nz = nzmax - nzmin;
                }

                //  vtsファイル名を読み込み
                pos_keyword = tmp.find("Source");
                if (pos_keyword != std::string::npos) {
                    std::string::size_type pos_dquot1 = tmp.find("\"", pos_keyword), pos_dquot2 = tmp.find("\"", pos_dquot1 + 1);
                    this->fname_vts.push_back(tmp.substr(pos_dquot1 + 1, pos_dquot2 - pos_dquot1 - 1));
                }
            }
        }
        VTKXMLImport(const VTKXMLImport&) = delete;
        ~VTKXMLImport() {}

        int GetNx() {   return this->nx;    };
        int GetNy() {   return this->ny;    };
        int GetNz() {   return this->nz;    };

        template<class T, template<class>class P>
        void GetPointScalar(P<T>& _p, std::string _label, T *_v);
        template<class T, template<class>class P>
        void GetPointVector(P<T>& _p, std::string _label, T *_vx, T *_vy);
        template<class T, template<class>class P>
        void GetPointVector(P<T>& _p, std::string _label, T *_vx, T *_vy, T *_vz);
        template<class T, template<class>class P>
        void GetPointTensor(P<T>& _p, std::string _label, T *_vxx, T *_vxy, T *_vyx, T *_vyy);
        template<class T, template<class>class P>
        void GetPointTensor(P<T>& _p, std::string _label, T *_vxx, T *_vxy, T *_vxz, T *_vyx, T *_vyy, T *_vyz, T *_vzx, T *_vzy, T *_vzz);

private:
        std::ifstream fin_pvts;
        int nx, ny, nz;
        std::string fpath;
        std::vector<std::string> fname_vts;

        template<class F>
        void GetPointBase(std::string _label, F _f);
    };

    template<class F>
    void VTKXMLImport::GetPointBase(std::string _label, F _f) {
        for (const auto& fname : this->fname_vts) {
            std::ifstream fin_vts(this->fpath + fname);
            std::string tmp;
            std::string::size_type pos_keyword;
            int nxmin, nxmax, nymin, nymax, nzmin, nzmax;
            while (std::getline(fin_vts, tmp)) {
                //  PE毎の格子点数を取得
                pos_keyword = tmp.find("Piece Extent");
                if (pos_keyword != std::string::npos) {
                    std::string::size_type pos_dquot1 = tmp.find("\"", pos_keyword), pos_dquot2 = tmp.find("\"", pos_dquot1 + 1);
                    tmp = tmp.substr(pos_dquot1 + 1, pos_dquot2 - pos_dquot1 - 1);
                    std::istringstream tmpstream(tmp);
                    tmpstream >> nxmin >> nxmax >> nymin >> nymax >> nzmin >> nzmax;
                }

                //  PE毎の値を読み込み
                pos_keyword = tmp.find("Name=\"" + _label + "\"");
                if (pos_keyword != std::string::npos) {
                    std::cout << "find out: " << _label << " in " << this->fpath + fname << std::endl;
                    for (int k = nzmin; k <= nzmax; ++k) {
                        for (int j = nymin; j <= nymax; ++j) {
                            for (int i = nxmin; i <= nxmax; ++i) {
                                _f(i, j, k, fin_vts);
                            }
                        }    
                    }
                }
            }
        }
    }

    template<class T, template<class>class P>
    void VTKXMLImport::GetPointScalar(P<T>& _p, std::string _label, T *_v) {
        this->GetPointBase(_label, [&](int _i, int _j, int _k, std::ifstream& _fin_vts) {
            int idx = _p.Index(_i, _j, _k);
            _fin_vts >> _v[idx];
        });
    }

    template<class T, template<class>class P>
    void VTKXMLImport::GetPointVector(P<T>& _p, std::string _label, T *_vx, T *_vy) {
        T tmp;
        this->GetPointBase(_label, [&](int _i, int _j, int _k, std::ifstream& _fin_vts) {
            int idx = _p.Index(_i, _j, _k);
            _fin_vts >> _vx[idx] >> _vy[idx] >> tmp;
        });
    }

    template<class T, template<class>class P>
    void VTKXMLImport::GetPointVector(P<T>& _p, std::string _label, T *_vx, T *_vy, T *_vz) {
        this->GetPointBase(_label, [&](int _i, int _j, int _k, std::ifstream& _fin_vts) {
            int idx = _p.Index(_i, _j, _k);
            _fin_vts >> _vx[idx] >> _vy[idx] >> _vz[idx];
        });
    }

    template<class T, template<class>class P>
    void VTKXMLImport::GetPointTensor(P<T>& _p, std::string _label, T *_vxx, T *_vxy, T *_vyx, T *_vyy) {
        T tmp;
        this->GetPointBase(_label, [&](int _i, int _j, int _k, std::ifstream& _fin_vts) {
            int idx = _p.Index(_i, _j, _k);
            _fin_vts >> _vxx[idx] >> _vxy[idx] >> tmp >> _vyx[idx] >> _vyy[idx] >> tmp >> tmp >> tmp >> tmp;
        });
    }

    template<class T, template<class>class P>
    void VTKXMLImport::GetPointTensor(P<T>& _p, std::string _label, T *_vxx, T *_vxy, T *_vxz, T *_vyx, T *_vyy, T *_vyz, T *_vzx, T *_vzy, T *_vzz) {
        this->GetPointBase(_label, [&](int _i, int _j, int _k, std::ifstream& _fin_vts) {
            int idx = _p.Index(_i, _j, _k);
            _fin_vts >> _vxx[idx] >> _vxy[idx] >> _vxz[idx] >> _vyx[idx] >> _vyy[idx] >> _vyz[idx] >> _vzx[idx] >> _vzy[idx] >> _vzz[idx];
        });
    }
}