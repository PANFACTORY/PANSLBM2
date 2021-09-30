#pragma once
#include <string>
#include <fstream>
#include <typeinfo>
#include <cassert>

namespace PANSLBM2 {
    template<class T>
    class VArray {
public:
        VArray() = delete;
        VArray(int _size, int _chunknum = 2) : 
            size(_size), chunknum(_chunknum), 
            chunksize((this->size + 1)/this->chunknum)
        {
            assert(this->size > 0 && this->chunknum > 1);
            this->data = new T[this->chunksize];
            this->chunkname = new std::string[this->chunknum];
            for (int chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                this->chunkname[chunkidx] = "varraychunk_" + typeid(T).name() + "_" + std::to_string(VArray<T>::objectid) + "_" + std::to_string(chunkidx);
            }
            this->chunkidxcurrent = 0;
            VArray<T>::objectid++;
        }
        VArray(const VArray<T>&) = delete;
        ~VArray() {
            delete[] this->data;
            for (int chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                remove(this->chunkname[chunkidx].c_str());
            }
            delete[] this->chunkname;
        }

        T &operator[](int _idx) {
            int chunkidx = _idx/this->chunksize;
            if (chunkidx != this->chunkidxcurrent) {
                this->Swap(chunkidx);
            }
            return this->data[_idx%this->chunksize];
        }

        const int size, chunknum, chunksize;
        
private:
        T *data;
        std::string *chunkname;
        int chunkidxcurrent;
        static int objectid;
        
        void Swap(int _chunkidx);
    };

    template<class T>
    int VArray<T>::objectid = 0;

    template<class T>
    void VArray<T>::Swap(int _chunkidx) {
        assert(0 <= _chunkidx && _chunkidx < this->chunknum);
        
        //  Write current chunk
        std::ofstream fout(this->chunkname[this->chunkidxcurrent], std::ios::out|std::ios::binary|std::ios::trunc);
        if (!fout) {}
        for (int idx = 0; idx < this->chunksize; ++idx) {
            fout.write((char*)&this->data[idx], sizeof(T));
        }
        fout.close();

        //  Read idx chunk
        std::ifstream fin(this->chunkname[_chunkidx], std::ios::in|std::ios::binary);
        if (fin) {
            for (int idx = 0; !fin.eof(); ++idx) {
                fin.read((char*)&this->data[idx], sizeof(T));
            }
            fin.close();
        }

        this->chunkidxcurrent = _chunkidx;
    }
}