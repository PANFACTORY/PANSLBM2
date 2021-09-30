#pragma once
#include <string>
#include <fstream>
#include <cassert>

namespace PANSLBM2 {
    template<class T>
    class VArray {
public:
        VArray() = delete;
        VArray(size_t _size, size_t _chunknum = 2) : 
            size(_size), chunknum(_chunknum), 
            chunksize((this->size + 1)/this->chunknum)
        {
            assert(this->size > 0 && this->chunknum > 1);
            this->data = new T[this->chunksize];
            this->chunkname = new std::string[this->chunknum];
            for (size_t chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                this->chunkname[chunkidx] = "varraychunk_" + std::to_string(VArray<T>::objectid) + "_" + std::to_string(chunkidx);
            }
            this->chunkidxcurrent = 0;
            VArray<T>::objectid++;
        }
        VArray(const VArray<T>&) = delete;
        ~VArray() {
            delete[] this->data;
            for (size_t chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                remove(this->chunkname[chunkidx].c_str());
            }
            delete[] this->chunkname;
        }

        T &operator[](size_t _idx) {
            size_t chunkidx = _idx/this->chunksize;
            if (chunkidx != this->chunkidxcurrent) {
                this->Swap(chunkidx);
            }
            return this->data[_idx%this->chunksize];
        }

        const size_t size, chunknum, chunksize;
        
private:
        T *data;
        std::string *chunkname;
        size_t chunkidxcurrent;
        static int objectid;
        
        void Swap(size_t _chunkidx);
    };

    template<class T>
    int VArray<T>::objectid = 0;

    template<class T>
    void VArray<T>::Swap(size_t _chunkidx) {
        assert(_chunkidx < this->chunknum);
        
        //  Write current chunk
        const int buffersize = 1000000;
        char buffer[buffersize];
        std::ofstream fout(this->chunkname[this->chunkidxcurrent], std::ios::out|std::ios::binary|std::ios::trunc);
        ofstream.rdbuf()->pubsetbuf(buffer, buffersize);
        if (!fout) {}
        for (size_t idx = 0; idx < this->chunksize; ++idx) {
            fout.write((char*)&this->data[idx], sizeof(T));
        }
        fout.close();

        //  Read idx chunk
        std::ifstream fin(this->chunkname[_chunkidx], std::ios::in|std::ios::binary);
        if (fin) {
            for (size_t idx = 0; !fin.eof(); ++idx) {
                fin.read((char*)&this->data[idx], sizeof(T));
            }
            fin.close();
        }

        this->chunkidxcurrent = _chunkidx;
    }
}