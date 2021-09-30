#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <cassert>

namespace PANSLBM2 {
    template<class T>
    class VArray {
public:
        VArray() = delete;
        VArray(int _size, int _chunknum = 2) : 
            size(_size), chunknum(_chunknum), 
            chunksize((this->size + 1)/this->chunksize),
            chunkname(this->chunknum)
        ) {
            assert(this->size > 0 && this->chunknum > 1);
            this->data = new T[this->size];
            for (int chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                chunkname[chunkidx] = "chunk" + std::to_string(chunkidx);
            }
            this->chunkidxcurrent = 0;
        }
        VArray(const VArray<T>&) = delete;
        ~VArray() {
            delete[] this->data;
            for (int chunkidx = 0; chunkidx < this->chunknum; ++chunkidx) {
                remove(this->chunkname[chunkidx]);
            }
        }

        void Load(int _chunkidx);
        void LoadNext() {}
        void LoadPrev() {}

        const int size, chunknum, chunksize;
        T *data;
        std::vector<std::string> chunkname;
private:
        int chunkidxcurrent;
    };

    template<class T>
    void VArray<T>::Load(int _chunkidx) {
        assert(0 <= _chunkidx && _chunkidx < this->chunknum);
        
        //  Write current chunk
        std::ofstream fout(this->chunkname[this->chunkidxcurrent], ios::out|ios::binary|ios::trunc);
        if (!fout) {}
        for (int idx = 0; idx < this->chunksize; ++idx) {
            fout.write((char*)&this->data[idx], sizeof(T));
        }
        fout.close();

        //  Read idx chunk
        std::fin(this->chunkname[_chunkidx], ios::in|ios::binary);
        if (fin) {
            for (int idx = 0; !fin.eof(); ++idx) {
                fin.read((char*)&this->data[idx], sizeof(T));
            }
            fin.close();
        }
    }
}