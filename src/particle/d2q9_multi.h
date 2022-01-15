#pragma once
#include <vector>

namespace PANSLBM2 {
    template<class T>
    class D2Q9_multi {
public:
        D2Q9_multi() = delete;
        D2Q9_multi(int _lx, int _ly) {

        }
        D2Q9_multi(const D2Q9_multi<T>&) = delete;
        ~D2Q9_multi() {}

        std::vector<D2Q9_multi<T>* > children;

        void GenerateChild(int _lx, int _ly) {
            this->children.push_back(new D2Q9_multi<T>(_lx, _ly));
        }
    };
}