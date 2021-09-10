#pragma once
#include <vector>

namespace PANSLBM2 {
    template<class T>
    class DensityFilter {
public:
        DensityFilter() = delete;
        DensityFilter(int _n) : n(_n) {
            this->neighbors(this->n);
            this->w(this->n);
        }
        DensityFilter(const DensityFilter<T>&) = delete;
        ~DensityFilter() {};

        std::vector<T> GetFilteredVariables(const std::vector<T> &_s);                                      //  Return filtered variables
        std::vector<T> GetFilteredSensitivitis(const std::vector<T> &_s, const std::vector<T> &_dfdrho);    //  Return filtered sensitivities

        std::vector<std::vector<int> > neighbors;   //  Index list of neighbor element
        std::vector<std::vector<T> > w;             //  Weight list of neighbor element    

private:
        const int n;                                //  Number of design variables 
    };

    template<class T>
    std::vector<T> DensityFilter<T>::GetFilteredVariables(const std::vector<T> &_s){
        std::vector<T> rho(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                rho[i] += this->w[i][j]*_s[this->neighbors[i][j]];
                wsum += this->w[i][j];
            }
            rho[i] /= wsum;
        }
        return rho;
    }

    template<class T>
    std::vector<T> DensityFilter<T>::GetFilteredSensitivitis(const std::vector<T> &_s, const std::vector<T> &_dfdrho){
        std::vector<T> dfds(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                dfds[i] += _dfdrho[this->neighbors[i][j]]*this->w[i][j];
                wsum += this->w[i][j];
            }
            dfds[i] /= wsum;       
        }
        return dfds;
    }
}