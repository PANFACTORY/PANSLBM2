#pragma once
#include <vector>


namespace PANSLBM2{
    //**********Heaviside filter class**********
    template<class T>
    class HeavisideFilter{
public:
        HeavisideFilter();
        ~HeavisideFilter();
        HeavisideFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w);
    

        void UpdateBeta(T _beta);
        std::vector<T> GetFilteredVariables(std::vector<T> _s);            //  Return filtered variables
        std::vector<T> GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho);    //  Return filtered sensitivities

    
private:
        const int n;                                //  Number of design variables 
        T beta;                                     //  Parameter of Heaviside filter
        std::vector<std::vector<int> > neighbors;   //  Index list of neighbor element
        std::vector<std::vector<T> > w;             //  Weight list of neighbor element    
    };


    template<class T>
    HeavisideFilter<T>::HeavisideFilter() : n(0){
        this->beta = 1.0;
    }


    template<class T>
    HeavisideFilter<T>::~HeavisideFilter(){}


    template<class T>
    HeavisideFilter<T>::HeavisideFilter(int _n, std::vector<std::vector<int> > _neighbors, std::vector<std::vector<T> > _w) : n(_n){
        this->beta = 1.0;
        this->neighbors = _neighbors;
        this->w = _w;
    }


    template<class T>
    void HeavisideFilter<T>::UpdateBeta(T _beta){
        this->beta = _beta;
    }


    template<class T>
    std::vector<T> HeavisideFilter<T>::GetFilteredVariables(std::vector<T> _s){
        std::vector<T> rho = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wssum = T();
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                wssum += this->w[i][j]*_s[this->neighbors[i][j]];
                wsum += this->w[i][j];
            }
            rho[i] = 0.5*(tanh(0.5*this->beta) + tanh(this->beta*(wssum/wsum - 0.5)))/tanh(0.5*this->beta);
        }
        return rho;
    }


    template<class T>
    std::vector<T> HeavisideFilter<T>::GetFilteredSensitivitis(std::vector<T> _s, std::vector<T> _dfdrho){
        std::vector<T> drhodstilde = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wssum = T();
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                wssum += this->w[i][j]*_s[this->neighbors[i][j]];
                wsum += this->w[i][j];
            }
            drhodstilde[i] = 0.5*this->beta*(1.0 - pow(tanh(this->beta*(wssum/wsum - 0.5)), 2.0))/tanh(0.5*this->beta);
        }
        
        std::vector<T> dfds = std::vector<T>(this->n);
        for(int i = 0; i < this->n; i++){
            T wsum = T();
            for(int j = 0; j < this->neighbors[i].size(); j++){
                dfds[i] += _dfdrho[this->neighbors[i][j]]*drhodstilde[this->neighbors[i][j]]*this->w[i][j];
                wsum += this->w[i][j];
            }
            dfds[i] /= wsum;       
        }
        return dfds;
    }
}