#include "OtherFunctions.h"

#include <sstream>
#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <cstdlib>

#include <iterator>
#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

std::random_device rd;
boost::random::mt19937 gen(rd());




//##### INITIALIZE VECTORS

std::vector<int> uniformVectorInt(const int& sizeVec, const int& value){
    std::vector<int> vec; 
    vec.reserve(sizeVec);
    for (int i = 0; i < sizeVec; ++i) {
        vec.push_back(value);
    }
    return vec;
}
std::vector<double> uniformVectorDouble(const int& sizeVec, const double& value){
    std::vector<double> vec; 
    vec.reserve(sizeVec);
    for (int i = 0; i < sizeVec; ++i) {
        vec.push_back(value);
    }
    return vec;
}



//##### PRINT TERMINAL

void print(std::vector<int> vec){
    for (auto iter : vec) {
        std::cout << iter << " ";
    }
    std::cout << std::endl << "------" << std::endl;
}
void print(std::vector<double> vec){
    for (auto iter : vec) {
        std::cout << iter << " ";
    }
    std::cout << std::endl << "------" << std::endl;
}
void print(std::vector<std::array<double,3>> vec){
    for (size_t i(0); i<3; ++i) {
        for (auto iter : vec) {
            std::cout << iter[i] << " ";
        }
        std::cout << std::endl;

    }
    std::cout << "------" << std::endl;
}
void print(std::vector<std::vector<double>> vec){
    for (auto iter : vec) {
        for (auto iter2 : iter) {
            std::cout << iter2 << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "------" << std::endl;
}
void print(double nb){
    std::cout << nb << std::endl;
}
void print(int nb){
    std::cout << nb << std::endl;
}
void print(size_t nb){
    std::cout << nb << std::endl;
}

void print(std::string str, std::vector<int> vec){
    std::cout << str << " : ";        
    for (auto iter : vec) {
        std::cout << iter << " ";
    }
    std::cout << std::endl;
}
void print(std::string str, std::vector<double> vec){
    std::cout << str << " : ";        
    for (auto iter : vec) {
        std::cout << iter << " ";
    }
    std::cout << std::endl;
}
void print(std::string str, std::vector<std::array<double,3>> vec){
    std::cout << str << " : " << std::endl;
    for (size_t i(0); i<3; ++i) {
        for (auto iter : vec) {
            std::cout << iter[i] << " ";
        }
        std::cout << std::endl;

    }
}
void print(std::string str, std::vector<std::vector<double>> vec){
    std::cout << str << " : " << std::endl;
    for (auto iter : vec) {
        for (auto iter2 : iter) {
            std::cout << iter2 << " ";
        }
        std::cout << std::endl;
    }
}
void print(std::string str, double nb){
    std::cout << str << " : " << nb << std::endl;
}
void print(std::string str, int nb){
    std::cout << str << " : " << nb << std::endl;
}
void print(std::string str, size_t nb){
    std::cout << str << " : " << nb << std::endl;
}
void print(std::string str){
    std::cout << str << std::endl;
}



//##### BASIC STATS

double mean(const std::vector<double>& vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}
double mean(const std::vector<int>& vec){
    return std::accumulate(vec.begin(), vec.end(), 0.0) / vec.size();
}
double sd(const std::vector<int>& vec){
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0) / vec.size() - mean(vec) * mean(vec));
}
double sd(const std::vector<double>& vec){
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0) / vec.size() - mean(vec) * mean(vec));
}



//##### RANDOM SAMPLING

double randomDouble(){
//     return ((double) rand() / (RAND_MAX));           // THIS LEADS TO A SLIGHT COMPUTING NUMERICAL ERROR; THE REASON IS UNKNOWN.... INCOMPATIBLE WITH RANDOM NUMBER GENERATOR?
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}
int nrand(int min, int max){   // [min, max]
    std::uniform_int_distribution<int> uni(min,max);
    return uni(gen);
}

std::vector<int> sampleFromDiscreteDistributionWithReplacement(const int& sampleSize, std::vector<double> weights){
        //New
//     // sum all weights    
//     double sum_of_weight = 0;
//     for(size_t i=0; i<weights.size(); i++) {
//         sum_of_weight += weights[i];
//     }
//     
//     // sample random numbers
//     std::vector<double> samplesRandWeight; 
//     for (int sampleIndex = 0; sampleIndex < sampleSize; sampleIndex++) {
//         samplesRandWeight.push_back(randomDouble()*sum_of_weight);
//     }
//     
//     std::vector<int> samples = uniformVectorInt(sampleSize, -1); 
//     int nbFound = 0;
//     int indexWeight = 0;
//     while (nbFound<sampleSize && indexWeight<weights.size()){
//         for (size_t j=0; j<samplesRandWeight.size(); j++) {
//             if (samples[j]==-1) {
//                 if(samplesRandWeight[j] < weights[indexWeight]) {
//                     samples[j] = indexWeight;
//                     ++nbFound;
//                 }
//                 samplesRandWeight[j] -= weights[indexWeight];
//             }
//         }
//         ++indexWeight;
//     }
    
    //New
//     // sum all weights    
//     double sum_of_weight = 0;
//     for(size_t i=0; i<weights.size(); i++) {
//         sum_of_weight += weights[i];
//     }
//     
//     // sample random numbers
//     std::vector<double> samplesRandWeight; 
//     for (int sampleIndex = 0; sampleIndex < sampleSize; sampleIndex++) {
//         samplesRandWeight.push_back(randomDouble()*sum_of_weight);
//     }
//     
//     std::vector<int> samples = uniformVectorInt(sampleSize, -1); 
//     for(int i=0; i<weights.size(); i++) {
//         for (size_t j=0; j<samplesRandWeight.size(); j++) {
//             if (samples[j]==-1) {
//                 if(samplesRandWeight[j] < weights[i])
//                     samples[j] = i;
//                 samplesRandWeight[j] -= weights[i];
//             }
//         }
//     }
//     
    //Old
    // total size vector taken from weights vector
    std::vector<int> samples; 
    samples.reserve(sampleSize);
    std::discrete_distribution<int> distribution(weights.begin(), weights.end()); 
    for (auto iter = 0; iter < sampleSize; iter++) {
        samples.push_back(distribution(gen));
    }

        
    return samples;
}

// std::vector<int> sampleFromDiscreteDistribution(const int& sampleSize, std::vector<double> weights){
//     // total size vector taken from weights vector
//     const size_t totalSize = weights.size();
//     std::vector<int> samples; 
//     
//     if ((int)totalSize <= sampleSize) {   // sample size is high; all numbers are picked
//         
//         samples.reserve(totalSize);
//         for (int i=0; i<totalSize; ++i) {
//             samples.push_back(i);
//         }
//         
//     }else{ // sample size is small; a subset is sample without replacement
// 
//         // from :
//         // https://stackoverflow.com/questions/53632441/c-sampling-from-discrete-distribution-without-replacement
//         // see also
//         // https://stackoverflow.com/questions/48591592/how-do-i-use-boost-random
//         
//         // sudo apt-get install libboost-all-dev
//         // g++ -std=c++14 -w testOtherFunctions.cpp OtherFunctions.cpp 
// 
//         boost::random::uniform_01<> dist;
//         std::vector<double> vals;
//         for (auto iter : weights) {
//             vals.push_back(std::pow(dist(gen), 1. / iter));
//         }
//         // Sorting vals, but retain the indices. 
//         // There is unfortunately no easy way to do this with STL.
//         std::vector<std::pair<int, double>> valsWithIndices;
//         for (size_t iter = 0; iter < vals.size(); iter++) {
//             valsWithIndices.emplace_back(iter, vals[iter]);
//         }
//         std::sort(valsWithIndices.begin(), valsWithIndices.end(), [](auto x, auto y) {return x.second > y.second; });
// 
//         for (auto iter = 0; iter < sampleSize; iter++) {
//             samples.push_back(valsWithIndices[iter].first);
//         }     
//     }
//     return samples;
// }

std::vector<int> sampleFromNormalDiscreteDistribution(double mean, double sd, int size){
    std::normal_distribution<float> noise(mean, sd);
    std::vector<int> samples;
    samples.reserve(size);
    for (int i=0; i<size; ++i) {
        samples.push_back((int) round(noise(gen)));
    }
    return samples;
}

std::vector<double> sampleFromNormalDistribution(double mean, double sd, int size){
    std::normal_distribution<float> noise(mean, sd);
    std::vector<double> samples;
    samples.reserve(size);
    for (int i=0; i<size; ++i) {
        samples.push_back((double) noise(gen));
    }
    return samples;
}


//##### FUNCTIONS INDIVIDUALS

// std::vector<int>






// int main()
// {   
//     
//     testRunTime();
// 
//     return 0;
// }
