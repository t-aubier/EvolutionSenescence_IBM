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


// #### UNCOMMENT TO USE THE BOOST PACKAGE; comments also the lines annotated 'boost' below
// #include <boost/random/uniform_01.hpp>
// #include <boost/random/mersenne_twister.hpp>
// #include <boost/random/discrete_distribution.hpp>
// std::random_device rd;
// boost::random::mt19937 gen(rd());


// #### UNCOMMENT NOT TO USE THE BOOST PACKAGE
//following lines were commented and edited:
//std::random_device rd;
//std::mt19937 gen(rd());
//https://stackoverflow.com/questions/65788866/stdrandom-on-windows-always-i-get-the-same-numbers
std::mt19937 gen(static_cast<unsigned int>(time(nullptr)));



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
double sd(const std::vector<double>& vec){
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0) / vec.size() - mean(vec) * mean(vec));
}
double sd(const std::vector<int>& vec){
    return std::sqrt(std::inner_product(vec.begin(), vec.end(), vec.begin(), 0.0) / vec.size() - mean(vec) * mean(vec));
}
double sgn(const double& val){
    if(val>0.0){
        return 1.0;
    }else if(val<0.0){
        return -1.0;
    }else{
        return 0.0;
    }
}
template<typename T> static inline double Lerp(T v0, T v1, T t){
    return (1 - t)*v0 + t*v1;
}
template<typename T> static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs){
    if (inData.empty()){
        return std::vector<T>();
    }

    if (1 == inData.size()){
        return std::vector<T>(1, inData[0]);
    }

    std::vector<T> data = inData;
    std::sort(data.begin(), data.end());
    std::vector<T> quantiles;

    for (size_t i = 0; i < probs.size(); ++i){
        T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);

        size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
        size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));

        T datLeft = data.at(left);
        T datRight = data.at(right);

        T quantile = Lerp<T>(datLeft, datRight, poi - left);

        quantiles.push_back(quantile);
    }
    return quantiles;
}

//##### RANDOM SAMPLING

double randomDouble(){
//     return ((double) rand() / (RAND_MAX));           // THIS LEADS TO A SLIGHT COMPUTING NUMERICAL ERROR; THE REASON IS UNKNOWN.... INCOMPATIBLE WITH RANDOM NUMBER GENERATOR?
    std::uniform_real_distribution<> dis(0, 1);
    return dis(gen);
}
double randomDouble(double minValue, double maxValue){
    std::uniform_real_distribution<> dis(minValue, maxValue);
    return dis(gen);
}
double randomDoubleExpDist(double meanValue){
    std::exponential_distribution<> expdist(1/meanValue);
    return expdist(gen);
}
int nrand(int min, int max){   // [min, max]
    std::uniform_int_distribution<int> uni(min,max);
    return uni(gen);
}

std::vector<int> sampleFromUniformDiscreteDistributionWithReplacement(const int& sampleSize, const int& maxInt){
    std::vector<int> samples;
    samples.reserve(sampleSize);
    std::uniform_int_distribution<> distribution(0, maxInt);
    for (auto iter = 0; iter < sampleSize; iter++) {
        samples.push_back(distribution(gen));
    }
    return samples;
}


std::vector<int> sampleFromDiscreteDistributionWithReplacement(const int& sampleSize, std::vector<double>& weights){

    // Here are old versions of this function; I keep them here just in case there is a way to optimize them
    // Isaac, you can give a try. This is a function that is very time-consuming; so if we manage to save time on this
    // one that would be great!

        // Version 1

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

        // Version 2

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
        // Version 3
//     std::vector<int> samples;
//     samples.reserve(sampleSize);
//     std::discrete_distribution<int> distribution(weights.begin(), weights.end());
//     for (auto iter = 0; iter < sampleSize; iter++) {
//         samples.push_back(distribution(gen));
//     }


        // Version 4
    std::vector<int> samples;
    samples.reserve(sampleSize);
//     boost::random::discrete_distribution<int> distribution(weights.begin(), weights.end()); // with the boost package
    std::discrete_distribution<int> distribution(weights.begin(), weights.end()); // without the boost package
    for (auto iter = 0; iter < sampleSize; iter++) {
        samples.push_back(distribution(gen));
    }


    return samples;
}

 // does not work but could be useful for later
std::vector<int> sampleFromDiscreteDistribution(const int& sampleSize, std::vector<double> weights){
    // total size vector taken from weights vector
    const size_t totalSize = weights.size();
    std::vector<int> samples;

    if ((int)totalSize <= sampleSize) {   // sample size is high; all numbers are picked

        samples.reserve(totalSize);
        for (int i=0; i<totalSize; ++i) {
            samples.push_back(i);
        }

    }else{ // sample size is small; a subset is sample without replacement

        samples.reserve(sampleSize);

        // from :
        // https://stackoverflow.com/questions/53632441/c-sampling-from-discrete-distribution-without-replacement
        // see also
        // https://stackoverflow.com/questions/48591592/how-do-i-use-boost-random

        // sudo apt-get install libboost-all-dev
        // g++ -std=c++14 -w testOtherFunctions.cpp OtherFunctions.cpp

       // boost::random::uniform_01<> dist;
        std::uniform_real_distribution<> dist(0,1);
        std::vector<double> vals;
        for (auto iter : weights) {
            vals.push_back(std::pow(dist(gen), 1. / iter));
        }
        // Sorting vals, but retain the indices.
        // There is unfortunately no easy way to do this with STL.
        std::vector<std::pair<int, double>> valsWithIndices;
        for (size_t iter = 0; iter < vals.size(); iter++) {
            valsWithIndices.emplace_back(iter, vals[iter]);
        }
        std::sort(valsWithIndices.begin(), valsWithIndices.end(), [](auto x, auto y) {return x.second > y.second; });

        for (auto iter = 0; iter < sampleSize; iter++) {
            samples.push_back(valsWithIndices[iter].first);
        }
    }
    return samples;
}

std::vector<int> sampleFromUniformDistribution0And1(const int& sampleSize){
    std::vector<int> samples;
    samples.reserve(sampleSize);
    for (int i=0; i<sampleSize; ++i) {
        samples.push_back(nrand(0, 1));
    }
    return samples;
}

std::vector<double> sampleFromUniformDistributionBetween0And1(const int& sampleSize){
    std::vector<double> samples;
    samples.reserve(sampleSize);
    for (int i=0; i<sampleSize; ++i) {
        samples.push_back(randomDouble());
    }
    return samples;
}

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

std::vector<double> sampleFromExponentialDistribution(double mean,int size){
    std::exponential_distribution<float> dist(1/mean);
    std::vector<double> samples;
    samples.reserve(size);
    for (int i=0; i<size; ++i) {
        samples.push_back((double) dist(gen));
    }
    return samples;
}

std::vector<int> sampleFromBinomialDiscreteDistribution(double prob, int n, int size){
    std::binomial_distribution<int> distribution(n,prob);
    std::vector<int> samples;
    samples.reserve(size);
    for (int i=0; i<size; ++i) {
        samples.push_back((int) round(distribution(gen)));
    }
    return samples;
}

// MAIN FILE  (COMMENTED TO RUN THE OTHER SCRIPTS)

// #include <chrono>
// using namespace std::chrono;
// int main()
// {
// //     print(randomDouble());
// //     print(nrand(0,2));
//
//     std::vector<double> weights = uniformVectorDouble(500,1.0);
//     std::vector<int> vec;
//     auto start = high_resolution_clock::now();
//     for(int i(0);i<10000;++i){
//         vec=sampleFromDiscreteDistributionWithReplacement(10,weights);
//     }
//     auto stop = high_resolution_clock::now();
//     auto duration = duration_cast<microseconds>(stop - start);
//     std::cout << duration.count() << std::endl;
//
//     return 0;
// }
