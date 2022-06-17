#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H

#include <vector>
#include <array>

//####################################################################
// Other functions
//####################################################################



//##### INITIALIZE VECTORS

std::vector<int> uniformVectorInt(const int& sizeVec, const int& value);
std::vector<double> uniformVectorDouble(const int& sizeVec, const double& value);

//##### PRINT TERMINAL

void print(std::vector<int> vec);
void print(std::vector<double> vec);
void print(std::vector<std::array<double,3>> vec);
void print(std::vector<std::vector<double>> vec);
void print(double nb);
void print(int nb);
void print(size_t nb);

void print(std::string str, std::vector<int> vec);
void print(std::string str, std::vector<double> vec);
void print(std::string str, std::vector<std::array<double,3>> vec);
void print(std::string str, std::vector<std::vector<double>> vec);
void print(std::string str, double nb);
void print(std::string str, int nb);
void print(std::string str, size_t nb);
void print(std::string str);

//##### BASIC STATS

double mean(const std::vector<double>& vec);
double mean(const std::vector<int>& vec);
double sd(const std::vector<int>& vec);
double sd(const std::vector<double>& vec);

//##### RANDOM SAMPLING

double randomDouble();
int nrand(int min, int max);

std::vector<int> sampleFromDiscreteDistributionWithReplacement(const int& sampleSize, std::vector<double> weights); // with replacement
// std::vector<int> sampleFromDiscreteDistribution(const int& sampleSize, std::vector<double> weights); // without replacement

std::vector<int> sampleFromNormalDiscreteDistribution(double mean, double sd, int size);
std::vector<double> sampleFromNormalDistribution(double mean, double sd, int size);

//##### RUNTIME FUNCTIONS

void testRunTime();

#endif
