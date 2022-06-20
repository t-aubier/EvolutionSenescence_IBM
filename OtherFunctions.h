#ifndef OTHER_FUNCTIONS_H
#define OTHER_FUNCTIONS_H

#include <vector>
#include <array>
//####################################################################
// Other functions
//####################################################################



//##### INITIALIZE VECTORS


//##   Create a vector of integers
//     @input  int sizeVec      :  the size of the vector
//     @input  int value        :  the int value of each element of the vector
//     @return std::vector<int> :  a uniform vector with the same int elements
std::vector<int> uniformVectorInt(const int& sizeVec, const int& value);

//##   Create a vector of double
//     @input  int sizeVec      :  the size of the vector
//     @input  double value     :  the double value of each element of the vector
//     @return std::vector<int> :  a uniform vector with the same double elements
std::vector<double> uniformVectorDouble(const int& sizeVec, const double& value);





//##### PRINT TERMINAL


//##   Print variables in the terminal
//##   (useful for debugging)
//     @input  the element to be printed; can be a vector, and array or a number
//     @return print these variables
//     note: if one wants to print another object, one can modify one of these functions
void print(std::vector<int> vec);
void print(std::vector<double> vec);
void print(std::vector<std::array<double,3>> vec);
void print(std::vector<std::vector<double>> vec);
void print(double nb);
void print(int nb);
void print(size_t nb);

//##   Print variables in the terminal with a given string before
//##   (useful for debugging; to annotate the line to be printed)
//     @input  the string to be printed
//     @input  the element to be printed; can be a vector, and array or a number
//     @return print this string and these variables
//     note: if one wants to print another object, one can modify one of these functions
void print(std::string str, std::vector<int> vec);
void print(std::string str, std::vector<double> vec);
void print(std::string str, std::vector<std::array<double,3>> vec);
void print(std::string str, std::vector<std::vector<double>> vec);
void print(std::string str, double nb);
void print(std::string str, int nb);
void print(std::string str, size_t nb);
void print(std::string str);





//##### BASIC STATS


//##   Compute the mean value of double elements stored in a vector
//     @input  std::vector<double> vec  :  the vector containing the double elements
//     @return double                   :  the mean values of those double elements
double mean(const std::vector<double>& vec);

//##   Compute the mean value of int elements stored in a vector
//     @input  std::vector<double> vec  :  the vector containing the int elements
//     @return double                   :  the mean values of those int elements
double mean(const std::vector<int>& vec);

//##   Compute the standard deviation of double elements stored in a vector
//     @input  std::vector<double> vec  :  the vector containing the double elements
//     @return double                   :  the standard deviation of those int elements
double sd(const std::vector<double>& vec);

//##   Compute the standard deviation of int elements stored in a vector
//     @input  std::vector<double> vec  :  the vector containing the int elements
//     @return double                   :  the standard deviation of those int elements
double sd(const std::vector<int>& vec);

//##   Compute the sign of a double value
//     @input  double val :  the double value
//     @return doubl      :  either -1.0 (if val<0), 0.0 (if val==0), or +1.0 (if val>0)
double sgn(const double& val);

//##   Function used in the Quantile() function (see below
template<typename T> static inline double Lerp(T v0, T v1, T t);

//##   Compute the quantile from a vector of numerics
//     @input  std::vector<T>& inData :  the vector containing the numeric elements
//     @input  std::vector<T>& probs  :  the quantiles
//     @return std::vector<T>         :  the numerics corresponding to the quantiles
template<typename T> static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs);



//##### RANDOM SAMPLING

//##   Compute a random double number in the range [0,1]
//     @return double                   :  the random number
double randomDouble();

//##   Compute a random double number following a uniform distribution in a given range
//     @input  double  minValue  :  the minimum value that can be drawn
//     @input  double  maxValue  :  the maximum value that can be drawn
//     @return double            :  the random number belonging to the range [minValue, maxValue]
double randomDouble(double minValue, double maxValue);

//##   Compute a random double number following an exponential distribution with a given mean
//     @input  double  meanValue :  the mean value describing the exponential distribution
//                                  (=1/lambda with lambda being the parameter that usually describes such distribution)
//     @return double            :  the random number drawn in this exponential distribution
double randomDoubleExpDist(double meanValue);

//##   Compute a random int number in a given range
//     @input  int  min  :  the minimum int value that can be drawn
//     @input  int  max  :  the maximum int value that can be drawn
//     @return int       :  the random number belonging to the range {minValue, minValue+1, ... maxValue}
int nrand(int min, int max);


//##   Compute a list of int number from 0 to maxInt following a uniform distribution; with replacement
//     @input  int  sampleSize  :   the number of int numbers to sample
//     @input  int  maxInt      :  the maximum int value that can be drawn
std::vector<int> sampleFromUniformDiscreteDistributionWithReplacement(const int& sampleSize, const int& maxInt);

//##   Compute a list of int numbers following a probability distribution *with* replacement,
//##   so that each int number can be sampled more than once
//     @input  int  sampleSize              :  the number of int numbers to sample
//     @input  std::vector<double> weights  :  a vector where each double element represent the weight associated with each element
//     @return std::vector<int>             :  a vector of size 'sampleSize' with each element drawn following its weight
//                                             for instance, if weight = {1.0,0.1,0.1} and sampleSize = 10,
//                                             each of the 10 elements of the output vector will be 0 with a probability  1  / 1.2
//                                                                                                  1 with a probability 0.1 / 1.2
//                                                                                                  2 with a probability 0.1 / 1.2
std::vector<int> sampleFromDiscreteDistributionWithReplacement(const int& sampleSize, std::vector<double>& weights);

//##   Compute a list of int numbers 0 and 1 following a uniform distribution,
//##   so that each int number can be sampled more than once
//     @input  int  sampleSize              :  the number of int numbers to sample
//     @return std::vector<int>             :  a vector of size 'sampleSize' with each element being 0 or 1, following a uniform distribution
std::vector<int> sampleFromUniformDistribution0And1(const int& sampleSize);

//##   Compute a list of double numbers between 0 and 1 following a uniform distribution,
//##   so that each int number can be sampled more than once
//     @input  int  sampleSize              :  the number of int numbers to sample
//     @return std::vector<double>          :  a vector of size 'sampleSize' with each element being between 0 or 1, following a uniform distribution
std::vector<double> sampleFromUniformDistributionBetween0And1(const int& sampleSize);

//##   Compute a list of int numbers following a probability distribution *without* replacement,
//##   so that each int number can only be sampled once
//     @input  int  sampleSize              :  the number of int numbers to sample; this number cannot be higher than the size of input weights
//     @input  std::vector<double> weights  :  a vector where each double element represent the weight associated with each element
//     @return std::vector<int>             :  a vector of size 'sampleSize' with each element drawn following its weight
//                                             for instance, if weight = {1.0,0.1,0.1} and sampleSize = 3,
//                                             the output will comprise all elements; e.g. {0,2,1}, or {2,1,0} because each element can
//                                             only be sampled once
std::vector<int> sampleFromDiscreteDistribution(const int& sampleSize, std::vector<double> weights);
// does not work but could be useful for later

//##   Compute a list of int numbers following a normal distribution, *with* replacement
//     @input  double  mean                 :  the mean characterizing the normal distribution
//     @input  double  sd                   :  the standard deviation characterizing the normal distribution
//     @input  int  size                    :  the number of int numbers to sample
//     @return std::vector<int>             :  a vector of size 'size' with each int element drawn following the normal distribution
std::vector<int> sampleFromNormalDiscreteDistribution(double mean, double sd, int size);

//##   Compute a list of double numbers following a normal distribution, *with* replacement
//     @input  double  mean                 :  the mean characterizing the normal distribution
//     @input  double  sd                   :  the standard deviation characterizing the normal distribution
//     @input  int  size                    :  the number of int numbers to sample
//     @return std::vector<double>          :  a vector of size 'size' with each double element drawn following the normal distribution
std::vector<double> sampleFromNormalDistribution(double mean, double sd, int size);

//##   Compute a list of double numbers following an exponential distribution, *with* replacement
//     @input  double  mean                 :  the mean characterizing the exponential distribution
//                                             (=1/lambda with lambda being the parameter that usually describes such distribution)
//     @input  int  size                    :  the number of int numbers to sample
//     @return std::vector<double>             :  a vector of size 'size' with each double element drawn following the exponential distribution
std::vector<double> sampleFromExponentialDistribution(double mean,int size);

//##   Compute a list of int numbers following a normal distribution, *with* replacement
//     @input  double  prob                 :  parameter p describing the binomial distribution
//     @input  int  n                       :  parameter N describing the binomial distribution
//     @input  int  size                    :  the number of int numbers to sample
//     @return std::vector<int>             :  a vector of size 'size' with each element drawn following the binomial distribution
std::vector<int> sampleFromBinomialDiscreteDistribution(double prob, int n, int size);



#endif
