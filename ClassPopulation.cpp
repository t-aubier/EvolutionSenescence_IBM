#include "ClassPopulation.h"
#include "ClassIndividual.h"
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

#include <boost/random/uniform_01.hpp>
#include <boost/random/mersenne_twister.hpp>

std::random_device rd2;
boost::random::mt19937 gen2(rd2());


//####################################################################
// Population Class
//####################################################################


    ///###### CONSTRUCTOR

Population::Population(){}


Population::Population( const int& carryingCapaxity,
                        const bool& densityDependenceSurvival,
                        const double& nbOffspringPerInd,
                        const double& alphaMax,
                        const double& rateAlphaFecundity,
                        const int& maxDamageConsidered,
                        const double& probSurvExtrinsicMortality,
                        const double& rateAccumul,
                        const std::string& typeAccumulation,
                        const double& probDeleteriousMutationPerAge,
                        const double& ratioReverseMutation,
                        const double& convertIntoYear,
                        const int& rangeEffectDeleteriousMutation,
                        const bool& boolReverseMutation,
                        const double& effectSurvDeleteriousMutation
                      ){
    
    _carryingCapaxity = carryingCapaxity;
    _densityDependenceSurvival = densityDependenceSurvival;
    _nbOffspringPerInd = nbOffspringPerInd;
    _alphaMax = alphaMax;
    _rateAlphaFecundity = rateAlphaFecundity;

    _maxDamageConsidered = maxDamageConsidered;

    _probDeleteriousMutationPerAge = probDeleteriousMutationPerAge;
    _ratioReverseMutation = ratioReverseMutation;
    _convertIntoYear = convertIntoYear;
    _rangeEffectDeleteriousMutation = rangeEffectDeleteriousMutation;
    _boolReverseMutation = boolReverseMutation;
    _effectSurvDeleteriousMutation = effectSurvDeleteriousMutation;
        
    _ListLivingInd.reserve(_carryingCapaxity*2);
    _ListDeadInd.reserve(_carryingCapaxity*2);
    _ListAllInd.reserve(_carryingCapaxity*2);

    Individual ind;
    for (int indexPopulation(0); indexPopulation<_carryingCapaxity; ++indexPopulation) {
        ind = Individual( indexPopulation,
                          maxDamageConsidered,
                          probSurvExtrinsicMortality,
                          rateAccumul,
                          typeAccumulation
                        );
        
        _ListInd.push_back(ind);
        _ListLivingInd.push_back(indexPopulation);
        _ListAllInd.push_back(indexPopulation);

    }
}



    ///###### DESTRUCTOR

Population::~Population(){};



    ///###### CLASS METHODS

void Population::updateNewTimeStep(){
    
    // update living state, damage accumulation, and survival probability
    
    if (_densityDependenceSurvival==true && _ListInd.size()<_carryingCapaxity){
        double adjustTerm = ((double)_ListInd.size()) / ((double) _carryingCapaxity);
        if(adjustTerm<0.95){
//             print((double)_ListInd.size());
            adjustTerm = 0.0;
        }
        for (auto& iter:_ListInd){
            iter.updateNewTimeStep(adjustTerm);
        }    
    }else{
        for (auto& iter:_ListInd){
            iter.updateNewTimeStep();
        }
    }

}

    
void Population::replaceDeadIndividuals(){
    
    // update living/dead and list of offspring
    _ListLivingInd.clear();
    _ListDeadInd.clear();
    _ListIndOfspring.clear();
    _ListIndOfspringIndex.clear();
    _ListIndUpdated.clear();
    for (auto& iter:_ListInd) {
        if (iter._livingState==true) {
            _ListLivingInd.push_back(iter._indexPopulation);
            _ListIndUpdated.push_back(iter);

            std::poisson_distribution<int> distribution(_nbOffspringPerInd*iter._fecundity);
            int nbOffspring = distribution(gen2);

            for (int nbOff(0); nbOff<nbOffspring; ++nbOff) {
                //_ListIndOfspring.push_back(iter);
                _ListIndOfspringIndex.push_back(iter._indexPopulation);  // HERE
            }
            
        }else{
            _ListDeadInd.push_back(iter._indexPopulation);
        }
    }
    
        // local competition among offspring; and mutation accumulation
    double sizePop = _ListLivingInd.size() + _ListIndOfspringIndex.size();
//     if(sizePop==0) {
//         std::cout << "EXTINCTION" << std::endl;   
//     }
    double probSurvOffspring = 0.0;
    if(_carryingCapaxity - (double) _ListLivingInd.size()>0){
        probSurvOffspring = 1.0 / ( 1.0 + ((double) _ListIndOfspringIndex.size() )/((double)  _carryingCapaxity - (double) _ListLivingInd.size() ));        
    }

    
//     for (auto& index:_ListIndOfspringIndex) {
//         if (randomDouble()<probSurvOffspring){
//             
//             for (auto& iter:_ListInd) { 
//                 if (iter._indexPopulation==index) {
//                     _ListIndOfspring.push_back(iter);
//                 }
//             }
//         }
//     }
    for (auto& index:_ListIndOfspringIndex) {
        if (randomDouble()<probSurvOffspring){
            _ListIndOfspring.push_back(_ListInd[index]);
            
//             print("");
//             print(index);
//             print(_ListIndOfspring[_ListIndOfspring.size()-1]._indexPopulation);
        }
    }
    
    
    
//     print("here");
//     print(_ListIndOfspringIndex.size());
//     print(_ListIndOfspring.size());
    for (auto& iter:_ListIndOfspring) {   //HERE
        //if (randomDouble()<probSurvOffspring){       // offspring survives  //HERE
            //mutation
//             int indexEnd(-1);
            int age(0);
//             while (age<iter._survCausedPerMutation.size() && indexEnd<0) {
            bool foundLowestAge = false;
            while (age<iter._survCausedPerMutation.size()) {
//                 if (iter._survCausedPerMutation[age]>1e-30 && randomDouble()<_probDeleteriousMutationPerAge){
//                 if (randomDouble()<_probDeleteriousMutationPerAge){
                
                if ((_boolReverseMutation==false && iter._nbLethalMutations[age]<1 &&  randomDouble()<_probDeleteriousMutationPerAge) || (_boolReverseMutation==true &&  randomDouble()<_probDeleteriousMutationPerAge)){ // HERE
                    
                    for(int age2(0);age2<_rangeEffectDeleteriousMutation;++age2){
                        if(age+age2<iter._survCausedPerMutation.size()){
                            iter._survCausedPerMutation[age+age2] -= _effectSurvDeleteriousMutation;
                            if(_boolReverseMutation==true){
                                iter._nbLethalMutations[age+age2] += 1.0;
                            }
                            if(iter._survCausedPerMutation[age+age2]<1e-30){
                                iter._survCausedPerMutation[age+age2] = 0.0;
                            }                                      
                        }              
                    }
                    
//                 }else if(_boolReverseMutation==true && iter._survCausedPerMutation[age]<1- 1e-30 && randomDouble()<_probDeleteriousMutationPerAge){
                }else if(_boolReverseMutation==true && iter._nbLethalMutations[age]>0 && randomDouble()<_probDeleteriousMutationPerAge*_ratioReverseMutation){

                    for(int age2(0);age2<_rangeEffectDeleteriousMutation;++age2){
//                         if(iter._nbLethalMutations[age+age2]>2.5){
//                             print(iter._nbLethalMutations[age+age2]);
//                             print(age+age2);
//                         }
                        
                        
                        
                        
                        if(age+age2<iter._survCausedPerMutation.size()-1){
                            iter._nbLethalMutations[age+age2] -= 1.0;
                            if(iter._nbLethalMutations[age+age2]<0.5){
                                iter._survCausedPerMutation[age+age2] += _effectSurvDeleteriousMutation;
                                if(iter._survCausedPerMutation[age+age2]>1- 1e-30){
                                    iter._survCausedPerMutation[age+age2] = 1.0;
                                }
                            }  
                        }
      
                        
                        
                        
                        
                        
                    }
                    
//                     print(iter._nbLethalMutations);
//                     print(iter._survCausedPerMutation);
                }
//                 if(iter._survCausedPerMutation[age]<1){
//                     print(age);
//                     print(iter._survCausedPerMutation[age]);
//                 }

                if(foundLowestAge == false && iter._survCausedPerMutation[age]<1e-30){
                    iter._fecundity = 1.0 + (_alphaMax - 1.0) * exp(-(_rateAlphaFecundity  * age));
                    foundLowestAge = true;
                }
                age+=_rangeEffectDeleteriousMutation;
                
                
                
                
//                 if (iter._survCausedPerMutation[age]>1e-30 && randomDouble()<_probDeleteriousMutationPerAge){
//                     if(age<iter._survCausedPerMutation.size()){
//                         iter._survCausedPerMutation[age] -= _effectSurvDeleteriousMutation;
//                         if(iter._survCausedPerMutation[age]<1e-30){
//                             iter._survCausedPerMutation[age] = 0.0;
//                         }                                      
//                     }              
//                 }
//                 ++age;
                
            }
//             if(indexEnd>0){
//                 std::vector<double> newSurvCausedPerMutation = {};
//                 for (int age(0); age<indexEnd+1; ++age) {
//                     newSurvCausedPerMutation.push_back(iter._survCausedPerMutation[age]);
//                 }
//                 iter._survCausedPerMutation = newSurvCausedPerMutation;
//             }

            iter.reset();
            _ListIndUpdated.push_back(iter);
        //} // HERE
    }
    
    
    // update population
    _ListInd = _ListIndUpdated;
    int index2 = 0;
    for (auto& iter:_ListInd) { 
        iter._indexPopulation=index2;
        ++index2;
    }
    
}


void Population::shortenMutationVector(bool boolReverseMutation){
    for (auto& iter:_ListInd)
        iter.shortenMutationVector(boolReverseMutation);
}


    ///###### ADDITIONAL METHODS
    
void Population::setProbDeleteriousMutationPerAge(const double& newProbDeleteriousMutationPerAge){
    _probDeleteriousMutationPerAge = newProbDeleteriousMutationPerAge;
}

void Population::reset(){
    for (auto& iter:_ListInd)
        iter.reset();
}
