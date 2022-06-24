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

    // We initialize parameters and variables
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

    // We initialize the individuals
    _ListLivingInd.reserve(_carryingCapaxity*2);
    _ListDeadInd.reserve(_carryingCapaxity*2);
    _ListAllInd.reserve(_carryingCapaxity*2);
    Individual ind;
    for (int indexPopulation(0); indexPopulation<_carryingCapaxity; ++indexPopulation) {
        ind = Individual( indexPopulation,
                          maxDamageConsidered,
                          probSurvExtrinsicMortality,
                          rateAccumul,
                          convertIntoYear,
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
    // We update individuals at each time step
    if (_densityDependenceSurvival==true && _ListInd.size()<_carryingCapaxity){             // With a density-dependent extrinsic mortality
        // We define adjustTerm that will cause density-dependent extrinsic mortality
        double adjustTerm = ((double)_ListInd.size()) / ((double) _carryingCapaxity);
        for (auto& iter:_ListInd){
            iter.updateNewTimeStep(adjustTerm);
        }
    }else{                                                                                  // Without a density-dependent extrinsic mortality
        for (auto& iter:_ListInd){
            iter.updateNewTimeStep();
        }
    }
}


void Population::replaceDeadIndividuals(){
    // We replace dead individuals by offspring individuals at each time step

    // We update living/dead and list of offspring
    // We add living individuals to the vector _ListIndUpdated
    _ListLivingInd.clear();
    _ListDeadInd.clear();
    _ListIndOfspring.clear();
    _ListIndOfspringIndex.clear();
    _ListIndUpdated.clear();
    for (auto& iter:_ListInd) {
        if (iter._livingState==true) {
            _ListLivingInd.push_back(iter._indexPopulation);
            _ListIndUpdated.push_back(iter);

            // We consider asexual reproduction and production of offspring
            std::poisson_distribution<int> distribution(_nbOffspringPerInd*iter._fecundity);
            int nbOffspring = distribution(gen2);
            for (int nbOff(0); nbOff<nbOffspring; ++nbOff) {
                _ListIndOfspringIndex.push_back(iter._indexPopulation);
            }

        }else{
            _ListDeadInd.push_back(iter._indexPopulation);
        }
    }

    // We consider local competition among offspring
    double sizePop = _ListLivingInd.size() + _ListIndOfspringIndex.size();
    double probSurvOffspring = 0.0;
    if(_carryingCapaxity - (double) _ListLivingInd.size()>0){
        probSurvOffspring = 1.0 / ( 1.0 + ((double) _ListIndOfspringIndex.size() )/((double)  _carryingCapaxity - (double) _ListLivingInd.size() ));
    }
    for (auto& index:_ListIndOfspringIndex) {
        if (randomDouble()<probSurvOffspring){
            _ListIndOfspring.push_back(_ListInd[index]);
        }
    }

    // We add living offspring to the vector _ListIndUpdated
    for (auto& iter:_ListIndOfspring) {
            // We implement mutations
            int age(0);
            bool foundLowestAge = false;
            while (age<iter._survCausedPerMutation.size()) {    // We loop over each damage categorie (i.e., each age categorie assuming non-stochastic linear accumulation)

                if ((_boolReverseMutation==false && iter._nbLethalMutations[age]<1 &&  randomDouble()<_probDeleteriousMutationPerAge) || (_boolReverseMutation==true &&  randomDouble()<_probDeleteriousMutationPerAge)){ // HERE
                    // Deleterious mutations
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

                }else if(_boolReverseMutation==true && iter._nbLethalMutations[age]>0 && randomDouble()<_probDeleteriousMutationPerAge*_ratioReverseMutation){
                    // Beneficial 'reverse' mutations
                    for(int age2(0);age2<_rangeEffectDeleteriousMutation;++age2){
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
                }
                // We consider a change in fecundity in the case of pleiotropic mutations
                if(foundLowestAge == false && iter._survCausedPerMutation[age]<1e-30){
                    iter._fecundity = 1.0 + (_alphaMax - 1.0) * exp(-(_rateAlphaFecundity  * age));
                    foundLowestAge = true;
                }
                age+=_rangeEffectDeleteriousMutation;
            }
            iter.reset();
            _ListIndUpdated.push_back(iter);
    }

    // We update the population
    _ListInd = _ListIndUpdated;
    int index2 = 0;
    for (auto& iter:_ListInd) {
        iter._indexPopulation=index2;
        ++index2;
    }
}

void Population::shortenMutationVector(bool boolReverseMutation){
    // We shorten the vector describing intrinsic mortality depending on the number of damages
    for (auto& iter:_ListInd)
        iter.shortenMutationVector(boolReverseMutation);
}


    ///###### ADDITIONAL METHODS

void Population::setProbDeleteriousMutationPerAge(const double& newProbDeleteriousMutationPerAge){
    // We change the probability of mutations
    _probDeleteriousMutationPerAge = newProbDeleteriousMutationPerAge;
}

void Population::reset(){
    // We reset the population
    for (auto& iter:_ListInd)
        iter.reset();
}
