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


//####################################################################
// Individual Class
//####################################################################


    ///###### CONSTRUCTOR

Individual::Individual(){}


Individual::Individual( const int& indexPopulation,
                        const int& maxDamageConsidered,
                        const double& probSurvExtrinsicMortality,
                        const double& rateAccumul,
                        const double& convertIntoYear,
                        const std::string& typeAccumulation
                      ){

    _indexPopulation = indexPopulation;
    _maxDamageConsidered = maxDamageConsidered;
    _probSurvExtrinsicMortality = probSurvExtrinsicMortality;
    _rateAccumul = rateAccumul;
    _convertIntoYear = convertIntoYear;
    _typeAccumulation = typeAccumulation;
    _fecundity = 1.0;

    for (int damage(0); damage<_maxDamageConsidered-1; ++damage) {
        _survCausedPerMutation.push_back(1.0);
    }
    _survCausedPerMutation.push_back(0.0);
    for (int damage(0); damage<_maxDamageConsidered-1; ++damage) {
        _nbLethalMutations.push_back(0.0);
    }
    _nbLethalMutations.push_back(1.0);

    reset();    // individual is newborn

}

    ///###### DESTRUCTOR

Individual::~Individual(){};

    ///###### CLASS METHODS

void Individual::reset(){
    _livingState = true;
    _causeDeath = -1;
    _age = 0;
    _damage = 0;
}

void Individual::updateNewTimeStep(){
    if (_livingState==true){

        if (randomDouble()<_probSurvExtrinsicMortality) {
            if ( (_damage<_survCausedPerMutation.size()   && randomDouble()<_survCausedPerMutation[_damage])  ||
                 (_damage>_survCausedPerMutation.size()-1 && randomDouble()<_survCausedPerMutation[_survCausedPerMutation.size()-1])    ) {
                ++_age;
                getLevelDamage();
            }else{
                _livingState = false;
                _causeDeath = 1;
            }
        }else{
            _livingState = false;
            _causeDeath = 0;
        }
    }
}

void Individual::updateNewTimeStep(double& adjustTerm){
    if (_livingState==true){
        if (adjustTerm<1e-20 || randomDouble()<1+(_probSurvExtrinsicMortality-1)*adjustTerm) {
            if ( (_damage<_survCausedPerMutation.size()   && randomDouble()<_survCausedPerMutation[_damage])  ||
                 (_damage>_survCausedPerMutation.size()-1 && randomDouble()<_survCausedPerMutation[_survCausedPerMutation.size()-1])    ) {
                ++_age;
                getLevelDamage();

            }else{
                _livingState = false;
                _causeDeath = 1;
            }
        }else{
            _livingState = false;
            _causeDeath = 0;
        }
    }
}

void Individual::getLevelDamage(){
    if(_typeAccumulation=="lin"){
        _damage = (int) (_age * _rateAccumul);
    }else if(_typeAccumulation=="randomlin"){
        double probAccum = 1-exp(-_rateAccumul/_convertIntoYear);
        if (randomDouble()<probAccum) {
            _damage++;
        }
    }
}

void Individual::shortenMutationVector(bool boolReverseMutation){
    int indexEnd(-1);
    int age(0);
    while (age<_survCausedPerMutation.size() && indexEnd<0) {
        if( (boolReverseMutation==false && _survCausedPerMutation[age]<1e-20 && _nbLethalMutations[age]>0) || (boolReverseMutation==true && _survCausedPerMutation[age]<1e-20 && _nbLethalMutations[age]>10)){
            indexEnd = age;
        }
        ++age;
    }
    if(indexEnd>0 && indexEnd<_survCausedPerMutation.size()-2){
        std::vector<double> newSurvCausedPerMutation = {};
        for (int age(0); age<indexEnd+1; ++age) {
            newSurvCausedPerMutation.push_back(_survCausedPerMutation[age]);
        }
        _survCausedPerMutation = newSurvCausedPerMutation;
        std::vector<double> newNbLethalMutations= {};
        for (int age(0); age<indexEnd+1; ++age) {
            newNbLethalMutations.push_back(_nbLethalMutations[age]);
        }
        _nbLethalMutations = newNbLethalMutations;
    }
}
