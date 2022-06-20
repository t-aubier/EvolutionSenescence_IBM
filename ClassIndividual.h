#ifndef CLASS_INDIVIDUAL_H
#define CLASS_INDIVIDUAL_H

#include "OtherFunctions.h"

#include <vector>
#include <array>


//####################################################################
// Class to describe individuals
//####################################################################

class Individual
{

public:

    ///###### CONSTRUCTORS

    Individual();

    Individual( const int& indexPopulation,
                const int& maxDamageConsidered,
                const double& probSurvExtrinsicMortality,
                const double& rateAccumul,
                const std::string& typeAccumulation
              );

    ///###### DESTRUCTOR

    ~Individual();

    ///###### CLASS MEMBERS

    // True for all individuals

    double _probSurvExtrinsicMortality;              // probability of survival considering only extrinsic mortality


    // Specific to each individual

        //change during lifetime

    int _age;                                       // age
    int _damage;                                    // damage
    double _rateAccumul;                            // rate of damage accumulation
    std::string _typeAccumulation;                  // type of accumulation: "lin", "exp" or "log"

    bool _livingState;                              // true: living; false: dead
    int _causeDeath;                                // -1: not dead; 0: extrinsic mortality; 1: intrinsic mortality

    std::vector<double> _survCausedPerMutation;     // probability of survival per damage level caused by the mutations
    std::vector<double> _nbLethalMutations;         // number of damage-dependent mutations

    int _maxDamageConsidered;                       // max damage considered;

        //constant during lifetime

    int _indexPopulation;                           // index of the individual in the population
    double _fecundity;                              // fecundity; accounts for pleiotropy

    
    ///###### CLASS METHODS

    void reset();
    void updateNewTimeStep();
    void updateNewTimeStep(double& adjustTerm);
    void getLevelDamage();
    void shortenMutationVector(bool boolReverseMutation);


};


#endif
