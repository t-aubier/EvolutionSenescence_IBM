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
                const double& convertIntoYear,
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
    double _convertIntoYear;                        // factor to convert to year
    std::string _typeAccumulation;                  // type of accumulation of damage: "lin", "exp", "log", or "randomlin"
                                                    // "lin" is implemented in most simulations (it reflects the age), "randomlin" is implemented when we consider stochastic
                                                    // change of somatic state

    bool _livingState;                              // true: living; false: dead
    int _causeDeath;                                // cause of death: -1: not dead; 0: extrinsic mortality; 1: intrinsic mortality

    std::vector<double> _survCausedPerMutation;     // probability of survival per damage level caused by the mutations
    std::vector<double> _nbLethalMutations;         // number of damage-dependent mutations

    int _maxDamageConsidered;                       // max damage considered;

        //constant during lifetime

    int _indexPopulation;                           // index of the individual in the population
    double _fecundity;                              // fecundity; accounts for pleiotropy


    ///###### CLASS METHODS

    // we reset the information of the individual (when replaced by a newborn individual)
    void reset();

    // update of the individual information at each time step; including extrinsic and intrinsic mortality
    void updateNewTimeStep();

    // same as updateNewTimeStep() but this time with adjustTerm that affects the level of extrinsic mortality dependening on population size
    void updateNewTimeStep(double& adjustTerm);

    // update the level of damage (either reflects the chronological age, or increase stochastically)
    void getLevelDamage();

    // shorten the vector _survCausedPerMutation depending on the overall effect of the mutations fixed in the populaiton
    void shortenMutationVector(bool boolReverseMutation);

};


#endif
