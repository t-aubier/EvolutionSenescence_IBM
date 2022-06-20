#ifndef CLASS_POPULATION_H
#define CLASS_POPULATION_H

#include "OtherFunctions.h"
#include "ClassIndividual.h"

#include <vector>
#include <array>


//####################################################################
// Class to describe populations
//####################################################################

class Population
{

public:

    ///###### CONSTRUCTORS

    Population();

    Population( const int& carryingCapaxity,
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
              );

    ///###### DESTRUCTOR

    ~Population();

    ///###### CLASS MEMBERS

    int _carryingCapaxity;                              // Carrying capacity
    bool _densityDependenceSurvival;                    // Whether survival is adjusted when N<CarryingCapacity after recruitment
    double _nbOffspringPerInd;                          // Nb of offspring per individuals
    double _alphaMax;                                   // Maximum pleiotropic effect (when mutation expressed at damage == 0)
    double _rateAlphaFecundity;                         // rate of pleiotropy
    int _maxDamageConsidered;                           // max damage considered;

    double _probDeleteriousMutationPerAge;              // probability of getting the deleterious mutation
    double _ratioReverseMutation;                       // proportion of beneficial reverse mutations relative to deleterious mutations

    double _convertIntoYear;                               // factor to convert to year
    int _rangeEffectDeleteriousMutation;                // discretization of the effect of the deleterious mutation
    bool _boolReverseMutation;                          // whether lethal mutation can be removed or not
    double _effectSurvDeleteriousMutation;              // effect of each deleterious mutation on survival


    std::vector<Individual> _ListInd;                   // Population of individuals
    std::vector<int> _ListLivingInd;                    // list of index of parental Ind
    std::vector<int> _ListDeadInd;                      // list of index of dead Ind
    std::vector<int> _ListAllInd;                       // list of all individuals

    std::vector<Individual> _ListIndUpdated;            // Population of individuals updated after replacement
    std::vector<Individual> _ListIndOfspring;           // Population of offspring individuals
    std::vector<int> _ListIndOfspringIndex;             // List of indexes of the parents for each offspring


    ///###### CLASS METHODS

    void updateNewTimeStep();
    void replaceDeadIndividuals();
    void shortenMutationVector(bool boolReverseMutation);

    ///###### ADDITIONAL METHODS

    void setProbDeleteriousMutationPerAge(const double& newProbDeleteriousMutationPerAge);
    void reset();

};


#endif
