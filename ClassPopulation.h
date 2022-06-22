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
    double _alphaMax;                                   // Maximum fecundity when pleiotropy (when mutation expressed at age == 0)
    double _rateAlphaFecundity;                         // Rate at which fecundity decreases with the age at which intrinsic mortality occurs (called gamma in the manuscript)
    int _maxDamageConsidered;                           // Maximum damage considered (correspond to maximum age in years);
                                                        // to determine the size of vectors (should be high enough)

    double _probDeleteriousMutationPerAge;              // Probability of getting a deleterious mutation expressed at a given age

    double _convertIntoYear;                            // Depending on the simulation: if time steps = months: 12 ; if time steps = days: 365
                                                        // Considering that time steps are months instead of days decreases simulation runtime

    int _rangeEffectDeleteriousMutation;                // Age at which deleterious mutations can have an effect (if=1: in months if convertIntoYear=12, and in days if convertIntoYear=365)

    bool _boolReverseMutation;                          // If true, beneficial 'reverse' mutations can occur
    double _ratioReverseMutation;                       // If boolReverseMutation==true, proportion of beneficial mutations relative to deleterious mutations

    double _effectSurvDeleteriousMutation;              // Survival probability reduced due to the expression of a single deleterious mutation (if ==1, lethal mutation)


    std::vector<Individual> _ListInd;                   // Population of individuals
    std::vector<int> _ListLivingInd;                    // List of index of parental individuals
    std::vector<int> _ListDeadInd;                      // List of index of dead individuals
    std::vector<int> _ListAllInd;                       // List of all individuals

    std::vector<Individual> _ListIndUpdated;            // Population of individuals updated after replacement
    std::vector<Individual> _ListIndOfspring;           // Population of offspring individuals
    std::vector<int> _ListIndOfspringIndex;             // List of indexes of the parents for each offspring


    ///###### CLASS METHODS

    // update of the individual information at each time step; including extrinsic and intrinsic mortality;
    // with or without density-dependent extrinsic mortality
    void updateNewTimeStep();

    // replace dead individuals by newborn individuals; include local competition among offspring, deleterious and beneficial mutations,
    // change in fecundity due to pleiotropy
    void replaceDeadIndividuals();

    // reduce the size of the vector storing mutations, when the overall effects of mutation makes that no individual can reach old ages
    void shortenMutationVector(bool boolReverseMutation);

    ///###### ADDITIONAL METHODS

    // change the probability of deleterious mutations
    void setProbDeleteriousMutationPerAge(const double& newProbDeleteriousMutationPerAge);

    // reset the information of the individual (when replaced by a newborn individual)
    void reset();

};


#endif
