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


    ///###### FUNCTIONS TO FIND QUANTILES

    //##   Compute the quantile from a vector of numerics
    //     @input  std::vector<T>& inData :  the vector containing the numeric elements
    //     @input  std::vector<T>& probs  :  the quantiles
    //     @return std::vector<T>         :  the numerics corresponding to the quantiles

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

    ///###### SENSITIVITY ANALYSIS

void SensitivityAnalysis        (   int carryingCapacity,
                                    bool densityDependenceSurvival,
                                    int maxDamageConsidered,
                                    double convertIntoYear,
                                    bool startAtEarlierLifeSpan,
                                    double quantileStart,
                                    double probDeleteriousMutationPerAge,
                                    std::string typeAccumulation,
                                    double rateAccumul,
                                    double ratioReverseMutation,
                                    int rangeEffectDeleteriousMutation,
                                    bool boolReverseMutation,
                                    double effectSurvDeleteriousMutation,
                                    int Tburnin,
                                    int Tmax,
                                    std::vector<double> vecBirthRatePerYear,
                                    double alphaMax,
                                    double rateAlphaFecundity,
                                    std::vector<double> vecMortRatePerYear,
                                    std::string Namefile
                           ){

    // We initialize the files where information will be stored

    std::fstream dataSensitivity;
    dataSensitivity.open ("Data/dataSensitivity_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataSensitivity << "BirthRatePerYear;" << "MortRatePerYear;" << "RateAccumul;" << "Time;" << "SizePop;"  << "ExtStatus;"  << "PercentOffspringDead;" << "propDeath;" << "propDeathExtrinsic;" << "propDeathMutation;" << "LifeSpan0025;" << "LifeSpan05;" << "LifeSpan0975;" << "LifeSpanMin;" << "LifeSpanMax;" << "Damage0025;" << "Damage05;" << "Damage0975;" << "DamageMin;" << "DamageMax"  << std::endl;

    for (auto& mortRatePerYear:vecMortRatePerYear) {        // We loop over the different mortality rates tested

        // We convert the rate into a probability per time step
        double probSurvExtrinsicMortality = exp(-mortRatePerYear/convertIntoYear);

        for (auto& birthRatePerYear:vecBirthRatePerYear) {  // We loop over the different birth rates tested

            // We convert the rate into a probability per time step
            double convertIntoYearDAY = 365;
            double probSurvExtrinsicMortalityDAY = exp(-mortRatePerYear/convertIntoYearDAY);
            double nbOffspringPerIndDAY = 1-exp(-birthRatePerYear/convertIntoYearDAY);
            double nbOffspringPerInd = nbOffspringPerIndDAY * convertIntoYearDAY/convertIntoYear  / probSurvExtrinsicMortality;

            // We initialize the population
            Population population = Population( carryingCapacity,
                                                densityDependenceSurvival,
                                                nbOffspringPerInd,
                                                1.0,                           // alphaMax == 1
                                                rateAlphaFecundity,
                                                maxDamageConsidered,
                                                probSurvExtrinsicMortality,
                                                rateAccumul,
                                                typeAccumulation,
                                                0.0,                            // prob of getting a deleterious mutation = 0.0 here
                                                ratioReverseMutation,
                                                convertIntoYear,
                                                rangeEffectDeleteriousMutation,
                                                boolReverseMutation,
                                                effectSurvDeleteriousMutation
                                            );

            // We run the burn-in phase
            bool ExtinctionPop = false;
            int t(0);
            while(t<Tburnin*convertIntoYear+1 && ExtinctionPop == false) {
                if(population._ListInd.size()==0){
                    ExtinctionPop = true;
                }else{
                    population.updateNewTimeStep();             // step 1 at each time step
                    population.replaceDeadIndividuals();        // step 2 at each time step
                }
                ++t;
            }

            if(ExtinctionPop==true){            // We record information if the population gets extinct even without mutation accumulation

                dataSensitivity << birthRatePerYear << ";"
                                << mortRatePerYear << ";"
                                << "NA" << ";"
                                << t << ";"
                                << "NA" << ";"
                                << "ini ext" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << ";"
                                << "NA" << std::endl;
            }else{                              // Else we record information at t = 0, after the burn-in phase
                t = 0;

                Population populationSave = population;
                std::vector<int> ageDeathVec = {};
                std::vector<int> damageDeathVec = {};
                double popSize(0.0);
                double propOffspringDead(0.0);
                double count(0.0);
                double count2(0.0);
                double count3(0.0);
                double count4(0.0);

                double propDeathExtrinsic(0.0);
                double propDeathMutation(0.0);
                double propDeath(0.0);
                for (int rep2(0); rep2<10; ++rep2){         // we 10 simulations lasting 1,000 generation to calculate mean values
                    population = populationSave;            // we consider the population after the burn-in phase
                    for (int t(0);t<1e3;++t) {

                        population.updateNewTimeStep();     // step 1 at each time step


                        // We measure the number of individuals, of dead individuals
                        // (those dying of extrinsic mortality and those dying from intrinsic mortality)
                        for (auto iter:population._ListInd) {
                            count4+=1.0;
                            if(iter._livingState==false){
                                propDeath+=1.0;
                                count3+=1.0;
                                ageDeathVec.push_back(iter._age);
                                damageDeathVec.push_back(iter._damage);
                                if(iter._causeDeath==0){            // death extrinsic mortality
                                    propDeathExtrinsic+=1.0;
                                }else if(iter._causeDeath==1){      // death mutation
                                    propDeathMutation+=1.0;
                                }
                            }
                        }

                        population.replaceDeadIndividuals();        // step 2 at each time step

                        // We measure population size (we consider the mean value later)
                        popSize += (double) population._ListInd.size();


                        // We measure the proportion of offspring that die before reaching sexual maturity
                        if(population._ListIndOfspringIndex.size()>0){
                            propOffspringDead += ((double) population._ListIndOfspringIndex.size() - ((double) population._ListInd.size() - population._ListLivingInd.size())) / ((double) population._ListIndOfspringIndex.size()) ;
                            count2+=1.0;
                        }
                        count+=1.0;
                    }

                }
                population = populationSave;                        // We consider the population after the burn-in phase

                // we measure the different quantiles (age of death, and level of damage, i.e., somatic state, at death)
                std::vector<double> ageDeathVecDouble(ageDeathVec.begin(), ageDeathVec.end());
                auto quantiles = Quantile<double>(ageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });
                std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });

                // In simulations were startAtEarlierLifeSpan==true, we implement a first lethal mutation and we thus re-initialize the population
                if(startAtEarlierLifeSpan==true){
                    // Age at which the first lethal mutation is expressed (defined by quantileStart)
                    int newMaxDamageConsidered = (int) (round(quantilesDamage[5])+1);

                    // population with limit in maxAgeConsidered
                    population = Population( carryingCapacity,
                                            densityDependenceSurvival,
                                            nbOffspringPerInd,
                                            alphaMax,
                                            rateAlphaFecundity,
                                            newMaxDamageConsidered,
                                            probSurvExtrinsicMortality,
                                            rateAccumul,
                                            typeAccumulation,
                                            1.0,                            // prob of getting a deleterious mutation = 0.0 here
                                            ratioReverseMutation,
                                             convertIntoYear,
                                            rangeEffectDeleteriousMutation,
                                            boolReverseMutation,
                                            effectSurvDeleteriousMutation
                                        );
                    // burn in phase
                    int t2(0);
                    while(t2<Tburnin*convertIntoYear+1 && ExtinctionPop == false) {
                        if(population._ListInd.size()==0){
                            ExtinctionPop = true;
                        }else{
                            population.updateNewTimeStep();             // step 1 at each time step
                            population.replaceDeadIndividuals();        // step 2 at each time step
                        }
                        ++t2;
                    }
                    if(ExtinctionPop==true){
                        ExtinctionPop = false;
                        population = populationSave;
                    }
                }

                // We consider a non-zero probability to have a deleterious mutation
                population.setProbDeleteriousMutationPerAge(probDeleteriousMutationPerAge);
            }

            // We run the actual simulation until all generations are considered or until population extinction
            bool popGettingExtinct = false;
            ++t;
            while(t<Tmax*convertIntoYear+1 && ExtinctionPop == false) {
                if(population._ListInd.size()==0){
                    ExtinctionPop = true;
                }else{
                    // To speed-up the simulation, we consider less and less age at which mutations can be expressed (once no individual reach this age)
                    if(t% ((int) 1e4)==0){
                        population.shortenMutationVector(boolReverseMutation);
                    }

                    // We record information when the population gets clost to extinction
                    if(((double) population._ListInd.size())/((double) carryingCapacity)<0.1 && popGettingExtinct==false){
                        popGettingExtinct=true;
                        ExtinctionPop = true;

                        // We save data at the population level
                        Population populationSave = population;
                        std::vector<int> ageDeathVec = {};
                        std::vector<int> damageDeathVec = {};
                        double popSize(0.0);
                        double propOffspringDead(0.0);
                        double count(0.0);
                        double count2(0.0);
                        double count3(0.0);
                        double count4(0.0);

                        double propDeathExtrinsic(0.0);
                        double propDeathMutation(0.0);
                        double propDeath(0.0);

                        for (int rep2(0); rep2<10; ++rep2){         // we 10 simulations lasting 1,000 generation to calculate mean values
                            population = populationSave;            // we consider each time the same population
                            for (int t(0);t<1e3;++t) {

                                population.updateNewTimeStep();     // step 1 at each time step


                                // We measure the number of individuals, of dead individuals
                                // (those dying of extrinsic mortality and those dying from intrinsic mortality)
                                for (auto iter:population._ListInd) {
                                    count4+=1.0;
                                    if(iter._livingState==false){
                                        propDeath+=1.0;
                                        count3+=1.0;
                                        ageDeathVec.push_back(iter._age);
                                        damageDeathVec.push_back(iter._damage);
                                        if(iter._causeDeath==0){
                                            propDeathExtrinsic+=1.0;
                                        }else if(iter._causeDeath==1){
                                            propDeathMutation+=1.0;
                                        }
                                    }
                                }

                                population.replaceDeadIndividuals();    // step 2 at each time step

                                // We measure population size (we consider the mean value later)
                                popSize += (double) population._ListInd.size();

                                // We measure the proportion of offspring that die before reaching sexual maturity
                                if(population._ListIndOfspringIndex.size()>0){
                                    propOffspringDead += ((double) population._ListIndOfspringIndex.size() - ((double) population._ListInd.size() - population._ListLivingInd.size())) / ((double) population._ListIndOfspringIndex.size()) ;
                                    count2+=1.0;
                                }

                                count+=1.0;
                            }

                        }
                        population = populationSave;    // We consider the population as it was before the analysis

                        // we measure the different quantiles (age of death, and level of damage, i.e., somatic state, at death)
                        std::vector<double> ageDeathVecDouble(ageDeathVec.begin(), ageDeathVec.end());
                        auto quantiles = Quantile<double>(ageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1  });
                        std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                        auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1  });

                        // We record the information at the population level
                        dataSensitivity << birthRatePerYear << ";"
                                        << mortRatePerYear << ";"
                                        << rateAccumul << ";"
                                        << t << ";"
                                        << popSize/count << ";"
                                        << "ext" << ";"
                                        << propOffspringDead/count2 << ";"
                                        << propDeath/count4 << ";"
                                        << propDeathExtrinsic/count3 << ";"
                                        << propDeathMutation/count3 << ";"
                                        << quantiles[1] << ";"
                                        << quantiles[2] << ";"
                                        << quantiles[3] << ";"
                                        << quantiles[0] << ";"
                                        << quantiles[4] << ";"
                                        << quantilesDamage[1] << ";"
                                        << quantilesDamage[2] << ";"
                                        << quantilesDamage[3] << ";"
                                        << quantilesDamage[0] << ";"
                                        << quantilesDamage[4] << std::endl;
                    }

                    population.updateNewTimeStep();         // step 1 at each time step
                    population.replaceDeadIndividuals();    // step 2 at each time step
                }
                ++t;
            }

            // We also record information at the end of the simulation
            if(ExtinctionPop == false){

                // Below all information is saved as before; see comments in the portion of the code above
                Population populationSave = population;
                std::vector<int> ageDeathVec = {};
                std::vector<int> damageDeathVec = {};
                double popSize(0.0);
                double propOffspringDead(0.0);
                double count(0.0);
                double count2(0.0);
                double count3(0.0);
                double count4(0.0);

                double propDeathExtrinsic(0.0);
                double propDeathMutation(0.0);
                double propDeath(0.0);

                for (int rep2(0); rep2<10; ++rep2){
                    population = populationSave;
                    for (int t(0);t<1e3;++t) {

                        population.updateNewTimeStep();

                        // life span
                        for (auto iter:population._ListInd) {
                            count4+=1.0;
                            if(iter._livingState==false){
                                propDeath+=1.0;
                                count3+=1.0;
                                ageDeathVec.push_back(iter._age);
                                damageDeathVec.push_back(iter._damage);
                                if(iter._causeDeath==0){            // death extrinsic mortality
                                    propDeathExtrinsic+=1.0;
                                }else if(iter._causeDeath==1){      // death mutation
                                    propDeathMutation+=1.0;
                                }
                            }
                        }

                        population.replaceDeadIndividuals();

                        // pop size
                        popSize += (double) population._ListInd.size();

                        // prop offspring
                        if(population._ListIndOfspringIndex.size()>0){
                            propOffspringDead += ((double) population._ListIndOfspringIndex.size() - ((double) population._ListInd.size() - population._ListLivingInd.size())) / ((double) population._ListIndOfspringIndex.size()) ;
                            count2+=1.0;
                        }
                        count+=1.0;
                    }

                }
                population = populationSave;

                std::vector<double> ageDeathVecDouble(ageDeathVec.begin(), ageDeathVec.end());
                auto quantiles = Quantile<double>(ageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1  });


                std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1  });

                dataSensitivity << birthRatePerYear << ";"
                                << mortRatePerYear << ";"
                                << rateAccumul << ";"
                                << t << ";"
                                << popSize/count << ";"
                                << "no ext" << ";"
                                << propOffspringDead/count2 << ";"
                                << propDeath/count4 << ";"
                                << propDeathExtrinsic/count3 << ";"
                                << propDeathMutation/count3 << ";"
                                << quantiles[1] << ";"
                                << quantiles[2] << ";"
                                << quantiles[3] << ";"
                                << quantiles[0] << ";"
                                << quantiles[4] << ";"
                                << quantilesDamage[1] << ";"
                                << quantilesDamage[2] << ";"
                                << quantilesDamage[3] << ";"
                                << quantilesDamage[0] << ";"
                                << quantilesDamage[4] << std::endl;
            }
        }
    }
    // We close the data set
    dataSensitivity.close();
}


    ///###### MAIN

int main()
{
   srand ( time(NULL) );

    // Parameters

    int carryingCapacity = 100;                     // Carrying capacity; default value: 500
    bool densityDependenceSurvival = false;         // Whether survival is adjusted when N<CarryingCapacity after recruitment

    double alphaMax = 1.0;                          // Maximum fecundity when pleiotropy (when mutation expressed at age == 0)
    double rateAlphaFecundityYear = 50     /365;    // Rate at which fecundity decreases with the age at which intrinsic mortality occurs (called gamma in the manuscript)

    double convertIntoYear = 12;                   // Discretization of one year into time steps
                                                    // Depending on the simulation: if time steps = months: 12 ; if time steps = days: 365
                                                    // Considering that time steps are months instead of days decreases simulation runtime

    int rangeEffectDeleteriousMutation = 1;         // Age at which deleterious mutations can have an effect (if=1: in months if convertIntoYear=12, and in days if convertIntoYear=365)

    int maxAgeConsideredYear = 60;                  // Maximum age considered (in years); to determine the size of vectors (should be high enough)
    if(abs(convertIntoYear-365.0)<1){               // To speed up simulations, this parameter is adjusted depending on the type of time steps considered
        maxAgeConsideredYear = 3;
    }else if(abs(convertIntoYear-12.0)<1){
        maxAgeConsideredYear = 30;
    }

    bool startAtEarlierLifeSpan = true;             // If true, implement a lethal mutation expressed at a late age
                                                    // It decreases simulation runtime because the first lethal mutation expressed at very late age
                                                    // take a long time to fixate
    double quantileStart = 0.5;                     // If startAtEarlierLifeSpan==true, it defines the quantile of lifespan at which the fist lethal mutation is expressed
                                                    // default value: 0.8

    double probDeleteriousMutationPerAge = 2e-3;    // Probability of getting a deleterious mutation expressed at a given age

    bool boolReverseMutation = false;               // if true, beneficial 'reverse' mutations can occur
    double ratioReverseMutation = 1.0;              // if boolReverseMutation==true, proportion of beneficial mutations relative to deleterious mutations

    double effectSurvDeleteriousMutation = 1.0;     // Survival probability reduced due to the expression of a single deleterious mutation (if ==1, lethal mutation)

    std::string typeAccumulation = "lin";           // Type of accumulation of damage: "lin", "exp", "log", or "randomlin"
                                                    // "lin" is implemented in most simulations (it reflects the age), "randomlin" is implemented when we consider stochastic
                                                    // change of somatic state

    double rateAccumul = 1;                         // if _typeAccumulation=="lin", link between somatic state and age (=1)
                                                    // if _typeAccumulation=="randomlin", rate of accumulation of somatic states (called 'damage' in the script

    // Characteristics sensitivity analysis

    int Tburnin = 2e3;                              // Burn-in period to assess population extinction before mutation accumulations (in years); default value: 2e3
    int Tmax = 1e5;                                 // Simulation sime (in years); default value: 1e6

    std::string Namefile = "exampleMonth_convertIntoYear365rangeEffect1";

    int rep = 1;                                    // replicate index ; in demo simulation, rep = 1 and 2; in manuscript, default value:rep = 1 to 10
    Namefile=Namefile+"_Rep"+std::to_string(rep);   // we increment the replicate index to the name of the file

    double nbValuesTested = 5;                     // Number of parameter values tested; default value: 20

    double minBirthRatePerYear = 0.1;               // Minimum birth rate (per year) tested
    double maxBirthRatePerYear = 2.0;               // Maximum birth rate (per year) tested

    double minMortRatePerYear = 0.1;                // Minimum adult mortality rate (per year) tested
    double maxMortRatePerYear = 2.0;                // Maximum adult mortality rate (per year) tested


    // Conversion from year to month (if each time step = month; if convertIntoYear==12) or to day (if each time step = day; if convertIntoYear==365)

    std::vector<double> vecBirthRatePerYear;
    double rateConst2 = minBirthRatePerYear;
    while (rateConst2<maxBirthRatePerYear+ 1e-5) {
        vecBirthRatePerYear.push_back(rateConst2);
        rateConst2 += (maxBirthRatePerYear-minBirthRatePerYear)/(nbValuesTested-1);
    }
    std::vector<double> vecMortRatePerYear;
    double rateConst = minMortRatePerYear;
    while (rateConst<maxMortRatePerYear+ 1e-5) {
        vecMortRatePerYear.push_back(rateConst);
        rateConst += (maxMortRatePerYear-minMortRatePerYear)/(nbValuesTested-1);
    }
    double rateAlphaFecundity = rateAlphaFecundityYear / convertIntoYear;
    int maxAgeConsidered = maxAgeConsideredYear*convertIntoYear;

    // From the maximum age considered (in time steps) to the maximum damage considered (equality is assumed given that we account for linear accumulation)
    int maxDamageConsidered(0);
    maxDamageConsidered = maxAgeConsidered;

    // Sensitivity analysis

    SensitivityAnalysis(            carryingCapacity,
                                    densityDependenceSurvival,
                                    maxDamageConsidered,
                                    convertIntoYear,
                                    startAtEarlierLifeSpan,
                                    quantileStart,
                                    probDeleteriousMutationPerAge,
                                    typeAccumulation,
                                    rateAccumul,
                                    ratioReverseMutation,
                                    rangeEffectDeleteriousMutation,
                                    boolReverseMutation,
                                    effectSurvDeleteriousMutation,
                                    Tburnin,
                                    Tmax,
                                    vecBirthRatePerYear,
                                    alphaMax,
                                    rateAlphaFecundity,
                                    vecMortRatePerYear,
                                    Namefile
                           );

    return 0;

}
