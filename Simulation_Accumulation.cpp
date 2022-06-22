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

#include <iomanip>

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

    ///###### SENSITIVITY ANALYSIS -- PROBABILITY OF INVASION OF DELETERIOUS MUTATIONS

void ReplicateSimulation        (   int carryingCapacity,
                                    bool densityDependenceSurvival,
                                    double nbOffspringPerInd,
                                    double alphaMax,
                                    double rateAlphaFecundity,
                                    int maxDamageConsidered,
                                    double convertIntoYear,
                                    bool startAtEarlierLifeSpan,
                                    double quantileStart,
                                    double probSurvExtrinsicMortality,
                                    double rateAccumul,
                                    std::string typeAccumulation,
                                    double probDeleteriousMutationPerAge,
                                    double ratioReverseMutation,
                                    int rangeEffectDeleteriousMutation,
                                    bool boolReverseMutation,
                                    double effectSurvDeleteriousMutation,
                                    int NbRep,
                                    int Tburnin,
                                    int Tmax,
                                    int Tstep,
                                    std::string Namefile
                           ){

    print(Namefile);

    // the files where information will be stored

    std::fstream dataPop;
    dataPop.open ("Data/dataPop_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataPop << "Rep;" << "Time;" << "SizePop;" << "PercentOffspringDead;" << "propDeath;" << "propDeathExtrinsic;" << "propDeathMutation;" << "LifeSpan0025;" << "LifeSpan05;" << "LifeSpan0975;" << "LifeSpanMin;" << "LifeSpanMax;" << "Damage0025;" << "Damage05;" << "Damage0975;" << "DamageMin;" << "DamageMax"  << std::endl;

    std::fstream dataMutationAccumulation;
    dataMutationAccumulation.open ("Data/dataMut_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataMutationAccumulation << "Rep;" << "Time;" << "Damage;" << "SurvivalMutation0025;" << "SurvivalMutation05;" << "SurvivalMutation0975;" << "MeanSurvivalMutation" << std::endl;

    int indCharacSaved_rep = 0;

    for(int rep(0); rep<NbRep; ++rep){

        // population without limit in maxDamageConsidered
        Population population = Population( carryingCapacity,
                                            densityDependenceSurvival,
                                            nbOffspringPerInd,
                                            alphaMax,
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

        // burn in phase
        for (int t(0);t<Tburnin*convertIntoYear;++t) {
            population.updateNewTimeStep();
            population.replaceDeadIndividuals();
        }

        int t = 0;
        print(t);

        // Save data pop at t=0
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
        auto quantiles = Quantile<double>(ageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });

        std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
        auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });

        dataPop << rep << ";" << t << ";"
                << popSize/count << ";"
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

        // Save data mutation accumulation

        int damage(0);
        while (damage<maxDamageConsidered){
            std::vector<double> SurvMutVec = {};
            for (auto iter:population._ListInd) {
                if(damage<iter._survCausedPerMutation.size()){
                    SurvMutVec.push_back(iter._survCausedPerMutation[damage]);
                }else{
                    if(iter._survCausedPerMutation[iter._survCausedPerMutation.size()-1]>1-0.001){
                        SurvMutVec.push_back(1.0);
                    }else{
                        SurvMutVec.push_back(0.0);
                    }
                }
            }

            auto quantiles2 = Quantile<double>(SurvMutVec, { 0.025, 0.5, 0.975 });
            if(mean(SurvMutVec)>1e-20){
                dataMutationAccumulation << rep << ";" << t << ";"
                                            << damage << ";"
                                            << quantiles2[0] << ";"
                                            << quantiles2[1] << ";"
                                            << quantiles2[2] << ";"
                                            << mean(SurvMutVec) << std::endl;
            }
            damage+=rangeEffectDeleteriousMutation;
        }
        int tCheck=(Tstep*convertIntoYear);

        if(startAtEarlierLifeSpan==true){
            // get median lifespan

            int newMaxDamageConsidered = (int) (round(quantilesDamage[5])+1);

            // population with limit in maxDamageConsidered
            population = Population( carryingCapacity,
                                     densityDependenceSurvival,
                                     nbOffspringPerInd,
                                     alphaMax,
                                     rateAlphaFecundity,
                                     newMaxDamageConsidered,
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
            // burn in phase
            for (int t(0);t<Tburnin*convertIntoYear;++t) {
                population.updateNewTimeStep();
                population.replaceDeadIndividuals();
            }

        }

        // add probability to have a deleterious mutation
        population.setProbDeleteriousMutationPerAge(probDeleteriousMutationPerAge);

        // actual simulation
        bool ExtinctionPop = false;
        bool popGettingExtinct = false;
        ++t;
        while(t<Tmax*convertIntoYear+1 && ExtinctionPop == false) {

            if(t% ((int) 1e4)==0){
                population.shortenMutationVector(boolReverseMutation);
            }

            if(population._ListInd.size()==0){

                std::cout << "EXTINCTION AT T = " <<  t<< std::endl;
                dataPop << rep << ";" << t << ";"
                        << population._ListInd.size() << ";"
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
                ExtinctionPop = true;

            }else if((((double) population._ListInd.size())/((double) carryingCapacity)<0.1 && popGettingExtinct==false) || t==tCheck){
                if(((double) population._ListInd.size())/((double) carryingCapacity)<0.1){
                    popGettingExtinct=true;
                }
                print(t);

                // Save data pop
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
                auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });

                dataPop << rep << ";" << t << ";"
                        << popSize/count << ";"
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

                // Save data mutation accumulation
                int damage(0);
                while (damage<maxDamageConsidered){
                    std::vector<double> SurvMutVec = {};
                    for (auto iter:population._ListInd) {
                        if(damage<iter._survCausedPerMutation.size()){
                            SurvMutVec.push_back(iter._survCausedPerMutation[damage]);
                        }else{
                            if(iter._survCausedPerMutation[iter._survCausedPerMutation.size()-1]>1-0.01){
                                SurvMutVec.push_back(1.0);
                            }else{
                                SurvMutVec.push_back(0.0);
                            }
                        }
                    }

                    auto quantiles2 = Quantile<double>(SurvMutVec, { 0.025, 0.5, 0.975 });
                    if(mean(SurvMutVec)>1e-20){

                        dataMutationAccumulation << rep << ";" << t << ";"
                                                    << damage << ";"
                                                    << quantiles2[0] << ";"
                                                    << quantiles2[1] << ";"
                                                    << quantiles2[2] << ";"
                                                    << mean(SurvMutVec) << std::endl;
                    }
                    damage+=rangeEffectDeleteriousMutation;
                }
                tCheck+=(Tstep*convertIntoYear);
            }
            population.updateNewTimeStep();
            population.replaceDeadIndividuals();
            ++t;
        }
    }
    dataPop.close();
    dataMutationAccumulation.close();
}



    ///###### MAIN

int main()
{
   srand ( time(NULL) );

    // Parameters

    int carryingCapacity = 500;                     // Carrying capacity
    bool densityDependenceSurvival = false;         // Whether survival is adjusted when N<CarryingCapacity after recruitment

    double birthRatePerYear = 1.5;                  // Birth rate per year; converted into a number of offspring per time step below

    double alphaMax = 1.0;                          // Maximum fecundity when pleiotropy (when mutation expressed at age == 0)
    double rateAlphaFecundityYear = 50     /365;    // Rate at which fecundity decreases with the age at which intrinsic mortality occurs (called gamma in the manuscript)

    double mortRatePerYear = 1.0;                   // Mortality rate per year; converted into a number of offspring per time step below
    int maxAgeConsideredYear = 30;                  // Maximum age considered (in years); to determine the size of vectors (should be high enough)

    double convertIntoYear = 12.0;                  // Discretization of one year into time steps
                                                    // Depending on the simulation: if time steps = months: 12 ; if time steps = days: 365
                                                    // Considering that time steps are months instead of days decreases simulation runtime

    bool startAtEarlierLifeSpan = true;             // If true, implement a lethal mutation expressed at a late age
                                                    // It decreases simulation runtime because the first lethal mutation expressed at very late age
                                                    // take a long time to fixate
    double quantileStart = 0.9;                     // If startAtEarlierLifeSpan==true, it defines the quantile of lifespan at which the fist lethal mutation is expressed

    double probDeleteriousMutationPerAge = 2e-3;    // Probability of getting a deleterious mutation expressed at a given age
    int rangeEffectDeleteriousMutation = 1;         // Age at which deleterious mutations can have an effect (if=1: in months if convertIntoYear=12, and in days if convertIntoYear=365)

    bool boolReverseMutation = false;               // if true, beneficial 'reverse' mutations can occur
    double ratioReverseMutation = 1;                // if boolReverseMutation==true, proportion of beneficial mutations relative to deleterious mutations

    double effectSurvDeleteriousMutation = 1;       // Survival probability reduced due to the expression of a single deleterious mutation (if ==1, lethal mutation)

    std::string typeAccumulation = "lin";           // Type of accumulation of damage: "lin", "exp", "log", or "randomlin"
                                                    // "lin" is implemented in most simulations (it reflects the age), "randomlin" is implemented when we consider stochastic
                                                    // change of somatic state

    double rateAccumul = 1;                         // if _typeAccumulation=="lin", link between somatic state and age (=1)
                                                    // if _typeAccumulation=="randomlin", rate of accumulation of somatic states (called 'damage' in the script)



    // Characteristics sensitivity analysis

    int NbRep = 1;                                  // Number of replicates

    int Tburnin = 2;                                // Burn-in period to assess population extinction before mutation accumulations (in years)
    int Tmax = 100000;                              // Simulation sime (in years) ; 3e6
    int Tstep = 1000;                               // Time at which information is recoreded (in years) ; 2e4

    std::string Namefile = "Run12somaticMonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin";

    // Conversion from year to month (if each time step = month; if convertIntoYear==12) or to day (if each time step = day; if convertIntoYear==365)

    double probSurvExtrinsicMortality = exp(-mortRatePerYear/convertIntoYear);
    double convertIntoYearDAY = 365;
    double probSurvExtrinsicMortalityDAY = exp(-mortRatePerYear/convertIntoYearDAY);
    double nbOffspringPerIndDAY = 1-exp(-birthRatePerYear/convertIntoYearDAY);
    double nbOffspringPerInd = nbOffspringPerIndDAY * convertIntoYearDAY/convertIntoYear  / probSurvExtrinsicMortality;
    int maxAgeConsidered = maxAgeConsideredYear*convertIntoYear;
    double rateAlphaFecundity = rateAlphaFecundityYear / convertIntoYear;

    // From the maximum age considered (in time steps) to the maximum damage considered (equality is assumed given that we account for linear accumulation)
    int maxDamageConsidered(0);
    maxDamageConsidered = maxAgeConsidered;

    // Sensitivity analysis

    time_t start, end;

    // Recording of start time
    time(&start);

    ReplicateSimulation(            carryingCapacity,
                                    densityDependenceSurvival,
                                    nbOffspringPerInd,
                                    alphaMax,
                                    rateAlphaFecundity,
                                    maxDamageConsidered,
                                    convertIntoYear,
                                    startAtEarlierLifeSpan,
                                    quantileStart,
                                    probSurvExtrinsicMortality,
                                    rateAccumul,
                                    typeAccumulation,
                                    probDeleteriousMutationPerAge,
                                    ratioReverseMutation,
                                    rangeEffectDeleteriousMutation,
                                    boolReverseMutation,
                                    effectSurvDeleteriousMutation,
                                    NbRep,
                                    Tburnin,
                                    Tmax,
                                    Tstep,
                                    Namefile
                           );

   // Recording end time.
    time(&end);

    // Calculating total time taken by the program.
    double time_taken = double(end - start);
    std::cout << "Time taken by program is : " << std::fixed
         << time_taken << std::setprecision(5);
    std::cout << " sec " << std::endl;

    return 0;

}
