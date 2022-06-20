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

    std::fstream dataPop;
    dataPop.open ("Data/dataPop_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataPop << "Rep;" << "Time;" << "SizePop;" << "PercentOffspringDead;" << "propDeath;" << "propDeathExtrinsic;" << "propDeathMutation;" << "LifeSpan0025;" << "LifeSpan05;" << "LifeSpan0975;" << "LifeSpanMin;" << "LifeSpanMax;" << "Damage0025;" << "Damage05;" << "Damage0975;" << "DamageMin;" << "DamageMax"  << std::endl;

    std::fstream dataMutationAccumulation;
    dataMutationAccumulation.open ("Data/dataMut_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataMutationAccumulation << "Rep;" << "Time;" << "Damage;" << "SurvivalMutation0025;" << "SurvivalMutation05;" << "SurvivalMutation0975;" << "MeanSurvivalMutation" << std::endl;

    int indCharacSaved_rep = 0;

    for(int rep(0); rep<NbRep; ++rep){

        // we assess the rate of damage accumulation that matches the median damage level when linear accumulation
        if(typeAccumulation!="lin"){

            int nbInd = 5e5;

            //linear accumulation
            std::string func = "lin";

            std::vector<int> damageDeathVec = {};
            int ind(0);
            int damage;
            while(ind<nbInd){
                ++ind;
                bool dead = false;
                int age = 0 ;
                while(dead==false){
                    if(func=="exp"){
                        damage = (int) exp(age * rateAccumul)-1;
                    }else if(func=="lin"){
                        damage = age;
                    }
                    damageDeathVec.push_back(damage);

                    if(randomDouble()>probSurvExtrinsicMortality){
                        dead = true;
                    }
                    age = age+1;
                }
            }

            std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
            auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0.5 });
            double MedianDamageLevelLinear = quantilesDamage[0];

            rateAccumul = 0;
            double rateAccumulInterv;
            if(typeAccumulation=="log"){
                rateAccumulInterv = 10^200;
            }else if(typeAccumulation=="exp"){
                rateAccumulInterv = 0.01;
            }
            double MedianDamageLevelLinearNonLinear = 0;
            double epsilon = 0.5;
            while(abs(MedianDamageLevelLinearNonLinear-MedianDamageLevelLinear)>epsilon){
                while(MedianDamageLevelLinearNonLinear<MedianDamageLevelLinear-epsilon){

                    rateAccumul = rateAccumul + rateAccumulInterv;

                    // Here we define the median damage level when deviation from linearity

                    std::string func = typeAccumulation;

                    std::vector<int> damageDeathVec = {};
                    int ind(0);
                    int damage;
                    while(ind<nbInd){
                        ++ind;
                        bool dead = false;
                        int age = 0 ;
                        while(dead==false){
                            if(func=="exp"){
                                damage = (int) exp(age * rateAccumul)-1;
                            }else if(func=="lin"){
                                damage = age;
                            }
                            damageDeathVec.push_back(damage);

                            if(randomDouble()>probSurvExtrinsicMortality){
                                dead = true;
                            }
                            age = age+1;
                        }
                    }

                    std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                    auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0.5 });
                    MedianDamageLevelLinearNonLinear = quantilesDamage[0];

                    print(rateAccumul);
                    print(MedianDamageLevelLinear);
                    print(MedianDamageLevelLinearNonLinear);
                    print(rateAccumulInterv);
                    print("");

                    if(MedianDamageLevelLinearNonLinear>MedianDamageLevelLinear){
                        rateAccumul -= rateAccumulInterv;
                        MedianDamageLevelLinearNonLinear = 0;
                        rateAccumulInterv = rateAccumulInterv/10;
                    }
                }
            }
        }
        print("rateAccumul :");
        print(rateAccumul);

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

    int carryingCapacity = 500;
    bool densityDependenceSurvival = false;

    double birthRatePerYear = 1.5;

    double alphaMax = 1.0;
    double rateAlphaFecundityYear = 50     /365;

    double mortRatePerYear = 1.0;
    int maxAgeConsideredYear = 30;                  // in years : typically  10

    double convertIntoYear = 12.0;                   // if time steps = months: 12 ; if time steps = days: 365

    bool startAtEarlierLifeSpan = true;
    double quantileStart = 0.9;     //0.8

    double probDeleteriousMutationPerAge = 2e-3;
    int rangeEffectDeleteriousMutation = 1;        // in months if convertIntoYear 12; in days if convertIntoYear 365
    bool boolReverseMutation = false;
    double ratioReverseMutation = 1;

    double effectSurvDeleteriousMutation = 1;
    
    std::string typeAccumulation = "lin";
    double rateAccumul = 1;                 // if _typeAccumulation=="randomlin", rate of accumulation of somatic states

    // Characteristics sensitivity analysis

    int NbRep = 1;

    int Tburnin = 2;        // in years ; 2
    int Tmax = 100000;         // in years ; 3e6
    int Tstep = 1000;        // in years ; 2e4

    std::string Namefile = "Run12somaticMonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin";

    // Conversion from year to month

    double probSurvExtrinsicMortality = exp(-mortRatePerYear/convertIntoYear);

    double convertIntoYearDAY = 365;
    double probSurvExtrinsicMortalityDAY = exp(-mortRatePerYear/convertIntoYearDAY);
    double nbOffspringPerIndDAY = 1-exp(-birthRatePerYear/convertIntoYearDAY);
    double nbOffspringPerInd = nbOffspringPerIndDAY * convertIntoYearDAY/convertIntoYear  / probSurvExtrinsicMortality;

    int maxAgeConsidered = maxAgeConsideredYear*convertIntoYear;

    double rateAlphaFecundity = rateAlphaFecundityYear / convertIntoYear;
    int maxDamageConsidered(0);
    maxDamageConsidered = maxAgeConsidered;

    // Sensitivity analysis

    time_t start, end;
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
