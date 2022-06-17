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

template<typename T>
static inline double Lerp(T v0, T v1, T t){
    return (1 - t)*v0 + t*v1;
}

template<typename T>
static inline std::vector<T> Quantile(const std::vector<T>& inData, const std::vector<T>& probs){
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
    
void SensitivityAnalysis        (   int carryingCapacity,
                                    bool densityDependenceSurvival,
                                    int maxDamageConsidered,
                                    double convertIntoYear,
                                    bool startAtEarlierLifeSpan,
                                    double quantileStart,
                                    double probDeleteriousMutationPerAge,
                                    std::string typeAccumulation,
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
    
//     print(Namefile);
    
    std::fstream dataSensitivity;
    dataSensitivity.open ("Data/dataSensitivity_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
    dataSensitivity << "BirthRatePerYear;" << "MortRatePerYear;" << "RateAccumul;" << "Time;" << "SizePop;"  << "ExtStatus;"  << "PercentOffspringDead;" << "propDeath;" << "propDeathExtrinsic;" << "propDeathMutation;" << "LifeSpan0025;" << "LifeSpan05;" << "LifeSpan0975;" << "LifeSpanMin;" << "LifeSpanMax;" << "Damage0025;" << "Damage05;" << "Damage0975;" << "DamageMin;" << "DamageMax"  << std::endl;

//     std::fstream dataMutationAccumulation;
//     dataMutationAccumulation.open ("Data/dataMut_"+Namefile+".csv", std::fstream::in | std::fstream::out | std::fstream::trunc);
//     dataMutationAccumulation << "Rep;" << "Time;" << "Age;" << "SurvivalMutation0025;" << "SurvivalMutation05;" << "SurvivalMutation0975;" << "MeanSurvivalMutation" << std::endl;
    
    
    double rateAccumul(1.0);
    
    
    for (auto& mortRatePerYear:vecMortRatePerYear) {
        
        double probSurvExtrinsicMortality = exp(-mortRatePerYear/convertIntoYear);
        
        
        

        // we determine the rate of accumulation if non linear accumulation.
        
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
            
            
            
            
            
            
            /*
            
            
            Population population = Population( carryingCapacity,
                                                densityDependenceSurvival,
                                                nbOffspringPerInd,
                                                maxDamageConsidered,
                                                probSurvExtrinsicMortality,
                                                1.0,
                                                "lin",
                                                0.0,                            // prob of getting a deleterious mutation = 0.0 here
                                                convertIntoYear,
                                                rangeEffectDeleteriousMutation,
                                                effectSurvDeleteriousMutation
                                            );
            // burn in phase
            for (int t(0);t<Tburnin*12;++t) {
                population.updateNewTimeStep();
                population.replaceDeadIndividuals();
            }
        

            // Save data pop at t=0
            int tmaxlinear=1e3;
            int nbRepLinear=20;
            Population populationSave = population;
            std::vector<int> damageDeathVec = {}; 
            for (int rep2(0); rep2<nbRepLinear; ++rep2){
                population = populationSave;
                std::vector<int> damageDeathVec2 = {}; 
                for (int t(0);t<tmaxlinear;++t) {
                    population.updateNewTimeStep();
                    population.replaceDeadIndividuals();
                }
                
                for (int t(0);t<tmaxlinear;++t) {
                    population.updateNewTimeStep();
                    // life span 
                    for (auto iter:population._ListInd) {
                        damageDeathVec2.push_back(iter._damage);
                        
                    }
                    population.replaceDeadIndividuals();
                }
                            
                std::vector<double> damageDeathVecDouble2(damageDeathVec2.begin(), damageDeathVec2.end());
                auto quantilesDamage2 = Quantile<double>(damageDeathVecDouble2, {0.5});
                
                damageDeathVec.push_back(quantilesDamage2[0]);
                
                
            }
            population = populationSave;
            
            std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
            auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0.5 });
            double MedianDamageLevelLinear = quantilesDamage[0];
            
            
            
            int tmax1 = tmaxlinear;
            int nbRep1 = nbRepLinear/10;
            rateAccumul = 0;
            double rateAccumulInterv;
            if(typeAccumulation=="log"){
                rateAccumulInterv = 10^200;
            }else if(typeAccumulation=="exp"){
                rateAccumulInterv = 0.01;
            }
            double MedianDamageLevelLinearNonLinear = 0;
            double epsilon = 1;
            while(abs(MedianDamageLevelLinearNonLinear-MedianDamageLevelLinear)>epsilon){
                while(MedianDamageLevelLinearNonLinear<MedianDamageLevelLinear-epsilon){
                        
                    rateAccumul = rateAccumul + rateAccumulInterv;
                    
                    // Here we define the median damage level when deviation from linearity
                    
                    Population population = Population( carryingCapacity,
                                                        densityDependenceSurvival,
                                                        nbOffspringPerInd,
                                                        maxDamageConsidered,
                                                        probSurvExtrinsicMortality,
                                                        rateAccumul,
                                                        typeAccumulation,
                                                        0.0,                            // prob of getting a deleterious mutation = 0.0 here
                                                        convertIntoYear,
                                                        rangeEffectDeleteriousMutation,
                                                        effectSurvDeleteriousMutation
                                                    );
                    // burn in phase
                    for (int t(0);t<Tburnin*12;++t) {
                        population.updateNewTimeStep();
                        population.replaceDeadIndividuals();
                    }

                    // Save data pop at t=0
                    Population populationSave = population;
                    
                    std::vector<int> damageDeathVec = {}; 
                    for (int rep2(0); rep2<nbRep1; ++rep2){
                        population = populationSave;
                        std::vector<int> damageDeathVec2 = {}; 
                        for (int t(0);t<tmax1;++t) {
                            population.updateNewTimeStep();
                            population.replaceDeadIndividuals();
                        }
                        
                        for (int t(0);t<tmax1;++t) {
                            population.updateNewTimeStep();
                            // life span 
                            for (auto iter:population._ListInd) {
                                damageDeathVec2.push_back(iter._damage);
                                
                            }
                            population.replaceDeadIndividuals();
                        }
                                    
                        std::vector<double> damageDeathVecDouble2(damageDeathVec2.begin(), damageDeathVec2.end());
                        auto quantilesDamage2 = Quantile<double>(damageDeathVecDouble2, {0.5});
                        
                        damageDeathVec.push_back(quantilesDamage2[0]);
                        
                        
                    }
                    population = populationSave;
                    
                    std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                    auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0.5  });
                    
                    MedianDamageLevelLinearNonLinear = quantilesDamage[0];     
                    
                    print(rateAccumul);
                    print(MedianDamageLevelLinear);
                    print(MedianDamageLevelLinearNonLinear);
                    print(rateAccumulInterv);
                    print(nbRep1);
                    print("");
                    
                    if(MedianDamageLevelLinearNonLinear>MedianDamageLevelLinear/2){
                       nbRep1 = nbRep1*2;
                       if(nbRep1>nbRepLinear){
                           nbRep1 = nbRepLinear;
                       }
                    }
                    if(abs(MedianDamageLevelLinearNonLinear-MedianDamageLevelLinear)<20){
                       nbRep1 = nbRepLinear;
                    }
                    if(MedianDamageLevelLinearNonLinear>MedianDamageLevelLinear){
                        rateAccumul -= rateAccumulInterv;
                        MedianDamageLevelLinearNonLinear = 0;
                        rateAccumulInterv = rateAccumulInterv/10;
                    }
                }
            }
            */

            
        }
/*             print("rateAccumul :");
        print(rateAccumul);     */           
        
        
        
        
        for (auto& birthRatePerYear:vecBirthRatePerYear) {
                    
            
//             double nbOffspringPerInd = 1-exp(-birthRatePerYear/convertIntoYear);   
            
            double convertIntoYearDAY = 365;
            double probSurvExtrinsicMortalityDAY = exp(-mortRatePerYear/convertIntoYearDAY);
            double nbOffspringPerIndDAY = 1-exp(-birthRatePerYear/convertIntoYearDAY);
            double nbOffspringPerInd = nbOffspringPerIndDAY * convertIntoYearDAY/convertIntoYear  / probSurvExtrinsicMortality;
            
            
            
//             print(birthRatePerYear);
//             print(mortRatePerYear);
        
            // population without limit in maxDamageConsidered
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
            // burn in phase
            bool ExtinctionPop = false;
            int t(0);
            while(t<Tburnin*convertIntoYear+1 && ExtinctionPop == false) {
                if(population._ListInd.size()==0){
                    ExtinctionPop = true;
                }else{
                    population.updateNewTimeStep();
                    population.replaceDeadIndividuals();
                }
                ++t;
            }
            
            if(ExtinctionPop==true){
//                 print("Initial Extinction");
//                 print("");

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
            }else{
                t = 0;
                
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
                        
                        // cause of death; A FAIRE
                        
                        count+=1.0;
                    }
                    
                }
                population = populationSave;
                
                std::vector<double> ageDeathVecDouble(ageDeathVec.begin(), ageDeathVec.end());
                auto quantiles = Quantile<double>(ageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });
                
                std::vector<double> damageDeathVecDouble(damageDeathVec.begin(), damageDeathVec.end());
                auto quantilesDamage = Quantile<double>(damageDeathVecDouble, { 0, 0.025, 0.5, 0.975, 1, quantileStart  });

    //             dataSensitivity << rep << ";" << t << ";" 
    //                     << popSize/count << ";" 
    //                     << propOffspringDead/count2 << ";" 
    //                     << propDeath/count4 << ";" 
    //                     << propDeathExtrinsic/count3 << ";" 
    //                     << propDeathMutation/count3 << ";" 
    //                     << quantiles[1] << ";"
    //                     << quantiles[2] << ";"
    //                     << quantiles[3] << ";" 
    //                     << quantiles[0] << ";" 
    //                     << quantiles[4] << std::endl;

                // Save data mutation accumulation
                        
                                        
    //             int age(0);
    //             while (age<maxDamageConsidered){
    //                 std::vector<double> SurvMutVec = {}; 
    //                 for (auto iter:population._ListInd) {
    //                     if(age<iter._survCausedPerMutation.size()){
    //                         SurvMutVec.push_back(iter._survCausedPerMutation[age]);
    //                     }else{                            
    //                         SurvMutVec.push_back(0.0);
    //                     }
    //                 }
    // 
    //                 auto quantiles2 = Quantile<double>(SurvMutVec, { 0.025, 0.5, 0.975 });
    //                 if(mean(SurvMutVec)>1e-20){
    //                     dataMutationAccumulation << rep << ";" << t << ";" 
    //                                                 << age << ";" 
    //                                                 << quantiles2[0] << ";" 
    //                                                 << quantiles2[1] << ";" 
    //                                                 << quantiles2[2] << ";" 
    //                                                 << mean(SurvMutVec) << std::endl;
    //                 }
    //                 age+=rangeEffectDeleteriousMutation;
    //             }


                
                
                
                
                    
                if(startAtEarlierLifeSpan==true){
                    // get median lifespan

                    int newMaxDamageConsidered = (int) (round(quantilesDamage[5])+1);
        //             print("new max age considered:");
        //             print(newMaxDamageConsidered);
                    
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
                            population.updateNewTimeStep();
                            population.replaceDeadIndividuals();
                        }
                        ++t2;
                    }
                    if(ExtinctionPop==true){
                        ExtinctionPop = false;
                        population = populationSave;
                    }
                                
                }
                
                
                
                // add probability to have a deleterious mutation
                population.setProbDeleteriousMutationPerAge(probDeleteriousMutationPerAge);
                
                
                
                
            
            }
            
            
            // actual simulation
            bool popGettingExtinct = false;
            ++t;
            while(t<Tmax*convertIntoYear+1 && ExtinctionPop == false) {
                if(population._ListInd.size()==0){
                    ExtinctionPop = true;
                }else{
                    if(t% ((int) 1e4)==0){
                        population.shortenMutationVector(boolReverseMutation);
                    }
                                        
                    if(((double) population._ListInd.size())/((double) carryingCapacity)<0.1 && popGettingExtinct==false){
                        popGettingExtinct=true;
                        ExtinctionPop = true;
//                         print("Final Extinction");
//                         print("");

                        
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
                                
                                // cause of death; A FAIRE
                                
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

                        // Save data mutation accumulation
    //                     int age(0);
    //                     while (age<maxAgeConsidered){
    //                         std::vector<double> SurvMutVec = {}; 
    //                         for (auto iter:population._ListInd) {
    //                             if(age<iter._survCausedPerMutation.size()){
    //                                 SurvMutVec.push_back(iter._survCausedPerMutation[age]);
    //                             }else{                            
    //                                 SurvMutVec.push_back(0.0);
    //                             }
    //                         }
    // 
    //                         auto quantiles2 = Quantile<double>(SurvMutVec, { 0.025, 0.5, 0.975 });
    //                         if(mean(SurvMutVec)>1e-20){
    //                             dataMutationAccumulation << rep << ";" << t << ";" 
    //                                                         << age << ";" 
    //                                                         << quantiles2[0] << ";" 
    //                                                         << quantiles2[1] << ";" 
    //                                                         << quantiles2[2] << ";" 
    //                                                         << mean(SurvMutVec) << std::endl;
    //                         }
    //                         age+=rangeEffectDeleteriousMutation;
    //                     }
                        
                    }
                    
                    
                    population.updateNewTimeStep();
                    population.replaceDeadIndividuals();
                                        
                }

                ++t;
            }
            
            
            if(ExtinctionPop == false){
                    
//                 print("Maintenance");
//                 print("");

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
                        
                        // cause of death; A FAIRE
                        
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

                // Save data mutation accumulation
//                 int age(0);
//                 while (age<maxAgeConsidered){
//                     std::vector<double> SurvMutVec = {}; 
//                     for (auto iter:population._ListInd) {
//                         if(age<iter._survCausedPerMutation.size()){
//                             SurvMutVec.push_back(iter._survCausedPerMutation[age]);
//                         }else{                            
//                             SurvMutVec.push_back(0.0);
//                         }
//                     }
// 
//                     auto quantiles2 = Quantile<double>(SurvMutVec, { 0.025, 0.5, 0.975 });
//                     if(mean(SurvMutVec)>1e-20){
//                         dataMutationAccumulation << rep << ";" << t << ";" 
//                                                     << age << ";" 
//                                                     << quantiles2[0] << ";" 
//                                                     << quantiles2[1] << ";" 
//                                                     << quantiles2[2] << ";" 
//                                                     << mean(SurvMutVec) << std::endl;
//                     }
//                     age+=rangeEffectDeleteriousMutation;
//                 }
            }
            
       
                        
            
        }

                        
    }
    
    dataSensitivity.close();
//     dataMutationAccumulation.close();
}



    ///###### MAIN

int main()
{   
   srand ( time(NULL) ); 
    
    // Parameters
    
    int carryingCapacity = 500;
    bool densityDependenceSurvival = false;
    
    double convertIntoYear = 365;                  // if time steps = months: 12 ; if time steps = days: 365
    int rangeEffectDeleteriousMutation = 1;        // in months if convertIntoYear 12; in days if convertIntoYear 365
    
    int maxAgeConsideredYear = 60;
    if(abs(convertIntoYear-365.0)<1){
        maxAgeConsideredYear = 3;
    }else if(abs(convertIntoYear-12.0)<1){
        maxAgeConsideredYear = 30;
    }
    
    
    bool startAtEarlierLifeSpan = true;
    double quantileStart = 0.8;
    
    double probDeleteriousMutationPerAge = 2e-3;
    bool boolReverseMutation = true;
    double ratioReverseMutation = 1.0/1.0;
    
    double effectSurvDeleteriousMutation = 1.0;

    std::string typeAccumulation = "lin";    
    
    // Characteristics sensitivity analysis
    
    int Tburnin = 2e3;      // in years ; 2
    int Tmax = 1e6;         // in years ; 5e6
    
    double alphaMax = 1.0;
    double rateAlphaFecundityYear = 50     /365;
    
    std::string Namefile = "RunNew03_scale365range1_Alpha1Rate50_reverseTratio1over1_K500DdepFQuantStart08_Mut2e3";
    
    int rep = 1;
    Namefile = Namefile+"_Rep"+std::to_string(rep);
    
    
    double nbValuesTested = 20;        //20

    double minBirthRatePerYear = 0.1;  //0.1
    double maxBirthRatePerYear = 2.0;  
    
    double minMortRatePerYear = 0.1;  
    double maxMortRatePerYear = 2.0;  
    

    
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
    
//     print(vecBirthRatePerYear);
//     print(vecMortRatePerYear);

    // Conversion from year to month
    
    int maxAgeConsidered = maxAgeConsideredYear*convertIntoYear;  
    
    int maxDamageConsidered(0);
    maxDamageConsidered = maxAgeConsidered;   
//     if(typeAccumulation=="lin"){
//         maxDamageConsidered = (int) (maxAgeConsidered * rateAccumul);        
//     }else if(typeAccumulation=="exp"){
//         maxDamageConsidered = (int) (exp(maxAgeConsidered * rateAccumul));
//     }else if(typeAccumulation=="log"){
//         maxDamageConsidered = (int) (log(maxAgeConsidered * rateAccumul));
//     }
    
    
    // Sensitivity analysis 
    
    SensitivityAnalysis(            carryingCapacity,
                                    densityDependenceSurvival,
                                    maxDamageConsidered,
                                    convertIntoYear,
                                    startAtEarlierLifeSpan,
                                    quantileStart,
                                    probDeleteriousMutationPerAge,
                                    typeAccumulation,
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
