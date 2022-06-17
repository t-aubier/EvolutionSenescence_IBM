setwd("~/Projets/Senescence_IBM/")

library(ggplot2)
library(cowplot)
library(plyr)
library(stringr)

theme_set(theme_cowplot())


####################### CONTINUOUS DAMAGE ACCUMULATION

nameFileVec = c("Run01b_Rep1_maxAge100QuantStart09_Mut1e3Range1K200DdepF")




#### MONTH AND YEAR ####

# Tmax 5e6 ; Range 1 12 24 ; K 200 500 1000 ; Rep 1 -> 10
# nameFileVec = "Run02_K500DdepFmaxAge50QuantStart09_Mut1e3Range12"
# divisionYear = 12
# K = 500

# Tmax 5e6 ; Range 1 12 ; K 500 ; Rep 1 -> 10 ; alpha 2 5 10
# nameFileVec = "Run02Alpha10Rate5_K500DdepFmaxAge50QuantStart09_Mut1e3Range12"
# divisionYear = 12
# K = 500

# Tmax 5e6 ; Range 1 12 ; K 500 ; Rep 1 -> 10 ; alpha 2 5 10
# nameFileVec = "Run02Alpha1Rate5_reverseTratio1over1_K500DdepFmaxAge50QuantStart099_Mut1e3Range1"
# divisionYear = 12
# K = 500


# Tmax 5e6 ; month bis ; Range 1 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; reverse Ratio 1/1      range 1 pour mois / 12 pour année !!!!!!!!
# nameFileVec = "Run02monthbisAlpha1Rate5_reverseTratio1over1_K500DdepFmaxAge50QuantStart099_Mut1e3Range12"
# divisionYear = 12
# K = 500

# Tmax 5e6 ; Range 1 12 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; reverse Ratio 1/1
# nameFileVec = "Run03Alpha1Rate5_reverseTratio1over1_K500DdepFmaxAge50QuantStart09_Mut1e3Range1"
# divisionYear = 12
# K = 500


# Tmax 5e6 ; month ; Range 1 12 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; Long simulations Tmax = 10e6 au lieu de 5e6
# nameFileVec = "Run02monthLongAlpha1Rate5_K500DdepFmaxAge50QuantStart09_Mut1e3Range12"
# divisionYear = 12
# K = 500


#### DAY ####

# Tmax 5e6 ; day ; Range 1 ; K 500 ; Rep 1 -> 10 ; alpha 2 5 10
# nameFileVec = "Run02dayAlpha5Rate10_K500DdepFmaxAge50QuantStart09_Mut1e3Range1"
# divisionYear = 365
# K = 500



# Tmax 5e6 ; day ; Range 1 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; reverse Ratio 1/1
# nameFileVec = "Run02dayAlpha1Rate5_reverseTratio1over1_K500DdepFmaxAge30QuantStart099_Mut1e3Range1"
# divisionYear = 365
# K = 500

# Tmax 5e6 ; day bis ; Range 1 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; reverse Ratio 1/1
# nameFileVec = "Run02daybisAlpha1Rate5_reverseTratio1over1_K500DdepFmaxAge30QuantStart099_Mut1e3Range1"
# divisionYear = 365
# K = 500

# Tmax 5e6 ; Range 1 12 ; K 500 ; Rep 1 -> 10 ; alpha 1 ; reverse Ratio 1/1
# nameFileVec = "Run03dayAlpha1Rate5_reverseTratio1over1_K500DdepFmaxAge50QuantStart09_Mut1e3Range1"
# divisionYear = 12
# K = 500





# Tmax 1e7 ; Range 1 12 ; K 200 500 ; Rep 1 2
# nameFileVec = c("Run03_K500DdepFmaxAge50QuantStart09_Mut1e3Range12")
# divisionYear = 12

# test avec deviation from linearity lin exp ;  Tmax 5e6 ; Range 1 ; K 200 500 1000; Rep 1 2
# nameFileVec = c("Run02b_K500DdepFmaxAge50QuantStart09_Mut1e3Range1_exp") 



# test avec density dependence T F; Range 1 12 ; K 200; Rep 1
# nameFileVec = c("Run02_K200DdepFmaxAge50QuantStart09_Mut1e3Range1_lin")
# nameFileVec = c("Run02_K200DdepTmaxAge50QuantStart09_Mut1e3Range1_lin_Rep1")



# nameFileVec = c("TESTRun02_K200DdepFmaxAge50QuantStart09_Mut1e3Range1_Rep1")






##### NEW SIMULATIONS

### A - Basic simulations 

# Tmax 1e6 ; 12range1 12range12 ; K 200 500 1000     (UZH cluster for K = 200 and K = 1000)
nameFileVec = "RunNew03_scale12range12_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 12
K = 500
alphaMax = 1
# delta = 1        # day
delta = 365/12   # month = ~30days
delta = 365      # year


# Tmax 1e6 ; 365range1 ; K 200 500 1000              (UZH cluster for K = 200 and K = 1000)
# nameFileVec = "RunNew03_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500


### B - Reverse mutations (K =500)

# Tmax 1e6 ; 365range1              (UZH cluster)
# nameFileVec = "RunNew03_scale365range1_Alpha1Rate50_reverseTratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500

# Tmax 1e6 ; 12range12              (UZH cluster)
# nameFileVec = "RunNew03_scale12range1_Alpha1Rate50_reverseTratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500


### C - Alpha 5 50                

# Tmax 1e6 ; 365range1              (UNC cluster)
# nameFileVec = "RunNew03_scale365range1_Alpha50Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500

# Tmax 1e6 ; 12range1 12range12             (UNC cluster)
# nameFileVec = "RunNew03_scale12range12_Alpha50Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500

### C - Somatic (rate 10 et rate 1 rep 1->5)           Faire rate 0.1 (pret), rate 100 ? 

# Tmax 1e6 ; 12range1              (UNC cluster)
nameFileVec = "RunNew03_SomaticRate10_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 12
K = 500

# nameFileVec = "RunNew03_SomaticRate1_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500



## Somatic (rate  1 5 10 100 rep 1->5) WITH MAX LIFE SPAN POSSIBLE THAT IS HIGHER; WITH MONTH
nameFileVec = "RunNew04_SomaticRate100_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 12
K = 500
# 
# nameFileVec = "RunNew04_SomaticRate01_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500

## Somatic (rate 10 100 rep 1->5) WITH DAY
nameFileVec = "RunNew04_SomaticRate10_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 365
K = 500


### C - Density dependence ;  K 500

## Without density dependence

nameFileVec = "RunNew05DD_F_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 365
K = 500

# nameFileVec = "RunNew05DD_F_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500


## With density dependence

nameFileVec = "RunNew05DD_T_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 365
K = 500

nameFileVec = "RunNew05DD_T_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 12
K = 500

## Weaker density dependence

# nameFileVec = "RunNew05DD_T09_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500

# nameFileVec = "RunNew05DD_T09_scale12range12_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500

##  density dependence without threshold

nameFileVec = "RunNew05DD_noT_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 365
K = 500

nameFileVec = "RunNew05DD_noT_scale12range12_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
divisionYear = 12
K = 500

### D - Long time series for K = 500 (UZH cluster)

# Tmax 2e6 ; 12range1 12range12     (UZH cluster)
# nameFileVec = "RunNew03Long_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500

  

### E - Small Mutation; effect on surv = 0.01 instead of 1

# Tmax 1e6 ; 365range1
# nameFileVec = "RunNew03_smallMut001_scale365range1_Alpha1Rate50_reverseTratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500

# Tmax 1e6 ; 12range1 12range12    
# nameFileVec = "RunNew03_smallMut001_scale12range12_Alpha1Rate50_reverseTratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500

# Tmax 1e6 ; 365range1; no rev mut
# nameFileVec = "RunNew03_smallMut001_scale365range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 365
# K = 500

# Tmax 1e6 ; 12range1 12range12; no rev mut
# nameFileVec = "RunNew03_smallMut001_scale12range12_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
# divisionYear = 12
# K = 500






#####
maxLifeSpanShown = 25
maxLifeSpanShown = 25
nbRep = 30
#######################

firstFile = 1
filesMissing <- c()
for (i in 1:length(nameFileVec)){
  
  nameOneFile = paste("Data/dataSensitivity_",nameFileVec[i],".csv",sep="")
  if(file.exists(nameOneFile)){
    dataSensitivityOne <- read.csv(nameOneFile, sep=";")
    dataSensitivityOne$rep = NA
    
    if(length(which(names(dataSensitivityOne)=="RateAccumul"))==0){
      dataSensitivityOne$RateAccumul = 1
    }
    if(length(which(names(dataSensitivityOne)=="Damage0025"))==0){
      dataSensitivityOne$Damage0025 = dataSensitivityOne$LifeSpan0025
    }
    if(length(which(names(dataSensitivityOne)=="Damage05"))==0){
      dataSensitivityOne$Damage05 = dataSensitivityOne$LifeSpan05
    }
    if(length(which(names(dataSensitivityOne)=="Damage0975"))==0){
      dataSensitivityOne$Damage0975 = dataSensitivityOne$LifeSpan0975
    }
    if(length(which(names(dataSensitivityOne)=="DamageMin"))==0){
      dataSensitivityOne$DamageMin = dataSensitivityOne$LifeSpanMin
    }
    if(length(which(names(dataSensitivityOne)=="DamageMax"))==0){
      dataSensitivityOne$DamageMax = dataSensitivityOne$LifeSpanMax
    }
    if(firstFile==1){
      dataSensitivity = dataSensitivityOne
      firstFile = 0  
    }else{
      dataSensitivity = rbind(dataSensitivity,dataSensitivityOne)
    }
  }else{
    filesMissing <- c(filesMissing,nameOneFile)
  }
  
  for (rep in 1:nbRep){
    nameOneFile = paste("Data/dataSensitivity_",nameFileVec[i],"_Rep",rep,".csv",sep="")
    if(file.exists(nameOneFile)){
      dataSensitivityOne <- read.csv(nameOneFile, sep=";")
      
      dataSensitivityOne$rep = rep
      
      if(length(which(names(dataSensitivityOne)=="RateAccumul"))==0){
        dataSensitivityOne$RateAccumul = 1
      }
      if(length(which(names(dataSensitivityOne)=="Damage0025"))==0){
        dataSensitivityOne$Damage0025 = dataSensitivityOne$LifeSpan0025
      }
      if(length(which(names(dataSensitivityOne)=="Damage05"))==0){
        dataSensitivityOne$Damage05 = dataSensitivityOne$LifeSpan05
      }
      if(length(which(names(dataSensitivityOne)=="Damage0975"))==0){
        dataSensitivityOne$Damage0975 = dataSensitivityOne$LifeSpan0975
      }
      if(length(which(names(dataSensitivityOne)=="DamageMin"))==0){
        dataSensitivityOne$DamageMin = dataSensitivityOne$LifeSpanMin
      }
      if(length(which(names(dataSensitivityOne)=="DamageMax"))==0){
        dataSensitivityOne$DamageMax = dataSensitivityOne$LifeSpanMax
      }
      if(firstFile==1){
        dataSensitivity = dataSensitivityOne
        firstFile = 0  
      }else{
        dataSensitivity = rbind(dataSensitivity,dataSensitivityOne)
      }
    }else{
      filesMissing <- c(filesMissing,nameOneFile)
    }
  }
}
summary(dataSensitivity)


convertRateYearToDay = 1/365
convertDayToYear = 1/365
# convert in year: expected lifespan based on analytical derivation
dataSensitivity$AnalyticalMinMaxLifeSpan = 
  1/(dataSensitivity$MortRatePerYear*convertRateYearToDay)   *
  log(
    (alphaMax * K - exp(-(dataSensitivity$MortRatePerYear*convertRateYearToDay) * delta) * (K - 1)) / (K * (alphaMax-1)+1)  
  ) * convertDayToYear

dataSensitivity$AnalyticalMaxMaxLifeSpan = 
  1/(dataSensitivity$MortRatePerYear*convertRateYearToDay) *
  log(
    (K*(alphaMax*exp((dataSensitivity$MortRatePerYear*convertRateYearToDay) * delta)-1)+1) / (K * (alphaMax-1)+1)
  ) * convertDayToYear

# convert in year: lifespan in simulations
dataSensitivity$LifeSpan0975day = dataSensitivity$LifeSpan0975 * delta
dataSensitivity$LifeSpan0975 <- dataSensitivity$LifeSpan0975/divisionYear
dataSensitivity$JussiX <- log(K)/dataSensitivity$MortRatePerYear
dataSensitivity$ratioMaxLifeSpanJussi <- dataSensitivity$LifeSpan0975/dataSensitivity$JussiX

dataSensitivity$diffAnalyticalMinMaxLifeSpan = dataSensitivity$LifeSpan0975 - dataSensitivity$AnalyticalMinMaxLifeSpan
dataSensitivity$diffAnalyticalMaxMaxLifeSpan = dataSensitivity$LifeSpan0975 - dataSensitivity$AnalyticalMaxMaxLifeSpan

# ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
#   geom_point(aes(AnalyticalMinMaxLifeSpan,LifeSpan0975))+
#   xlim(0,40)+
#   ylim(0,40)+
#   geom_abline(intercept = 0)
# ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
#   geom_point(aes(AnalyticalMaxMaxLifeSpan,LifeSpan0975))+
#   xlim(0,40)+
#   ylim(0,40)+
#   geom_abline(intercept = 0)

# ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
#   geom_point(aes(LifeSpan0975,AnalyticalMinMaxLifeSpan),color="blue")+
#   geom_point(aes(LifeSpan0975,AnalyticalMaxMaxLifeSpan),color="red",size=0.5)+
#   xlim(0,40)+
#   ylim(0,40)+
#   geom_abline(intercept = 0)


ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
  geom_point(aes(AnalyticalMinMaxLifeSpan,LifeSpan0975),color="red")+
  xlim(0,40)+
  ylim(0,40)+
  geom_abline(intercept = 0,lty=2)
  

dataSensitivitySum <- 
  ddply(dataSensitivity, 
        .(BirthRatePerYear,MortRatePerYear), 
        summarise, 
        LifeSpan0975 = mean(LifeSpan0975, na.rm=T),
        JussiX = mean(JussiX, na.rm=T),
        ratioMaxLifeSpanJussi= mean(ratioMaxLifeSpanJussi, na.rm=T)
  )
# dataSensitivity <-ddply(dataSensitivity, 
#                         .(MortRatePerYear, BirthRatePerYear,PriorProbAskMutant,PriorProbGiveMutant,Age,
#                                       Param), summarize, 
#                         Psurv = mean(Psurv,na.rm=T), 
#                         Psurvtoage = mean(Psurvtoage,na.rm=T),
#                         diffPsurv = mean(diffPsurv,na.rm=T))
#   
# dataSensitivity$LifeSpan0975 <- dataSensitivity$LifeSpan0975/divisionYear
# dataSensitivitySum$LifeSpan0975 <- dataSensitivitySum$LifeSpan0975/divisionYear

# dataSensitivitySum$JussiX <- log(K)/dataSensitivitySum$MortRatePerYear
# dataSensitivitySum$ratioMaxLifeSpanJussi <- dataSensitivitySum$LifeSpan0975/dataSensitivitySum$JussiX

summary(dataSensitivity)
dataSensitivity$ExtStatus2 <- dataSensitivity$ExtStatus
dataSensitivity$ExtStatus2[dataSensitivity$ExtStatus2=="no ext"] <- NA

dataSensitivity$ExtStatus <- factor(dataSensitivity$ExtStatus, levels=c("no ext","ext","ini ext"))
dataSensitivity$ExtStatus2 <- factor(dataSensitivity$ExtStatus2, levels=c("ext","ini ext"))

step1 = (unique(dataSensitivity$MortRatePerYear)[2]-unique(dataSensitivity$MortRatePerYear)[1])/2
step2 = (unique(dataSensitivity$BirthRatePerYear)[2]-unique(dataSensitivity$BirthRatePerYear)[1])/2

step1 = 0.05
step2 = 0.05

colVec <- c("blue","gray50","black")
p1 <- ggplot(dataSensitivity)+
  geom_rect(aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
                fill=ExtStatus,color=ExtStatus))+
  scale_color_manual(name="",breaks=c("no ext","ext","ini ext"),values=colVec,
                     labels=c("no ext.","ext. during mut. accumul.","ext. before mut. accumul."))+
  scale_fill_manual(name="",breaks=c("no ext","ext","ini ext"),values=colVec,
                    labels=c("no ext.","ext. during mut. accumul.","ext. before mut. accumul."))+
  scale_x_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  # xlim(0.05,2.05)+
  
  # ylim(0.05,2.05)+      
  xlab("Extrinsic mortality rate")+
  ylab("Birth rate")




p2 <- ggplot(dataSensitivity)+
  geom_rect(aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
                fill=LifeSpan0975,color=LifeSpan0975))+
  scale_color_continuous("Max life span (in years)",limits=c(0,30),na.value=colVec[3])+
  scale_fill_continuous("Max life span (in years)",limits=c(0,30),na.value=colVec[3])+
  scale_x_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  # xlim(0.05,2.05)+
  # ylim(0.05,2.05)+
  xlab("Extrinsic mortality rate")+
  ylab("Birth rate")

dataSensitivity2 <- dataSensitivity

dataSensitivityExtIni = dataSensitivity[dataSensitivity$ExtStatus=="ini ext",]
dataSensitivityExtIni$coord = paste(dataSensitivityExtIni$BirthRatePerYear,dataSensitivityExtIni$MortRatePerYear,sep="_")
dataSensitivityExtIni2 = data.frame(BirthRatePerYear=numeric(),MortRatePerYear=numeric())
for(i in unique(dataSensitivityExtIni$coord)){
  dataSensitivityExtIni2=rbind(dataSensitivityExtIni2,data.frame(BirthRatePerYear=as.numeric(str_split(i,"_")[[1]][1]),MortRatePerYear=as.numeric(str_split(i,"_")[[1]][2])))
}

dataSensitivityExt = dataSensitivity[dataSensitivity$ExtStatus=="ext",]
dataSensitivityExt$coord = paste(dataSensitivityExt$BirthRatePerYear,dataSensitivityExt$MortRatePerYear,sep="_")
dataSensitivityExt2 = data.frame(BirthRatePerYear=numeric(),MortRatePerYear=numeric())
for(i in unique(dataSensitivityExt$coord)){
  dataSensitivityExt2=rbind(dataSensitivityExt2,data.frame(BirthRatePerYear=as.numeric(str_split(i,"_")[[1]][1]),MortRatePerYear=as.numeric(str_split(i,"_")[[1]][2])))
}

p2b <- ggplot()+
  geom_rect(data=dataSensitivitySum[is.na(dataSensitivitySum$LifeSpan0975)==FALSE,],
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2,
                fill=LifeSpan0975,color=LifeSpan0975))+
  scale_color_continuous("Max life span (in years)",type = "viridis",limits=c(0,maxLifeSpanShown),na.value=colVec[3])+
  scale_fill_continuous("Max life span (in years)",type = "viridis",limits=c(0,maxLifeSpanShown),na.value=colVec[3])+
  geom_rect(data=dataSensitivityExt2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[2],color=colVec[2])+
  geom_rect(data=dataSensitivityExtIni2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[3],color=colVec[3])+
  scale_x_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  # xlim(0.05,2.05)+
  # ylim(0.05,2.05)+
  xlab("Extrinsic mortality rate")+
  ylab("Birth rate")

p2b_Jussi <- ggplot()+
  geom_rect(data=dataSensitivitySum,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
                fill=JussiX,color=JussiX))+
  scale_color_continuous("Max life span (in years)",type = "viridis",limits=c(0,maxLifeSpanShown*2.5),na.value=colVec[3])+
  scale_fill_continuous("Max life span (in years)",type = "viridis",limits=c(0,maxLifeSpanShown*2.5),na.value=colVec[3])+
  geom_rect(data=dataSensitivityExt2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[2],color=colVec[2])+
  geom_rect(data=dataSensitivityExtIni2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[3],color=colVec[3])+
  scale_x_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  # xlim(0.05,2.05)+
  # ylim(0.05,2.05)+
  xlab("Extrinsic mortality rate")+
  ylab("Birth rate")


p2b_diffJussi <- ggplot()+
  geom_rect(data=dataSensitivitySum,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
                fill=ratioMaxLifeSpanJussi,color=ratioMaxLifeSpanJussi))+
  scale_color_continuous(name=expression(paste("Max life span\nrelative to\n",ln(Ne)/mu,sep="")),type = "viridis",option="plasma",limits=c(0,0.63),na.value=colVec[3])+
  scale_fill_continuous(name=expression(paste("Max life span\nrelative to\n",ln(Ne)/mu,sep="")),type = "viridis",option="plasma",limits=c(0,0.63),na.value=colVec[3])+
  geom_rect(data=dataSensitivityExt2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[2],color=colVec[2])+
  geom_rect(data=dataSensitivityExtIni2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[3],color=colVec[3])+
  scale_x_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0.05,2.05),breaks = c(0.1,0.5,1,1.5,2))+
  # xlim(0.05,2.05)+
  # ylim(0.05,2.05)+
  xlab(expression(paste("Extrinsic mortality rate (",mu,")",sep="")))+
  ylab("Birth rate")


tail(dataSensitivity)
plotAll <- plot_grid(p1,p2,p2b,nrow=3)

save_plot("Graphes/singlePlotSensitivity.pdf",
          p2b+theme(legend.position="none"),
          base_height = 5, base_aspect_ratio = 1)


save_plot("Graphes/singlePlotSensitivityJussi.pdf",
          p2b_diffJussi+theme(legend.position="none"),
          base_height = 5, base_aspect_ratio = 1)


save_plot("Graphes/singlePlotSensitivityLegend.pdf",
          p2b,
          base_height = 5, base_aspect_ratio = 1)


save_plot("Graphes/singlePlotSensitivityJussiLegend.pdf",
          p2b_diffJussi,
          base_height = 5, base_aspect_ratio = 1)


save_plot("Graphes/plotSensitivity.pdf",
          plotAll,
          base_height = 4, base_aspect_ratio = 3)
# 
# p3 <- ggplot(dataSensitivity)+
#   geom_rect(aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
#                 fill=RateAccumul,color=RateAccumul))+
#   xlab("Extrinsic mortality rate")+
#   ylab("Birth rate")+
#   ylim(0,2.1)

# plot_grid(p2b,p2b_Jussi,p2b_diffJussi,nrow=3)

max(dataSensitivity$LifeSpan0975[dataSensitivity$ExtStatus=="no ext"] , na.rm=T)
tail(dataSensitivityOne)
filesMissing

max(dataSensitivity$LifeSpan0975,na.rm=T)
dim(dataSensitivityOne)[1] 
( (20*20*10)-dim(dataSensitivityOne)[1]) 
dim(dataSensitivityOne)[1] / (20*20*10)

plotAll
  