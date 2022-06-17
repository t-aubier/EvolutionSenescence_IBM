setwd("~/Projets/Senescence_IBM/")

library(ggplot2)
library(cowplot)
library(plyr)
library(stringr)

theme_set(theme_cowplot())


####################### 

nameFileVec = c("RunNew03_scale12range1_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3","RunNew03_scale12range12_Alpha1Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3",
                "RunNew03_scale12range1_Alpha1Rate50_reverseFratio1over1_K1000DdepFQuantStart08_Mut2e3","RunNew03_scale12range12_Alpha1Rate50_reverseFratio1over1_K1000DdepFQuantStart08_Mut2e3",
                "RunNew03_scale12range1_Alpha5Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3","RunNew03_scale12range12_Alpha5Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3",
                "RunNew03_scale12range1_Alpha10Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3","RunNew03_scale12range12_Alpha10Rate50_reverseFratio1over1_K500DdepFQuantStart08_Mut2e3"
                
                )
divisionYear = c(12,12,
                 12,12,
                 12,12,
                 12,12)
K=c(500,500,
    1000,1000,
    500,500,
    500,500)
alphaMax=c(1,1,
           1,1,
           5,5,
           10,10)
delta=c(365/12,365,
        365/12,365,
        365/12,365,
        365/12,365)




#####
maxLifeSpanShown = 25
# maxLifeSpanShown = 40
nbRep = 10
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
    
    
    
    
    dataSensitivityOne$AnalyticalMinMaxLifeSpan = 
      1/dataSensitivityOne$MortRatePerYear *log(
        (alphaMax[i] * K[i] - exp(-dataSensitivityOne$MortRatePerYear*delta[i]) * (K[i] - 1)) / (K[i] * (alphaMax[i]-1)+1)    / 365
      )
    dataSensitivityOne$AnalyticalMaxMaxLifeSpan = 
      1/dataSensitivityOne$MortRatePerYear *log(
        (K[i]*(alphaMax[i]*exp(dataSensitivityOne$MortRatePerYear*delta[i])-1)+1) / (K[i] * (alphaMax[i]-1)+1)             / 365
      )    
    
    dataSensitivityOne$LifeSpan0975 <- dataSensitivityOne$LifeSpan0975/divisionYear[i]
    
    
    
    
    
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
      
      
      
      dataSensitivityOne$AnalyticalMinMaxLifeSpan = 
        1/dataSensitivityOne$MortRatePerYear *log(
          (alphaMax[i] * K[i] - exp(-dataSensitivityOne$MortRatePerYear*delta[i]) * (K[i] - 1)) / (K[i] * (alphaMax[i]-1)+1)    / 365
        )
      dataSensitivityOne$AnalyticalMaxMaxLifeSpan = 
        1/dataSensitivityOne$MortRatePerYear *log(
          (K[i]*(alphaMax[i]*exp(dataSensitivityOne$MortRatePerYear*delta[i])-1)+1) / (K[i] * (alphaMax[i]-1)+1)             / 365
        )    
      
      dataSensitivityOne$LifeSpan0975 <- dataSensitivityOne$LifeSpan0975/divisionYear[i]
      
      
      
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

ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
  geom_point(aes(LifeSpan0975,AnalyticalMinMaxLifeSpan),color="blue")+
  geom_point(aes(LifeSpan0975,AnalyticalMaxMaxLifeSpan),color="red")+
  xlim(0,40)+
  ylim(0,40)+
  geom_abline(intercept = 0)

ggplot(dataSensitivity[dataSensitivity$ExtStatus=="no ext",])+
  geom_point(aes(diffAnalyticalMinMaxLifeSpan,diffAnalyticalMaxMaxLifeSpan))+
  # xlim(-40,40)+
  # ylim(-40,40)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)

# save_plot("Graphes/plotSensitivity.pdf",
#           plotAll,
#           base_height = 4, base_aspect_ratio = 3)
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

