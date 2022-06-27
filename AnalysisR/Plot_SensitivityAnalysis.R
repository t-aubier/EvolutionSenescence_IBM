##### Directory

setwd("/home/taubier/Projets/Github/EvolutionSenescence_IBM")     # change the working directory

##### Packages

library(scales)
library(ggplot2)
library(cowplot)
library(plyr)
library(stringr)
theme_set(theme_cowplot())

####################### PLOT THE SENSITIVITY ANALYSIS ##########  


##### Parameters

# simulation with monthly-expressed mutations (see parameters implemented in Sensitivity_Analysis.cpp):
nameFileVec = c("exampleMonth_convertIntoYear365rangeEffect1")               # Name of the simulation output
divisionYear = 12                                                            # Discretization of one year into time steps, as in the simulation output
                                                                             # Depending on the simulation: if time steps = months: 12 ; if time steps = days: 365

# simulation with daily-expressed mutations (same parameters as in Sensitivity_Analysis.cpp, but with convertIntoYear = 365):
# nameFileVec = c("exampleDay_convertIntoYear365rangeEffect1")
# divisionYear = 365

K = 200                                                                       # Carrying capacity as in the simulation

maxLifeSpanShown = 25                                                         # The maximum lifespan shown in the the plot (color scale)
nbRep = 3                                                                     # Maximum number of replicates tested (that is fine if the replicate simulations do not exist)


##### Process data and plot

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

dataSensitivity$LifeSpan0975 <- dataSensitivity$LifeSpan0975/divisionYear                       # convert in year: lifespan in simulations
dataSensitivity$LehtonenX <- log(K)/dataSensitivity$MortRatePerYear
dataSensitivity$ratioMaxLifeSpanLehtonen <- dataSensitivity$LifeSpan0975/dataSensitivity$LehtonenX

dataSensitivitySum <- 
  ddply(dataSensitivity, 
        .(BirthRatePerYear,MortRatePerYear), 
        summarise, 
        LifeSpan0975 = mean(LifeSpan0975, na.rm=T),
        LehtonenX = mean(LehtonenX, na.rm=T),
        ratioMaxLifeSpanLehtonen= mean(ratioMaxLifeSpanLehtonen, na.rm=T)
      )

dataSensitivity$ExtStatus2 <- dataSensitivity$ExtStatus
dataSensitivity$ExtStatus2[dataSensitivity$ExtStatus2=="no ext"] <- NA
dataSensitivity$ExtStatus <- factor(dataSensitivity$ExtStatus, levels=c("no ext","ext","ini ext"))
dataSensitivity$ExtStatus2 <- factor(dataSensitivity$ExtStatus2, levels=c("ext","ini ext"))

step2 = (unique(dataSensitivity$BirthRatePerYear)[2]-unique(dataSensitivity$BirthRatePerYear)[1])/2
# step1 = (unique(dataSensitivity$MortRatePerYear)[2]-unique(dataSensitivity$MortRatePerYear)[1])/2   # bug if sensitivity analyses not done
step1 = step2


colVec <- c("blue","gray50","black")

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

plot <- ggplot()+
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
  scale_x_continuous(limits=c(0-step1*1.01,2+step1*1.01),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0-step2*1.01,2+step2*1.01),breaks = c(0.1,0.5,1,1.5,2))+
  xlab("Adult mortality rate")+
  ylab("Birth rate")

plot_diffLehtonen <- ggplot()+
  geom_rect(data=dataSensitivitySum,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2, 
                fill=ratioMaxLifeSpanLehtonen,color=ratioMaxLifeSpanLehtonen))+
  scale_color_continuous(name=expression(paste("Max life span\nrelative to\n",ln(Ne)/mu,sep="")),type = "viridis",option="plasma",limits=c(0,0.63),na.value=colVec[3])+
  scale_fill_continuous(name=expression(paste("Max life span\nrelative to\n",ln(Ne)/mu,sep="")),type = "viridis",option="plasma",limits=c(0,0.63),na.value=colVec[3])+
  geom_rect(data=dataSensitivityExt2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[2],color=colVec[2])+
  geom_rect(data=dataSensitivityExtIni2,
            aes(xmin=MortRatePerYear-step1, xmax=MortRatePerYear+step1, ymin=BirthRatePerYear-step2, ymax=BirthRatePerYear+step2),
            fill=colVec[3],color=colVec[3])+
  scale_x_continuous(limits=c(0-step1*1.01,2+step1*1.01),breaks = c(0.1,0.5,1,1.5,2))+
  scale_y_continuous(limits=c(0-step2*1.01,2+step2*1.01),breaks = c(0.1,0.5,1,1.5,2))+
  xlab(expression(paste("Adult mortality rate (",mu,")",sep="")))+
  ylab("Birth rate")


plotAll <- plot_grid(plot,plot_diffLehtonen,ncol=1)

save_plot("Graphes/singlePlotSensitivity.pdf",
          plot,
          base_height = 5, base_aspect_ratio = 1.5)


save_plot("Graphes/singlePlotSensitivityLehtonen.pdf",
          plot_diffLehtonen,
          base_height = 5, base_aspect_ratio = 1.5)

plotAll

  