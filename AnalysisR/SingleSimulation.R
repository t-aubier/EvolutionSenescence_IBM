setwd("~/Projets/Senescence_IBM/")

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


####################### CONTINUOUS DAMAGE ACCUMULATION

# PC
nameFileVec = c("Run02_maxAge6MedT_Mut1e6K200DdepF_extrMort1birth2")
nameFileVec = c("Run02_maxAge6MedT_Mut1e4K200DdepT_extrMort1birth2")


# IMAC   Mut 1e4 1e5  ;  K 200 500 ; dep T F ; b=2
# IMAC   Mut 1e4 1e5  ;  K 200 ; dep T F ; b = 4 10
nameFileVec = c("Run02_maxAge6QuantStart09_Mut1e4K200DdepF_extrMort1birth4_Rep1")



nameFileVec = c("Run01_Rep1_maxAge6QuantStart09_Mut1e5K200DdepT_extrMort1birth4")

nameFileVec = c("Run01_Rep1_maxAge6QuantStart05_Mut1e5K1000DdepF_extrMort1birth4")
# nameFileVec = c("Run01_Rep1_maxAge6QuantStart04_Mut1e5K1000DdepF_extrMort1birth4")

nameFileVec = c("Run04_Rep1_maxAge50QuantStart09_Mut1e4Range12K200DdepF_extrMort01birth1")
nameFileVec = c("Run04_Rep1_maxAge50QuantStart09_Mut1e3Range1K200DdepF_extrMort01birth1")

nameFileVec = c("Run04_Rep1_maxAge50QuantStart09_Mut1e3Range1K200DdepF_extrMort01birth05")

nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth1extrMort0.5exp")    # rateAccumul  0.0505
nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth1extrMort0.5bexp")    # rateAccumul 0.051

nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth1extrMort05exp")    # rateAccumul  ?
# nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth1extrMort05lin")


nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort06exp")    
nameFileVec = c("Run05_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort06lin")


nameFileVec = c("Run06_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort06lin")
divisionYear = 12
divisionTime = 1000

nameFileVec = c("Run06dayRangeDay_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort06lin")
nameFileVec = c("Run06dayRangeMonth_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort06lin")
divisionYear = 365
divisionTime = 1000

nameFileVec = c("Run06dayRangeDay_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort02lin")
nameFileVec = c("Run06dayRangeMonth_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort02lin")
divisionYear = 365
divisionTime = 1000


nameFileVec = c("Run06MonthRangeMonthAlphaMax1_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort02lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax50_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort02lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax100_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort02lin")

# nameFileVec = c("Run06MonthRangeMonthAlphaMax1_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax50_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax100_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin")


nameFileVec = c("Run06MonthRangeMonthAlphaMax1_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
nameFileVec = c("Run06MonthRangeMonthAlphaMax10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
nameFileVec = c("Run06MonthRangeMonthAlphaMax50_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax100_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin")


nameFileVec = c("Run06MonthRangeMonthAlphaMax1Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax2Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax5Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")
# nameFileVec = c("Run06MonthRangeMonthAlphaMax10Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin")

nameFileVec ="Run07MonthRangeMonthAlphaMax1Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin"
nameFileVec ="Run07bMonthRangeMonthAlphaMax1Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin"

nameFileVec ="Run08b50MonthRangeMonthAlphaMax1Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin"
nameFileVec ="Run08b50DayRangeDayAlphaMax1Rate5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort1lin"

nameFileVec ="Run08MonthRangeMonthAlphaMax1Rate5_ratioRev10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
nameFileVec ="Run08MonthRangeMonthAlphaMax5Rate5_ratioRev10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
# nameFileVec ="Run08MonthRangeMonthAlphaMax2Rate5_ratioRev10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"


nameFileVec ="Run08MonthRangeMonthAlphaMax5Rate5_ratioRev5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
nameFileVec ="Run08MonthRangeMonthAlphaMax5Rate5_ratioRev10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"

nameFileVec ="Run08MonthRangeMonthAlphaMax1Rate5_ratioRev5_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
# nameFileVec ="Run08MonthRangeMonthAlphaMax1Rate5_ratioRev10_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"

nameFileVec ="Run08b1MonthRangeYearAlphaMax1Rate5_ratioRev1_Rep1_maxAge50QuantStart099_Mut1e3Range1K500DdepF_birth2extrMort01lin"





nameFileVec = "Run08b2MonthRangeYearAlphaMax1Rate5_ratioRev10000_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
divisionYear = 12
divisionTime = 1000


nameFileVec = "Run08b2MonthRangeYearAlphaMax1Rate5_ratioRev1_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
divisionYear = 12
divisionTime = 1000
# 
nameFileVec = "Run08b2MonthRangeMonthAlphaMax1Rate5_ratioRev1_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort01lin"
divisionYear = 12
divisionTime = 1000

    
nameFileVec = "Run10DayRangeDayAlphaMax10Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin"
divisionYear = 365
divisionTime = 1000
# 
nameFileVec = "Run10MonthRangeMonthAlphaMax50Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut1e3Range1K500DdepF_birth2extrMort04lin"
divisionYear = 12
divisionTime = 1000

nameFileVec = "Run11bis2MonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth03extrMort01lin"
nameFileVec = "Run11bis2MonthRangeMonthAlphaMax1Rate50_ratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth03extrMort01lin"
nameFileVec = "Run11bis2MonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# nameFileVec = "Run11bis2MonthRangeMonthAlphaMax1Rate50_ratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# nameFileVec = "Run11bis2MonthRangeYearAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth03extrMort01lin"
  
# nameFileVec ="Run11bis2MonthRangeMonthAlphaMax10Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# nameFileVec ="Run11bis2MonthRangeMonthAlphaMax1Rate50_ratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"

# nameFileVec = "RunSmallMut01MonthRangeMonthAlphaMax1Rate50_ratioRev1_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# nameFileVec = "RunSmallMut01MonthRangeMonthAlphaMax1Rate50_ratioRev10_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# nameFileVec = "RunSmallMut01MonthRangeMonthAlphaMax1Rate50_ratioRev100_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"
# 
# nameFileVec = "RunSmallMut01MonthRangeYearAlphaMax1Rate50_ratioRev1_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth05extrMort01lin"


nameFileVec = "Run12MonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin"

nameFileVec = "Run12MonthRangeYearAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin"

nameFileVec = "Run12somaticMonthRangeYearAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin"
nameFileVec = "Run12somaticMonthRangeMonthAlphaMax1Rate50_NoratioRev_Rep1_maxAge50QuantStart09_Mut2e3Range1K500DdepF_birth15extrMort1lin"


divisionYear = 12
divisionTime = 1000
  
# divisionYear = 12
# divisionTime = 1000
# divisionYear = 12
# divisionTime = 1000 

#######################

firstFile = 1
for (i in 1:length(nameFileVec)){
  dataPopOne <- read.csv(paste("Data/dataPop_",nameFileVec[i],".csv",sep=""), sep=";")
  dataMutOne <- read.csv(paste("Data/dataMut_",nameFileVec[i],".csv",sep=""), sep=";")
  if(firstFile==1){
    dataPop = dataPopOne
    dataMut = dataMutOne
    firstFile = 0  
  }else{
    dataPop = rbind(dataPop,dataPopOne)
    dataMut = rbind(dataMut,dataMutOne)
  }
}
dataPop$propDeathMutation2 = dataPop$propDeathMutation/(dataPop$propDeathMutation+dataPop$propDeathExtrinsic)
summary(dataPop)


dataMut = dataMut[dataMut$Time %in% 
                    c(min(dataMut$Time),unique(dataMut$Time)[round(length(unique(dataMut$Time))/2)],max(dataMut$Time))   ,]

summary(dataMut)
dataMutOne = dataMut[1,]
dataMutOne[1,4:7] = 0
for(rep in unique(dataMut$Rep)){
  for(time in unique(dataMut$Time)){
    for(damage in c(unique(dataMut$Damage),max(dataMut$Damage)+1)){
        if(length(dataMut$MeanSurvivalMutation[dataMut$Rep==rep & dataMut$Time==time & dataMut$Damage==damage])==0){
        dataMutOne$Rep = rep
        dataMutOne$Time = time
        dataMutOne$Damage = damage
        dataMut <- rbind(dataMut,
                         dataMutOne)
      }
    }
  }
}

plotMutSurvMean <- ggplot(dataMut)+
  geom_line( aes(Damage/divisionYear,1-MeanSurvivalMutation,color=as.factor(round(Time/divisionYear/divisionTime))))+
  geom_point( aes(Damage/divisionYear,1-MeanSurvivalMutation,color=as.factor(round(Time/divisionYear/divisionTime))))+
  # geom_hline(yintercept = 0.5)+
  scale_color_discrete(name=expression(paste("Time (x",10^3," years)",sep="")))+
  theme(legend.position = c(0.8, 0.4))+
  ylab("Proportion of individuals carrying\nlethal mutations expressed at each age")+
  xlab("Age at which mutations are expressed (years)")+
  ylim(0,1)+
  xlim(0,max(dataPop$LifeSpan0975/divisionYear,na.rm=T))

plotMutSurvMed <- ggplot(dataMut)+
  geom_line( aes(Damage/divisionYear,1-SurvivalMutation05,color=as.factor(round(Time/divisionYear/divisionTime))))+
  geom_point( aes(Damage/divisionYear,1-SurvivalMutation05,color=as.factor(round(Time/divisionYear/divisionTime))))+
  scale_color_discrete(name=expression(paste("Time (x",10^3," years)",sep="")))+
  theme(legend.position = c(0.8, 0.4))+
  ylab("Probability of dying\nbecause of mutations (med)")+
  xlab("Age (years)")+
  ylim(0,1)+
  xlim(0,max(dataPop$LifeSpan0975/divisionYear,na.rm=T))

plotLifeSpan <- ggplot(dataPop)+
  geom_line( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  geom_point( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  # geom_ribbon(aes(Time/divisionYear,ymin=LifeSpan0025/divisionYear, ymax=LifeSpan0975/divisionYear), alpha=0.1)+
  ylim(0,1.2*max(dataPop$LifeSpan05/divisionYear))+
  ylab("Median life span (years)")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))

dataPop$propDeathMutation3 = dataPop$propDeathMutation*dataPop$propDeath
dataPop$propDeathExtrinsic3 = dataPop$propDeathExtrinsic*dataPop$propDeath
dataPop$propDeathYear = 0
dataPop$propDeathMutationYear = 0
# dataPop$propDeathExtrinsicYear = 0
for(i in 1:divisionYear){
  dataPop$propDeathYear = dataPop$propDeathYear + dataPop$propDeath*(1-dataPop$propDeath)^(i-1)
  dataPop$propDeathMutationYear = dataPop$propDeathMutationYear + dataPop$propDeathMutation3*(1-dataPop$propDeath)^(i-1)
  # dataPop$propDeathExtrinsicYear = dataPop$propDeathExtrinsicYear + dataPop$propDeathExtrinsic3*(1-dataPop$propDeath)^(i-1)
}

dataPop2 <- dataPop[,c("Time","propDeathYear","propDeathMutationYear")]
dataPop2$Up = dataPop2$propDeathYear
dataPop2$Bottom = dataPop2$propDeathMutationYear
dataPop2$cause = "Extrinsic mortality"
dataPop2b <- dataPop[,c("Time","propDeathYear","propDeathMutationYear")]
dataPop2b$Up = dataPop2$propDeathMutationYear
dataPop2b$Bottom = 0
dataPop2b$cause = "Intrinsic mortality"
dataPop2 <- rbind(dataPop2,dataPop2b)
dataPop2$cause <- as.factor(dataPop2$cause)

plotCausesDeath <- ggplot(dataPop2)+
  geom_ribbon(aes(Time/divisionYear/divisionTime, ymin=Bottom, ymax=Up, fill=cause))+
  scale_fill_manual(name="Cause of death",values = c("#bdbdbd","#494949"))+
  geom_line(data=dataPop, aes(Time/divisionYear/divisionTime,propDeathYear))+
  geom_point(data=dataPop, aes(Time/divisionYear/divisionTime,propDeathYear))+
  ylab("Probability of dying each year")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  ylim(0,1)+
  theme(legend.position = c(0.05, 0.9))


dataPop2 <- dataPop[,c("Time","propDeath","propDeathMutation3")]
dataPop2$Up = dataPop2$propDeath
dataPop2$Bottom = dataPop2$propDeathMutation3
dataPop2$cause = "Extrinsic mortality"
dataPop2b <- dataPop[,c("Time","propDeath","propDeathMutation3")]
dataPop2b$Up = dataPop2$propDeathMutation3
dataPop2b$Bottom = 0
dataPop2b$cause = "Intrinsic mortality"
dataPop2 <- rbind(dataPop2,dataPop2b)
dataPop2$cause <- as.factor(dataPop2$cause)

plotCausesDeathMonth <- ggplot(dataPop2)+
  geom_ribbon(aes(Time/divisionYear/divisionTime, ymin=Bottom, ymax=Up, fill=cause))+
  scale_fill_manual(name="Cause of death",values = c("#bdbdbd","#494949"))+
  geom_line(data=dataPop, aes(Time/divisionYear/divisionTime,propDeath))+
  geom_point(data=dataPop, aes(Time/divisionYear/divisionTime,propDeath))+
  ylab("Proportion of individuals\ndying each month")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  ylim(0,0.05)+
  theme(legend.position = c(0.07, 0.8))


dataPop2 <- dataPop[,c("Time","propDeath","propDeathExtrinsic3")]
dataPop2$Up = dataPop2$propDeath
dataPop2$Bottom = dataPop2$propDeathExtrinsic3
dataPop2$cause = "Intrinsic mortality"
dataPop2b <- dataPop[,c("Time","propDeath","propDeathExtrinsic3")]
dataPop2b$Up = dataPop2$propDeathExtrinsic3
dataPop2b$Bottom = 0
dataPop2b$cause = "Extrinsic mortality"
dataPop2 <- rbind(dataPop2,dataPop2b)
dataPop2$cause <- factor(dataPop2$cause,levels=c("Intrinsic mortality","Extrinsic mortality"))

plotCausesDeathMonth <- ggplot(dataPop2)+
  geom_ribbon(aes(Time/divisionYear/divisionTime, ymin=Bottom, ymax=Up, fill=cause))+
  scale_fill_manual(name="Cause of death",values = c("#616161","#C7C7C7"))+
  geom_line(data=dataPop, aes(Time/divisionYear/divisionTime,propDeath))+
  geom_point(data=dataPop, aes(Time/divisionYear/divisionTime,propDeath))+
  ylab("Proportion of individuals\ndying each month")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  ylim(0,max(0.05,max(dataPop2$Up)*1.3))+
  theme(legend.position = c(0.07, 0.8))

plotLifeSpanRibbon <- ggplot(dataPop)+
  geom_line( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  # geom_point( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  geom_ribbon(aes(Time/divisionYear/divisionTime,ymin=LifeSpan0025/divisionYear, ymax=LifeSpan0975/divisionYear), alpha=0.1)+
  ylab("Life span (years)")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  theme(text = element_text(size=12),
        axis.text =  element_text(size=10))

plotLifeSpanRibbon2 <- ggplot(dataPop)+
  geom_line( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  geom_point( aes(Time/divisionYear/divisionTime,LifeSpan05/divisionYear))+
  geom_ribbon(aes(Time/divisionYear/divisionTime,ymin=LifeSpan0025/divisionYear, ymax=LifeSpan0975/divisionYear), alpha=0.1)+
  geom_line( aes(Time/divisionYear/divisionTime,LifeSpan0025/divisionYear),lty=2, alpha=0.5)+
  geom_line( aes(Time/divisionYear/divisionTime,LifeSpan0975/divisionYear),lty=2, alpha=0.5)+
  ylab("Life span (years)")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))

# Specify position of plot2 (in percentages of plot1)
g2 <- ggplotGrob(plotLifeSpanRibbon)
# minLifeSpanPlot =  min(dataPop$LifeSpan05/divisionYear,na.rm=T) 
minLifeSpanPlot =  0
maxLifeSpanPlot =  1.2 * max(dataPop$LifeSpan05/divisionYear,na.rm=T)
ymax  = maxLifeSpanPlot
ymin  = minLifeSpanPlot + (maxLifeSpanPlot-minLifeSpanPlot)/5
xmin  = max(dataPop$Time/divisionYear/divisionTime,na.rm=T)
xmax  = min(dataPop$Time/divisionYear/divisionTime,na.rm=T) + (max(dataPop$Time/divisionYear/divisionTime,na.rm=T)-min(dataPop$Time/divisionYear/divisionTime,na.rm=T))/1.7
plotLifeSpanInset = plotLifeSpan + annotation_custom(grob = g2, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

# ou bien
# plotLifeSpanInset = plotLifeSpan

plotOffspringDying <- ggplot(dataPop)+
  geom_line( aes(Time/divisionYear/divisionTime,PercentOffspringDead))+
  geom_point( aes(Time/divisionYear/divisionTime,PercentOffspringDead))+
  # geom_hline(yintercept = 0, lty=2)+
  # geom_hline(yintercept = 1, lty=2)+
  ylab("Proportion of offspring dying\nbecause of local competition")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  ylim(c(0,1))

plotPopSize <- ggplot(dataPop)+
  geom_line( aes(Time/divisionYear/divisionTime,SizePop))+
  geom_point( aes(Time/divisionYear/divisionTime,SizePop))+
  # geom_hline(yintercept = max(dataPop$SizePop,na.rm=TRUE), lty=2)+
  # geom_hline(yintercept = 0.95*max(dataPop$SizePop,na.rm=TRUE), lty=2)+
  ylab("Population size")+
  xlab(expression(paste("Time (x",10^3," years)",sep="")))+
  ylim(c(0,NA))

  
#### adjust margin


PlotMutSurvMean     <- ggplotGrob(plotMutSurvMean)
PlotMutSurvMed      <- ggplotGrob(plotMutSurvMed)
PlotLifeSpan        <- ggplotGrob(plotLifeSpanRibbon2)
PlotLifeSpanInset   <- ggplotGrob(plotLifeSpanInset)
PlotCausesDeath     <- ggplotGrob(plotCausesDeath)
PlotCausesDeathMonth  <- ggplotGrob(plotCausesDeathMonth)
PlotOffspringDying  <- ggplotGrob(plotOffspringDying)
PlotPopSize         <- ggplotGrob(plotPopSize)
maxWidth = grid::unit.pmax(PlotMutSurvMean$widths[2:5],
                           PlotMutSurvMed$widths[2:5],
                           PlotLifeSpan$widths[2:5],
                           PlotLifeSpanInset$widths[2:5],
                           PlotCausesDeath$widths[2:5],
                           PlotCausesDeathMonth$widths[2:5],
                           PlotOffspringDying$widths[2:5],
                           PlotPopSize$widths[2:5])
PlotMutSurvMean$widths[2:5] <- as.list(maxWidth)
PlotMutSurvMed$widths[2:5] <- as.list(maxWidth)
PlotLifeSpan$widths[2:5] <- as.list(maxWidth)
PlotLifeSpanInset$widths[2:5] <- as.list(maxWidth)
PlotCausesDeath$widths[2:5] <- as.list(maxWidth)
PlotCausesDeathMonth$widths[2:5] <- as.list(maxWidth)
PlotOffspringDying$widths[2:5] <- as.list(maxWidth)
PlotPopSize$widths[2:5] <- as.list(maxWidth)

plotAll <- 
  plot_grid(
    PlotMutSurvMean,
    plot_grid(
      PlotLifeSpan,   #PlotLifeSpanInset
      PlotPopSize,
      PlotCausesDeathMonth,
      PlotOffspringDying,
      labels = c("b","c","d","e"),
      rel_heights = c(1, 1),
      ncol = 2
    ),
    labels = c("a",NA),
    rel_heights = c(1,2),
    ncol = 1
  )

save_plot("Graphes/plotSimulation.pdf",
          plotAll,
          base_height = 10, base_aspect_ratio = 1)
plotAll

