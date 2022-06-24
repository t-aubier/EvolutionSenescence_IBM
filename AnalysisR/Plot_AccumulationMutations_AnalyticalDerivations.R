##### Directory

setwd("/home/taubier/Projets/Github/EvolutionSenescence_IBM")     # change the working directory


##### Packages

library(scales)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())


##########  PLOT THE ACCUMULATION OF LETHAL MUTATIONS AND THE MAXIMUM LIFESPAN FOR DIFFERENT GRAINS OF AGE DEPENDENCE OF MUTATION EXPRESSION  ##########  


##### Parameters

nbMut = 1000                    # Number of mutations that accumulate
Ne0 = 1000                      # Baseline effective population size
alpha0 = 1                      # Baseline maximum fecundity when pleiotropy (here no pleiotropy)
mu0 = 0.2                       # Baseline extrinsic mortality rate
maxYPlot = 100                  # Maximum value on the vertical axis

deltaVec = seq(0,30,by=0.2)     # Grain of age dependence tested



##### Process data and plot for variation in effective population size

paramVec = c(100,1000,10000)    # Effective population size tested
Ne = Ne0
alpha = alpha0
mu = mu0
vecMut = c()
vecParam = c()
vecX = c()
for(param in paramVec){
  Ne = param
  X=1e19
  for(i in 1:nbMut){
    x=1/mu *log(alpha/(Ne*(alpha-1)+1+exp(-mu*X)*(Ne-1))*Ne)
    X = x
    if(i<100 || ( i<500 && i%%5==0)|| ( i<1000 && i%%15==0) || ( i<2000 && i%%30==0) ||( i<10000 && i%%150==0)){
      vecMut = c(vecMut,i)
      vecParam = c(vecParam,param)
      vecX = c(vecX,x)
    }
  }
}
dataNe=data.frame(vecMut,vecX,vecParam)
plotNe = ggplot(dataNe)+
  geom_point(aes(vecMut,vecX, color = as.factor(vecParam)))+
  scale_x_continuous(trans='log10')+
  scale_colour_grey(start = 0.8,end = 0.2)+
  xlab("Sequence of mutations fixed")+
  ylab("Maximum age at death")+
  scale_y_continuous(limits=c(0,maxYPlot),breaks=c(0,50,100))


Ne = Ne0
alpha = alpha0
mu = mu0
vecDelta = c()
vecParam = c()
vecX = c()
vecXmax = c()
for(param in paramVec){
  Ne = param
  for(delta in deltaVec){
    
    
    ##x=1/mu *log(alpha/(Ne*(alpha-1)+1+exp(-mu*X)*(Ne-1))*Ne)
    x=1/mu *log(
      (alpha * Ne - exp(-mu*delta) * (Ne - 1)) / (Ne * (alpha-1)+1)
    )
    xmax=1/mu *log(
      (Ne*(alpha*exp(mu*delta)-1)+1) / (Ne * (alpha-1)+1)
    )   
    
    vecDelta = c(vecDelta,delta)
    vecX = c(vecX,x)
    vecXmax = c(vecXmax,xmax)
    vecParam = c(vecParam,param)
  }
}
dataDeltaNe=data.frame(vecX,vecXmax,vecDelta,vecParam)
plotDeltaNe = ggplot(dataDeltaNe)+
  geom_line(aes(vecDelta,vecX,color=as.factor(vecParam)),size=1.5)+
  geom_line(aes(vecDelta,vecXmax,color=as.factor(vecParam)),size=1.5)+
  geom_ribbon(aes(x=vecDelta,ymin=vecX,ymax=vecXmax,fill=as.factor(vecParam)), alpha=0.5) +
  scale_colour_grey(start = 0.8,end = 0.2)+
  scale_fill_grey(start = 0.8,end = 0.2)+
  xlab("Grain of age dependence")+
  ylab("Maximum age at death at evolutionary equilibrium")+
  scale_y_continuous(limits=c(-0.0001,maxYPlot),breaks=c(0,50,100))



##### Process data and plot for variation in maximum fecundity when pleiotropy

paramVec = c(0.1,0.2,0.4)       # Maximum fecundity tested
Ne = Ne0
alpha = alpha0
mu = mu0
vecMut = c()
vecParam = c()
vecX = c()
for(param in paramVec){
  mu = param
  X=1e19
  for(i in 1:nbMut){
    x=1/mu *log(alpha/(Ne*(alpha-1)+1+exp(-mu*X)*(Ne-1))*Ne)
    X = x
    if(i<100 || ( i<500 && i%%5==0)|| ( i<1000 && i%%15==0) || ( i<2000 && i%%30==0) ||( i<10000 && i%%150==0)){
      vecMut = c(vecMut,i)
      vecParam = c(vecParam,param)
      vecX = c(vecX,x)
    }
  }
}
dataMu=data.frame(vecMut,vecX,vecParam)
plotMu = ggplot(dataMu)+
  geom_point(aes(vecMut,vecX, color = as.factor(vecParam)))+
  scale_x_continuous(trans='log10')+
  scale_colour_grey(start = 0.8,end = 0.2)+
  xlab("Sequence of mutations fixed")+
  ylab("Maximum age at death")+
  scale_y_continuous(limits=c(0,maxYPlot),breaks=c(0,50,100))


Ne = Ne0
alpha = alpha0
mu = mu0
vecDelta = c()
vecParam = c()
vecX = c()
vecXmax = c()
for(param in paramVec){
  mu = param
  for(delta in deltaVec){
    x=1/mu *log(
      (alpha * Ne - exp(-mu*delta) * (Ne - 1)) / (Ne * (alpha-1)+1)
    )
    xmax=1/mu *log(
      (Ne*(alpha*exp(mu*delta)-1)+1) / (Ne * (alpha-1)+1)
    )   
    vecDelta = c(vecDelta,delta)
    vecX = c(vecX,x)
    vecXmax = c(vecXmax,xmax)
    vecParam = c(vecParam,param)
  }
}
dataDeltaMu=data.frame(vecX,vecXmax,vecDelta,vecParam)
plotDeltaMu = ggplot(dataDeltaMu)+
  geom_line(aes(vecDelta,vecX,color=as.factor(vecParam)),size=1.5)+
  geom_line(aes(vecDelta,vecXmax,color=as.factor(vecParam)),size=1.5)+
  geom_ribbon(aes(x=vecDelta,ymin=vecX,ymax=vecXmax,fill=as.factor(vecParam)), alpha=0.5) +
  scale_colour_grey(start = 0.8,end = 0.2)+
  scale_fill_grey(start = 0.8,end = 0.2)+
  xlab("Grain of age dependence")+
  ylab("Maximum age at death at evolutionary equilibrium")+
  scale_y_continuous(limits=c(-0.0001,maxYPlot),breaks=c(0,50,100))



##### Process data and plot for variation in extrinsic mortality rate

paramVec = c(1,1.05,1.1)        # Extrinsic mortality rate tested
Ne = Ne0
alpha = alpha0
mu = mu0
vecMut = c()
vecParam = c()
vecX = c()
for(param in paramVec){
  alpha = param
  X=1e19
  for(i in 1:nbMut){
    x=1/mu *log(alpha/(Ne*(alpha-1)+1+exp(-mu*X)*(Ne-1))*Ne)
    X = x
    
    if(i<100 || ( i<500 && i%%5==0)|| ( i<1000 && i%%15==0) || ( i<2000 && i%%30==0) ||( i<10000 && i%%150==0)){
      vecMut = c(vecMut,i)
      vecParam = c(vecParam,param)
      vecX = c(vecX,x)
    }
  }
}
dataAlpha=data.frame(vecMut,vecX,vecParam)
plotAlpha = ggplot(dataAlpha)+
  geom_point(aes(vecMut,vecX, color = as.factor(vecParam)))+
  scale_x_continuous(trans='log10')+
  scale_colour_grey(start = 0.8,end = 0.2)+
  xlab("Sequence of mutations fixed")+
  ylab("Maximum age at death")+
  scale_y_continuous(limits=c(0,maxYPlot),breaks=c(0,50,100))


Ne = Ne0
alpha = alpha0
mu = mu0
vecDelta = c()
vecParam = c()
vecX = c()
vecXmax = c()
for(param in paramVec){
  alpha = param
  for(delta in deltaVec){
    x=1/mu *log(
      (alpha * Ne - exp(-mu*delta) * (Ne - 1)) / (Ne * (alpha-1)+1)
    )
    xmax=1/mu *log(
      (Ne*(alpha*exp(mu*delta)-1)+1) / (Ne * (alpha-1)+1)
    )    
    
    vecDelta = c(vecDelta,delta)
    vecX = c(vecX,x)
    vecXmax = c(vecXmax,xmax)
    vecParam = c(vecParam,param)
  }
}
dataDeltaAlpha=data.frame(vecX,vecXmax,vecDelta,vecParam)
plotDeltaAlpha = ggplot(dataDeltaAlpha)+
  geom_line(aes(vecDelta,vecX,color=as.factor(vecParam)),size=1.5)+
  geom_line(aes(vecDelta,vecXmax,color=as.factor(vecParam)),size=1.5)+
  geom_ribbon(aes(x=vecDelta,ymin=vecX,ymax=vecXmax,fill=as.factor(vecParam)), alpha=0.5) +
  scale_colour_grey(start = 0.8,end = 0.2)+
  scale_fill_grey(start = 0.8,end = 0.2)+
  xlab("Grain of age dependence")+
  ylab("Maximum age at death at evolutionary equilibrium")+
  scale_y_continuous(limits=c(-0.0001,maxYPlot),breaks=c(0,50,100))



##### Compilation of all plots

plotAll <- plot_grid(
  plotNe+theme(legend.position="none"),
  plotMu+theme(legend.position="none"),
  plotAlpha+theme(legend.position="none")
)
plotAll

plotAllDelta <- plot_grid(
  plotDeltaNe+theme(legend.position="none"),
  plotDeltaMu+theme(legend.position="none"),
  plotDeltaAlpha+theme(legend.position="none")
)
plotAllDelta

save_plot("Graphes/Figure_AccumulationMutations_AnalyticalResults1.pdf",
          plotAll,
          base_height = 8, base_aspect_ratio = 1)

save_plot("Graphes/Figure_AccumulationMutations_AnalyticalResults2.pdf",
          plotAllDelta,
          base_height = 8, base_aspect_ratio = 1)

