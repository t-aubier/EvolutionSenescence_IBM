library(scales)


convertIntoYears = 12
ageVec = seq(0,30*convertIntoYears,by=1)
gammaP =  50   /365
alphamaxVec = c(1,5,50)

maxFecundityPlot = 50

AlphamaxVec = c()
AgeVec = c()
GammaVec = c()
for(alphamax in alphamaxVec){
  GammaOutput = 1+(alphamax-1)*exp(-gammaP/convertIntoYears* ageVec)
  
  AlphamaxVec = c(AlphamaxVec,rep(alphamax,length(ageVec)))
  AgeVec = c(AgeVec,ageVec)
  GammaVec = c(GammaVec,GammaOutput)
  
}

dataFunctions <- data.frame(age=AgeVec,gamma=GammaVec,alphamax=AlphamaxVec)

sizeTextAxis = 12

plotFunctions <- 
  ggplot(dataFunctions)+
  geom_line(aes(age/convertIntoYears,gamma,color=as.factor(alphamax)),size=1.5)+
  scale_color_manual(name=expression(alpha[max]),
                     values=c("orange","red","purple"))+
  theme(legend.position="right",
        plot.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
        axis.text = element_text(size=sizeTextAxis))+
  scale_x_continuous(breaks=c(0,15,30))+
  scale_y_continuous(limits=c(0,maxFecundityPlot),breaks=c(1,5,50))+
  xlab(expression("Age at which intrinsic mortality occurs ("~a[die]~", converted in years)"))+
  ylab(expression("Factor change in fecundity ("~Gamma~")"))



save_plot("Graphes/plotPleiotropyFunction.pdf",
          plotFunctions,
          base_height = 4.8, base_aspect_ratio = 1.4)

plotFunctions


