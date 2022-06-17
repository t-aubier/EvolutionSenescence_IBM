library(scales)



ageVec = seq(0,100,by=1)



rateAccumulExp = 0.3


extrinsicMort = 0.5


Nind = 10000

maxDamagePlot = 300


damageVecLin = ageVec
damageVecExp = floor(exp(ageVec * rateAccumulExp)-1);

dataFunctions <- data.frame(age=ageVec,damage=damageVecLin,func="lin")
dataFunctions <- rbind(dataFunctions,data.frame(age=ageVec,damage=damageVecExp,func="exp"))
dataFunctions$func <- factor(dataFunctions$func)

probSurvExtrinsicMortality = exp(-extrinsicMort/12.0)
damagesIndLIN <-
  base::replicate(Nind, {
    func = "lin"
    damageVec <- c()
    dead = FALSE
    age = 0 
    while(dead==FALSE){
      if(func=="exp"){
        damage = floor(exp(age * rateAccumulExp)-1)
      }else if(func=="lin"){
        damage = age
      }
      damageVec <-c(damageVec,damage)
      
      if(runif(1)>probSurvExtrinsicMortality){
        dead = TRUE
      }
      age = age+1;
    }
    damageVec
  })
damagesIndLIN <- unlist(damagesIndLIN)
medLin <- median(damagesIndLIN)
maxDamageLin <- quantile(damagesIndLIN,probs=0.975)[[1]]


damagesIndEXP <-
  base::replicate(Nind, {
  func = "exp"
  damageVec <- c()
  dead = FALSE
  age = 0 
  while(dead==FALSE){
    if(func=="exp"){
      damage = floor(exp(age * rateAccumulExp)-1)
    }else if(func=="lin"){
      damage = age
    }
    damageVec <-c(damageVec,damage)
    if(runif(1)>probSurvExtrinsicMortality){
      dead = TRUE
    }
    age = age+1;
  }
  damageVec
})
damagesIndEXP <- unlist(damagesIndEXP)
medExp <- median(damagesIndEXP)
maxDamageExp <- quantile(damagesIndEXP,probs=0.975)[[1]]

damagesInd <- rbind(data.frame(damage=damagesIndLIN,fun="lin"),data.frame(damage=damagesIndEXP,fun="exp"))

sizeTextAxis = 15.1290493

plotFunctions <- 
  ggplot(dataFunctions)+
  geom_line(aes(age,damage,color=func),size=1.5)+
  scale_color_manual(name="",values=c("orange","red"))+
  theme(legend.position="none",
        plot.margin = margin(0, -0.2, 0, 0, "cm"),
        axis.text = element_text(size=sizeTextAxis))+
  scale_y_continuous(limits=c(0,maxDamagePlot))+
  xlab("Age")+
  ylab("Damage")


plotHist1 <- 
  ggplot(data=damagesInd[damagesInd$fun=="lin",],aes(damage))+
  geom_histogram(aes(y = (..count..)/sum(..count..)),color="black", fill="orange", alpha=1, position = 'identity',bins=20)+ 
  geom_vline(xintercept = medLin,color="orange",size=1.2,lty=2)+
  scale_y_continuous(limits=c(0,1),labels = percent,breaks = c(0,1))+
  scale_x_continuous(limits=c(NA,maxDamagePlot))+
  theme(legend.position="none",
        plot.margin = margin(0, -0., 0, -0., "cm"),
        axis.text = element_text(size=sizeTextAxis))+
  coord_flip()+
  xlab("")+
  ylab("")

plotHist2 <- 
  ggplot(data=damagesInd[damagesInd$fun=="exp",],aes(damage))+
  geom_histogram(aes(y = (..count..)/sum(..count..)),color="black", fill="red", alpha=1, position = 'identity',bins=20)+ 
  geom_vline(xintercept = medExp,color="red",size=1.2,lty=2)+
  scale_y_continuous(limits=c(0,1),labels = percent,breaks = c(0,1))+
  scale_x_continuous(limits=c(NA,maxDamagePlot))+
  theme(legend.position="none",
        plot.margin = margin(0, -0., 0, -0., "cm"),
        axis.text = element_text(size=sizeTextAxis))+ 
  coord_flip()+
  xlab("")+
  ylab("")

plotAll <- plot_grid(plotFunctions,plotHist1,plotHist2, rel_widths = c(3.5,1,1),ncol=3)

save_plot("Graphes/plotAccumulFunctions.pdf",
          plotAll,
          base_height = 4, base_aspect_ratio = 2.5)

plotAll
print(medLin)
print(medExp)


