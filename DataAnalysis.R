## ----knitr_option, cache=FALSE, echo=FALSE, results='hide'---------------
library(knitr)

## set global chunk options
opts_chunk$set(echo=FALSE, fig.align='center',fig.show="asis", out.width="0.95\\textwidth",  par=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=55), warning=FALSE, error=FALSE)


## ----R_lab,warning=FALSE, results='hide', echo=FALSE---------------------
#####
# Libraries 
#####
library("RColorBrewer")
library(plotrix)
library(graphics)
library(ggplot2)
library(MGLM)
library(xtable)
library(mvtnorm)
library(combinat)

## ----Data----------------------------------------------------------------
source('Data i.o.R')

## ----Multinom_reg, echo=1:2, results = 'asis', out.width=6---------------
Multinom_mod <- MGLMreg(as.matrix(RotaVirusBB[,agNames])~t+sin+cos,
        data=RotaVirusBB, dist="MN") 

test <- do.call("cbind", list(Multinom_mod$coeff,Multinom_mod$test))
xtable(test, align=c("l|cccccc"), digits=3, caption="Multinomial regression parameter estimate and Wald tests for each parameter.", label="tab:MultinomTab")

## ----Oddsdevelop, echo=1:3, fig.cap=c("Odds for age category ''00-04'' plotted over time, grey area indicates values of odds less than one"), fig.pos="H", fig.height=4----
beta_hat <- Multinom_mod$coeff
Design_mat <- t(Multinom_mod$data$X)
Odds <- exp(beta_hat[,1]%*%Design_mat) # choose category "00-04"
linearOdds <- exp(c(beta_hat[1,1],beta_hat[2,1])%*%matrix(c(Design_mat[1,],Design_mat[2,]),ncol=144,byrow=TRUE))

LTO <- which(Odds < 1) # which odds are less than one 
v<- which(RotaVirusBB$time=="2009-01-01") 

#####
#Plot
#####
plot(RotaVirusBB$time,Odds, type="l", lty=2, lwd=1, ylab="Odds", xlab="Time (Months)")
lines(RotaVirusBB$time, linearOdds, type="l", lty=1, lwd=2, ylab="Odds", xlab="Time (Months)")
axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
abline(v=RotaVirusBB$time[v],col="black",lty=5, lwd=2)
rect(RotaVirusBB$time[min(LTO)],0,RotaVirusBB$time[max(LTO)],max(Odds), border=NA, col=adjustcolor("grey60", alpha.f = 0.2))
title("Multinomial logit regression, Development of odds over time")
legend("bottomleft",legend=c("Odds without harmonic covariates", "Odds with harmonic covariates"), lty=c(1,2), lwd=c(2,1), pt.cex=1,cex=0.8)

## ----Dirichletfit, echo=1:2, results='asis'------------------------------
Dirichlet_mod <- MGLMreg(as.matrix(RotaVirusBB[,agNames])~t+sin+cos,
        data=RotaVirusBB, dist="DM") 

test2 <- do.call("cbind", list(Dirichlet_mod$coeff,Dirichlet_mod$test))
xtable(test2, 
       caption="Estimated coefficients of Dirichlet-multinomial model and Wald tests for regression parameters", align=c("l|ccccccc"), digits=3,label="tab:DirMultinomTab")

## ----PlotMultinom_mod, echo=1:2, fig.pos='H', fig.cap=c("Fitted mean values age-categories ''00-04'', ''05-09'' ''10-14'', ''15-69'' and  ''70+''"), fig.height=7----
pred_Multinom <- predict(Multinom_mod, newdata=t(Design_mat))
pred_dirich <- predict(Dirichlet_mod, newdata=t(Design_mat))

grey <- adjustcolor("grey20", alpha.f = 0.9)
Titles <- c("Fitted values for category ''00-04''","Fitted values for category ''05-09''","Fitted values for category ''10-14''","Fitted values for category ''15-69''","Fitted values for category ''70+''")
Place <- c("bottomleft", NA, "topright", NA, "topleft")
plotFit <- function() {
  par(mfrow=c(2,1))
    for (i in 1:(length((RotaVirusBB[,agNames])))) {
        plot(RotaVirusBB$time, pred_Multinom[,i], type="l", ylim=c(0,max(RotaVirusBB[,agNames[i]]
              /RotaVirusBB[,"totalAGs"])), lwd=2, lty=1,xlab="Time (months)", ylab="")
         lines(RotaVirusBB$time, pred_dirich[,i], type="l", ylim=c(0,max(RotaVirusBB[,agNames[i]]/RotaVirusBB[,"totalAGs"])), lwd=2, lty=3,xlab="Time (months)", ylab="")
        points(RotaVirusBB$time,RotaVirusBB[,agNames[i]]/RotaVirusBB[,"totalAGs"],lwd=1,
             col=grey, pch=1)
  axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
  title(Titles[i],cex.main=0.8 ,ylab="Proportions of RotaVirus cases and models fit", cex.lab=0.7)
     if ((i==1|i==3|i==5)) {
  legend(Place[i],legend=c("Multinomial logit regression", "Dirichlet-multinomial regression", "Propotions of reported rotavirus cases"), lty=c(1,3,NA), pch=c(NA,NA,1), border="black",pt.cex=1,cex=0.77, lwd=c(1,2,NA)) }
   }
}
plotFit()

## ----AICBIC, results='asis'----------------------------------------------
AIC_BIC<- data.frame(x=matrix(c(Multinom_mod$AIC,Multinom_mod$BIC),ncol=1),y=matrix(c( Dirichlet_mod$AIC,Dirichlet_mod$BIC),ncol=1), row.names=c("AIC","BIC"))
colnames(AIC_BIC) <- c("Multinomial regression", "Dirichlet-multinomial regression") 
xtable(AIC_BIC, caption="Goodness-of-fit measures for the two regression models representing a tradeoff between fit and model complexity", 
       align=c("l|cc"), label="AICBIC")

## ----PredicitveInterval0004, echo=2:14, fig.cap="95 percent Predicitive intervals for multinomial and Dirichlet-multinomial, age category ''00-04''", fig.pos="H",fig.lp="fig:", fig.height=4.9----
source('~/Desktop/Bachelor Multinomial time series modelling/ReportLayout 06052014/DataAnalysis-DirMultinom-timeseries/MGLM-Sampling-fun.R') # call: MGLMpredInt()

N <- 1000 # Number of replicates. 

New_propMultinom <- replicate(N,MGLMpredInt(Multinom_mod)) 
New_propDirichlet <- replicate(N,MGLMpredInt(Dirichlet_mod))


oneAGPI <- function(ageIndex, predictive) {
 res <- matrix(apply(X= predictive[,ageIndex,],1, FUN=quantile, probs=c(0.025,0.975)),ncol=2, byrow=TRUE)
 return(res)
}

ageInd <- c("00-04"=1,"05-09"=2,"10-14"=3, "15-69"=4,"70+"=5)
piMulti <- lapply(ageInd, oneAGPI, predictive=New_propMultinom)
piDirMulti <- lapply(ageInd, FUN=oneAGPI, predictive=New_propDirichlet)

# Category 00-04 
plot(RotaVirusBB$time, pred_Multinom[,1], type="l", ylim=c(0,1), lwd=1, lty=1,xlab="Time (months)", ylab="")
 lines(RotaVirusBB$time, pred_dirich[,1], type="l", ylim=c(0,max(RotaVirusBB[,agNames[1]]/RotaVirusBB[,"totalAGs"])), lwd=1, lty=3,xlab="Time (months)", ylab="")
points(RotaVirusBB$time,RotaVirusBB[,"00-04"]/RotaVirusBB[,"totalAGs"],lwd=1,col=grey, pch=1)

lines(RotaVirusBB$time, piMulti$`00-04`[,1], lwd=2)
lines(RotaVirusBB$time, piMulti$`00-04`[,2], lwd=2)
polygon(c(RotaVirusBB$time,rev(RotaVirusBB$time)),c(piMulti$`00-04`[,2],rev(piMulti$`00-04`[,1])), col=adjustcolor("grey20", alpha.f = 0.3))

lines(RotaVirusBB$time, piDirMulti$`00-04`[,1], lwd=2)
lines(RotaVirusBB$time, piDirMulti$`00-04`[,2], lwd=2)  
polygon(c(RotaVirusBB$time,rev(RotaVirusBB$time)),c(piDirMulti$`00-04`[,2],rev(piDirMulti$`00-04`[,1])), col=adjustcolor("grey20", alpha.f = 0.15))

axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
title("Predicitive intervals")
legend("bottomleft", legend=c("Dirichelt-multinomial predicitve interval","Multinomial predicitve interval"), col=c(adjustcolor("grey20", alpha.f = 0.15),adjustcolor("grey20", alpha.f = 0.3)), pch=c(22,22), pt.bg=c(adjustcolor("grey20", alpha.f = 0.15),adjustcolor("grey20", alpha.f = 0.3)),pt.cex=1)


## ----Predictiveinterval70, echo=FALSE, fig.cap="95 percent Predictive intervals  for multinomial and Dirichlet-multinomial, age category ''70+''",fig.pos="H",fig.lp="fig:", fig.height=4.9----
# Plot prop and fit of age-category 70+
plot(RotaVirusBB$time, pred_Multinom[,5], type="l", ylim=c(0,max(piDirMulti$`70+`[,2])), lwd=1, lty=1,xlab="Time (months)", ylab="")
lines(RotaVirusBB$time, pred_dirich[,5], type="l", ylim=c(0,max(RotaVirusBB[,"70+"]/RotaVirusBB[,"totalAGs"])), lwd=1, lty=3,xlab="Time (months)", ylab="")
points(RotaVirusBB$time,RotaVirusBB[,"70+"]/RotaVirusBB[,"totalAGs"],lwd=1,col=grey)
# Multinom predInt
lines(RotaVirusBB$time, piMulti$`70+`[,1], lwd=2)
lines(RotaVirusBB$time, piMulti$`70+`[,2], lwd=2)
polygon(c(RotaVirusBB$time,rev(RotaVirusBB$time)),c(piMulti$`70+`[,2],rev(piMulti$`70+`[,1])), col=adjustcolor("grey20", alpha.f = 0.3))
# Dirichlet predInt
lines(RotaVirusBB$time, piDirMulti$`70+`[,1], lwd=2)
lines(RotaVirusBB$time, piDirMulti$`70+`[,2], lwd=2)  
polygon(c(RotaVirusBB$time,rev(RotaVirusBB$time)),c(piDirMulti$`70+`[,2],rev(piDirMulti$`70+`[,1])), col=adjustcolor("grey20", alpha.f = 0.15))
#esthetics and legend
axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
title("Predicitive intervals")
legend("topleft", legend=c("Dirichelt-multinomial predicitve interval","Multinomial predicitve interval"), col=c(adjustcolor("grey20", alpha.f = 0.15),adjustcolor("grey20", alpha.f = 0.3)), pch=c(22,22), pt.bg=c(adjustcolor("grey20", alpha.f = 0.15),adjustcolor("grey20", alpha.f = 0.3)),pt.cex=1)


