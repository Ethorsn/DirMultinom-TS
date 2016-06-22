## ----knitr_option, cache=FALSE, echo=FALSE, results='hide'---------------
library(knitr)

## set global chunk options
opts_chunk$set(echo=FALSE, fig.align='center',fig.show="asis", 
               out.width="0.95\\textwidth", fig.height=6,
               par=TRUE, tidy=TRUE, tidy.opts=list(width.cutoff=56), warning=FALSE, error=FALSE)

## ----R_lab,warning=FALSE-------------------------------------------------
#####
# Libraries 
#####
library("RColorBrewer")
library(plotrix)
library(graphics)
library(ggplot2)
library(MGLM)

## ----Data----------------------------------------------------------------
source('Data i.o.R')

## ----TotalcountPlot, fig.cap="Number of reported cases over time, vertical arrow indicates the recommendation of the vaccine by the state of Brandenburg", fig.pos="H",fig.env='figure', fig.height=5----
#####
# Plot, total counts over time 
#####
 v<- which(RotaVirusBB$time=="2009-01-01")
plot(RotaVirusBB$time,RotaVirusBB$totalAGs, type="l",col="grey60",
     xlab="Time (months)",ylab="Reported cases",lwd=2)
   axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
arrows(x0=RotaVirusBB$time[v],y0=200,y1=0,col="black",lty=1, lwd=2)
title("Absolute values of reported cases over time")

## ----MonthPlots, fig.cap="Subsequent periods of September to August. In our dataset, there is great seasonal variation." ,fig.pos="H",fig.lp="fig:", fig.height=7----
col <- c(rep_len(c("grey30", "black", "grey40"), length.out=11))
Type <- c(rep_len(c(1, 2, 3, 4), 11))
par(mar = c(7.3, 4, 4, 2),new=FALSE)
plotPeriod <- function() {
    for (i in 1:(length(unique(RotaVirusBB[,1]))-1)) {
      fun <- if (i==1) plot else lines
        fun(1:12, subset(RotaVirusBB$totalAGs,
         Begin.Date[i]<=RotaVirusBB$time & RotaVirusBB$time<=End.Date[i],
          select=RotaVirusBB$TotalAGs),
          lty=Type[i], type="l",
          col=col[i],lwd=2, ylim=c(0,max(RotaVirusBB$totalAGs)),
          ylab="No. reported cases, each period",xlab="",xaxt="n")
  }
  title("Numbers of reported rotavirus cases")
  axis(1,tck=-0.01, label=FALSE,at=1:12)
  staxlab(1,1:12 ,nlines=12, labels=Months, srt=45)
}
plotPeriod()
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", legend=Periods ,col=col, lty=Type, lwd=2, bg="white", title="Periods",ncol=4, xpd=TRUE, )

## ----MaxMin, fig.cap="Min max and meadian aggregated over months" ,fig.pos="H",fig.lp="fig:",fig.height=5----

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2), new = FALSE)
plotTot <- function() {
  plot(1:12, TotalAGs$med, lty=1, type="l",col="black",lwd=2,
       ylim=c(0,max(RotaVirusBB$totalAGs)) ,
       ylab="No. reported cases, each month",xlab="",xaxt="n")
  lines(1:12, TotalAGs$ymin, lty=4, type="l", lwd=1)
  lines(1:12, TotalAGs$ymax, lty=2, type="l", lwd=1)
  title("Max, median and minimal value in each month")
  axis(1,tck=-0.01, label=FALSE,at=1:12)
  staxlab(1,1:12 ,nlines=12, labels=Months, srt=45)
  legend("left", legend=c("Max","Median", "Min"), lty=c(2,1,4), lwd=2, bg="grey96", title="Lines",ncol=1, xpd=TRUE, border="white")
}
plotTot()

## ----BoxPlots, fig.cap="Boxplots of age-categories in data. (left) all age-categories (right) reduced amount of age-categories \\ \textit{Note, age-categories do not have the same width.}", fig.pos='H', fig.lp="fig:", fig.height=5----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2), new = FALSE)
par(mfrow=c(1,2))
boxplot(RotaVirusBB[,3:17], las=2, ylab="No. reported cases",ymax=c(0,800))
title("All age-categories")
boxplot(RotaVirusBB[,agNames], las=2, ylab="No. reported cases")
title("Reduced age-categories")

## ----propPlots, fig.lp="fig:", fig.cap="Proportions of reduced age-categories. Notice the age-shift betweeen age-category ''00-04'' and ''70+''.", fig.pos='H', fig.height=6----
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(5, 4, 4, 2), new = FALSE)
pal <- c("black","grey20","grey20","grey20","grey10")
Lwd <- c(2,1,1,1,2)
LTY <- c(1,2,2,3,6)
plotProp <- function() {
  for (i in 1:length(agNames)) {
    fun <- if (i==1) plot else lines
    fun(RotaVirusBB$time,RotaVirusBB[,agNames[i]]/RotaVirusBB[,"totalAGs"], type="l",xlab="Time (months)",ylab="Proportion of reported cases", ylim=c(0,max(RotaVirusBB[,agNames]/RotaVirusBB[,"totalAGs"])),col=pal[i],lwd=Lwd[i], lty=LTY[i])
  }
  #Add legend
  axis(1,at=as.numeric(RotaVirusBB$time),label=NA,tck=-0.01)
  legend(x="left",agNames,col=pal,lty=c(1,2,3,4,5),lwd=Lwd,bg="white")
  title("Proportions of reduced age-categories change over time")
}

plotProp()

