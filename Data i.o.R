
RotaVirusBB <- read.csv("RotavirusBBmonth.csv", skip=1)

RotaVirusBB$time <- as.Date(with(RotaVirusBB,paste(Jahr,"-",sprintf("%0.2d",Monat),"-01",sep="")))
colnames(RotaVirusBB)[1:2] <- c("Year", "Month") 
colnames(RotaVirusBB) <- gsub("\\.\\.","-",gsub("^70$","70+",gsub("^X","",colnames(RotaVirusBB))))

agNamesJung   <- c("00-00", "01-01", "02-02", "03-03","04-04")
agNamesMiddle <- c("15-19", "20-24", "25-29","30-39", "40-49", "50-59", "60-69")
RotaVirusBB$"00-04" <- apply(RotaVirusBB[,agNamesJung],1,sum)
RotaVirusBB$"15-69" <- apply(RotaVirusBB[,agNamesMiddle],1,sum)
agNames <- c("00-04","05-09", "10-14","15-69","70+")

RotaVirusBB$totalAGs <- apply(RotaVirusBB[,agNames],1,sum)


RotaVirusBB$t<- seq(from=1,length=length(RotaVirusBB$time))
RotaVirusBB$sin <- with(RotaVirusBB,sin(2*pi/12*t))
RotaVirusBB$cos <- with(RotaVirusBB,cos(2*pi/12*t))

Begin.Date <- (seq(as.Date("01/09/2002", format = "%d/%m/%Y"),by = "year", 
                   length = length(unique(RotaVirusBB[,1]))-1 ))
End.Date <- (seq(as.Date("01/08/2003", format = "%d/%m/%Y"),by = "year", 
                 length = (length(unique(RotaVirusBB[,1]))-1 )))

Months <-c("Sep","Okt","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July","Aug")
years <- 2001:2013
Periods <- paste(head(years,n=-1),"-",tail(years, n=-1),sep="")

TotalAGs <- matrix(,nrow=12,ncol=11)
for (i in 1:11){
  TotalAGs[,i] <- subset(RotaVirusBB$totalAGs, Begin.Date[i]<=RotaVirusBB$time & RotaVirusBB$time<=End.Date[i],select=RotaVirusBB$TotalAGs)
} 

TotalAGs <- data.frame(TotalAGs, x=factor(c(1:12),labels=Months) ,               
                       ymin = apply(TotalAGs,1,FUN=min) ,
                       ymax = apply(TotalAGs,1,FUN=max) ,
                       med = apply(TotalAGs,1,FUN=median)) 
